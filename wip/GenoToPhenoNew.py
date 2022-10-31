#!/usr/bin/python3
#================================================================
# GenoToPheno.py - extract genotypes and predict phenotypes
#================================================================
# DESCRIPTION
#   This script will take a canine VCF, either single sample or
#   multi-sample, and extract genetic variants of interest,
#   using the consensus genotypes to predict phenotypes for
#   well-studied, simple genetic traits.
#================================================================
# Initialization
#================================================================
import sys, argparse
import os
import numpy as np
import scipy
import scipy.stats
import pandas as pd
import allel
import math
import re
import gzip
import json
workdir=os.getcwd()

arrLen = np.vectorize(len)

parser = argparse.ArgumentParser(description="This script extracts and interprets genotype calls from a canine VCF into phenotype predictions for well-studied physical traits.")

parser.add_argument("-F", "--vcf",
                    help="input variant call file",
                    default="input.vcf.gz")
parser.add_argument("-V", "--variants",
                    help="input variants of interest, comma-delimited file",
                    default="bin/DogAgingProject_VariantsOfInterest.csv")
parser.add_argument("-I", "--imputation",
                    help="If imputed, then imputation panel version? OPTS: impute-v1, impute-v2",
                    default="impute-v2")
parser.add_argument("-R", "--reference",
                    help="Which reference genome? OPTS: CanFam3.1, UU_Cfam_GSD_1.0, Dog10K_Boxer_Tasha",
                    default="CanFam3.1")
parser.add_argument("-O", "--output",
                    help="output filepath and prefix",
                    default="output")

args = parser.parse_args()

if args.vcf:
    print("Variant call file (VCF): %s" % args.vcf)

if args.variants:
    print("Variants of interest: %s" % args.variants)

if args.imputation:
    print("Asserting that imputation panel is %s" % args.imputation)

if args.reference:
    print("Asserting that reference panel is %s" % args.reference)

#================================================================
# Extract Genotypes
#================================================================
def extractGenotypes(vcf,variantTable):
    "loops through variant call file to extract information for variants of interest"

    # initialize lists
    varVAR = []
    varUSE = []
    varSMP = []
    varCHR = []
    varPOS = []
    varTAB = []
    varTRAIT = []
    varLOCUS = []
    varVARIANT = []
    varMODULE = []
    varAlleles = []
    varGT = []
    varGP = []
    varSamples = []

    print("Loading variants from VCF (ignore warnings about missing GP field)...")

    for index, row in variantTable.iterrows():
        print("Extracting information for variant:")
        print(row)

        try:
            thisCall = allel.read_vcf(vcf, fields = ['samples','variants/CHROM','variants/POS','variants/REF','variants/ALT','calldata/GT','calldata/GP'], rename_fields = {'variants/CHROM': 'variants/CHR'}, numbers = {'ALT': 10}, region = str(row['CHR']) +":"+ str(row['POS']) +"-"+ str(row['POS']))
        except RuntimeError:
            print("Update the damned index file!")
            sys.exit()

        if thisCall is None: # if output = None, then chromosome ID may contain prefix "chr", so try adding "chr" to region:
            thisCall = allel.read_vcf(vcf, fields = ['samples','variants/CHROM','variants/POS','variants/REF','variants/ALT','calldata/GT','calldata/GP'], rename_fields = {'variants/CHROM': 'variants/CHR'}, numbers = {'ALT': 10}, region = "chr" + str(row['CHR']) +":"+ str(row['POS']) +"-"+ str(row['POS']))

            if thisCall is None: # if output = None still, then skip variant:
                varVAR.append(index)
                varUSE.append(False)
                continue
            else:
                if "chr" in str(thisCall['variants/CHR'][0]): # remove "chr" from CHR and change type to integer
                    thisCall['variants/CHR'] = np.asarray([int(re.sub("chr","",i)) for i in thisCall['variants/CHR']])

        if thisCall['calldata/GP'].dtype != 'float32':
            print("The GP field item is not float32, likely missing. The GP field item will be removed.")
            thisCall.pop('calldata/GP', None)
        elif thisCall['calldata/GP'].shape[0] > 1:
            print("The GP field item is present, and it contains more than one site. Split multiallelic site detected!")

        # SAMPLES
        # 'samples'
        #   samples
        #   ndim = 1
        #   shape = (# dogs)
        thisCall['samples'] = np.array(thisCall['samples'], dtype = 'str')
        #print("Genotype calls exist for the following " + str(thisCall['samples'].shape[0]) + " samples:")
        #print(thisCall['samples'])

        # CHROMOSOME
        # 'variants/CHR'
        #   chromosome
        #   ndim = 1
        #   shape = (# records, )
        thisCall['variants/CHR'] = np.array(thisCall['variants/CHR'], dtype = 'int')
        #print("variants/CHR")
        #print(thisCall['variants/CHR'].shape)
        #print(thisCall['variants/CHR'].ndim)

        # POSITION (1-based)
        # 'variants/POS'
        #   position
        #   ndim = 1
        #   shape = (# records, )
        thisCall['variants/POS'] = np.array(thisCall['variants/POS'], dtype = 'int')
        #print("variants/POS")
        #print(thisCall['variants/POS'].shape)
        #print(thisCall['variants/POS'].ndim)

        # REFERENCE ALLELE
        # 'variants/REF'
        #   reference allele
        #   ndim = 1
        #   shape = (# records, )
        thisCall['variants/REF'] = np.array(thisCall['variants/REF'], dtype = 'str')
        #print("variants/REF")
        #print(thisCall['variants/REF'].shape)
        #print(thisCall['variants/REF'].ndim)

        # ALTERNATE ALLELE(S)
        # 'variants/ALT'
        #   alternate allele(s)
        #   ndim = 2 = (# records, # alleles)
        thisCall['variants/ALT'] = np.array(thisCall['variants/ALT'], dtype = 'str')
        #print("variants/ALT")
        #print(thisCall['variants/ALT'].shape)
        #print(thisCall['variants/ALT'].ndim)

        # be aware for split multiallelic sites
        # need to get consensus
        # 18	20443726   A   AGCGCGCAGGGAGGCGCGCACCGCTCCGC   0/.
        # 18	20443726   A   AGCGCGC 0/1
        # 18	20443726   A   AGCGCGCAGGG 0/.
        # would be A/AGCGCGC

        # ALLELES
        # 'alleles'
        #   all alleles collapsed per record + "missing" allele ("NA")
        #   ndim = 2
        #   shape = (# records, # alleles)
        thisCall['alleles'] = np.concatenate((np.expand_dims(thisCall['variants/REF'],axis=1),thisCall['variants/ALT'],np.full(shape = (thisCall['variants/REF'].shape[0], 1), fill_value = 'NA')), axis=1)

        # ALL ALLELES
        # 'alleles_all'
        #   all unique alleles collapsed across records, drop empties
        #   ndim = 1
        #   shape = (# alleles, )
        thisCall['alleles_all'] = np.array([x for x in np.unique(np.ndarray.flatten(thisCall['alleles'])) if x], dtype = 'str')
        thisCall['alleles_all'] = thisCall['alleles_all'][thisCall['alleles_all'] != 'NA']

        # GENOTYPE PROBABILITIES
        # 'calldata/GP'
        #   genotype probability
        #   ndim = 3
        #   shape = (# records, # dogs, 3 genotypes)
        if 'calldata/GP' in thisCall.keys():
            thisCall['calldata/GP'] = np.array(thisCall['calldata/GP'])

        # PROBABILITY AVERAGES
        # 'probability_averages'
        #   max genotype probability averaged across records
        #   ndim = 1
        #   shape = (# dogs)
        if 'calldata/GP' in thisCall.keys():
            thisCall['probability_averages'] = np.mean(np.amax(thisCall['calldata/GP'], axis = 2), axis = 0)
        else:
            # set all average GP to 1, same shape of # dogs
            thisCall['probability_averages'] = np.full(thisCall['samples'].shape,1)

        # GENOTYPE CALLS
        # 'calldata/GT'
        #   genotypes
        #   ndim = 3
        #   shape = (# records, # dogs, 2 chromosomes)
        thisCall['calldata/GT'] = np.array(thisCall['calldata/GT'])

        # GENOTYPE ALLELES
        # 'genotypes_alleles'
        #   genotypes in the form of allele strings
        #   ndim = 3
        #   shape = (# dogs, # records, 2 chromosomes)
        thisCall['genotypes_alleles'] = np.array([[thisCall['alleles'][varInd,thisCall['calldata/GT'][varInd,varDog,:]] for varInd in np.arange(thisCall['alleles'].shape[0])] for varDog in np.arange(thisCall['calldata/GT'].shape[1])], dtype = 'str')

        # RESOLVED Genotype
        # 'genotypes_resolved'
        #   genotypes in the form of allele strings, resolved for multiallelics
        #   ndim = 2
        #   shape = (# dogs, 2 chromosomes)
        gi = []
        for varDog in np.arange(thisCall['genotypes_alleles'].shape[0]):
            thisSet = np.unique(thisCall['genotypes_alleles'][varDog,:,:], axis = 0)
            if thisSet.shape[0] == 1:
                gi.append(thisSet)
            else:
                # old way: removing hom ref records:
                #refArray = np.array([thisCall['variants/REF'][0],thisCall['variants/REF'][0]])
                #gi.append(np.delete(thisSet,np.where((thisSet == refArray).all(axis=1)), axis = 0))
                # new way: keeping the longest of the records:
                lengths = arrLen(thisSet).flatten()
                idx = (-lengths).argsort()[:2]
                gi.append(thisSet.flatten()[idx][np.newaxis,:])

        thisCall['genotypes_resolved'] = np.array(gi)

        varVAR.append(index)
        varUSE.append(True)
        varCHR.append(np.unique(thisCall['variants/CHR'])[0])
        varPOS.append(np.unique(thisCall['variants/POS'])[0])
        varTAB.append(row['tab'])
        varTRAIT.append(row['trait'])
        varLOCUS.append(row['locus'])
        varVARIANT.append(row['variant'])
        varMODULE.append(row['module'])
        varAlleles.append(thisCall['alleles_all'])
        varGT.append(thisCall['genotypes_resolved'])
        varGP.append(thisCall['probability_averages'])
        varSamples.append(thisCall['samples'])

    variantCalls = {'VAR': np.array(varVAR),
                    'USE': np.array(varUSE),
                    'CHR': np.array(varCHR),
                    'POS': np.array(varPOS),
                    'tab': np.array(varTAB),
                    'trait': np.array(varTRAIT),
                    'locus': np.array(varLOCUS),
                    'variant': np.array(varVARIANT),
                    'module': np.array(varMODULE),
                    'alleles': np.array(varAlleles),
                    'samples': np.array(varSamples)[0],
                    'genotypes_alleles': np.array(varGT),
                    'probs': np.array(varGP)}
    return variantCalls

#================================================================
# Prediction Modules
#================================================================
def mod_liver(variantCalls,dogIndex,variantTable):
    "autosomal recessive, 'liver' genetics morph black eumelanin pigment in the coat and skin into a brown shade"
    thisModule = "mod_liver"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "Y"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            recCount = recCount + 1
    if result == "N" and recCount >= 2:
        result = "Y"

    return result

def mod_cocoa(variantCalls,dogIndex,variantTable):
    "autosomal recessive, 'cocoa' genetics morph black eumelanin pigment in the coat and skin into a dark brown shade"
    thisModule = "mod_cocoa"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "Y"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            recCount = recCount + 1
    if result == "N" and recCount >= 2:
        result = "Y"

    return result

def mod_dilute(variantCalls,dogIndex,variantTable):
    "autosomal recessive, 'dilute' genetics morph eumelanin pigment into a lighter shade"
    thisModule = "mod_dilute"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "Y"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            recCount = recCount + 1
    if result == "N" and recCount >= 2:
        result = "Y"

    return result

def mod_dom_black(variantCalls,dogIndex,variantTable):
    "Dominant Black: prevents the expression of pheomelanin in the coat."
    thisModule = "mod_dom_black"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount >= 1:
        result = "Y"
    return result

def mod_rec_red(variantCalls,dogIndex,variantTable):
    "Recessive Red: prevents the expression of eumelanin in the coat."
    thisModule = "mod_rec_red"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount >= 2:
        result = "Y"
    elif recCount == 1:
        result = "C"
    return result

def mod_rec_black(variantCalls,dogIndex,variantTable):
    "Recessive Black: prevents the expression of pheomelanin in the coat."
    thisModule = "mod_rec_black"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount >= 2:
        result = "Y"
    return result

def mod_sable(variantCalls,dogIndex,variantTable):
    "Sable: causes shading of eumelanin pigment among pheomelanin pigment"
    thisModule = "mod_sable"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount >= 1:
        result = "Y"
    return result

def mod_tan_points(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "mod_tan_points"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount > 1:
        result = "Y"
    return result

def mod_mask(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "mod_mask"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount >= 1:
        result = "Y"
    return result

def mod_grizzle(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "mod_grizzle"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount == 2:
        result = "Y"
    elif recCount == 1:
        result = "C"
    return result

def mod_domino(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "mod_domino"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount == 2:
        result = "Y"
    elif recCount == 1:
        result = "C"
    return result

def mod_ticking(variantCalls,dogIndex,variantTable):
    "autosomal recessive; associated marker"
    thisModule = "mod_ticking"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "Y"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            recCount = recCount + 1
    if result == "N" and recCount >= 2:
        result = "Y"

    return result

def mod_brindle(variantCalls,dogIndex,variantTable):
    "Brindle: mosaic expression of eumelanin and pheomelanin in the coat caused by an unstable genetic locus and masked by any dominant black allele."
    thisModule = "mod_brindle"

    result = "N"
    recCount = 0
    counts = []

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            counts.append(1)
        else:
            counts.append(0)
    if np.sum(np.array(counts)) == np.array(counts).shape[0]:
        result = "Y"
    elif recCount >= 2:
        result = "Y"
    return result

def mod_merle(variantCalls,dogIndex,variantTable):
    "add a description here"
    thisModule = "mod_merle"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "double"
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            result = "single"

    return result

def mod_harlequin(variantCalls,dogIndex,variantTable):
    "add a description here"
    thisModule = "mod_harlequin"

    result = "N"

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            result = "Y"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            result = "EL"

    return result

def mod_short_legs(variantCalls,dogIndex,variantTable):
    "Checks marker(s) for retrogene insertion(s) of FGF4 that result in shortened legs."
    thisModule = "mod_short_legs"

    result = "N"

    # insertion start
    for varIndex in variantTable.index[variantTable['locus'] == "fgf4c18"].tolist():
        if np.sum(arrLen(variantCalls['genotypes_alleles'][varIndex,dogIndex,:])) > 2:
            result = "Y"

    # G>A marker
    #if result != "Y":
    #    for varIndex in variantTable.index[variantTable['locus'] == "fgf4c18alt"].tolist():
    #        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) >= 1:
    #            result = "M"
    # not ideal, trouble with GSDs

    # A>G marker from Darwin's Ark GWAS
    if result != "Y":
        for varIndex in variantTable.index[variantTable['locus'] == "fgf4c18gwas"].tolist():
            if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) >= 1:
                result = "M"

    return result

def mod_high_altitude_adaptation(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "AltitudeAdaptation"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) >= 1:
            recCount = recCount + 1

    if recCount == 4:
        result = "Y"
    return result

def mod_shedding_propensity(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "Shedding"

    result = "N"
    recCount = 0

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) >= 1:
            recCount = recCount + 1

    if recCount >= 1:
        result = "Y"

    if result == "N":
        curlyRecCount = 0
        for varIndex in variantTable.index[variantTable['module'] == "mod_curly_coat"].tolist():
            curlyRecCount = curlyRecCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])
        if curlyRecCount >= 1:
            result = "Y"
    return result

#================================================================
# Mulit-Module Trait Predictions
#================================================================
def trait_black(variantCalls,dogIndex,variantTable):
    "The trait reported as 'black' describes the way eumelanin is expressed in the base coat color of a dog. This takes into consideration the modules: mod_liver(), mod_dilute(), mod_dom_black(), mod_rec_black(), mod_rec_red(), mod_cocoa()."

    eumelaninHex = {"black": "#17181D",
                    "dilute_black": "#746F76",
                    "liver": "#442825",
                    "dilute_liver": "#947D80",
                    "cocoa": "#21100e",
                    "dilute_cocoa": "#947D80"}

    eumelaninStr = {"black": "Black",
                    "dilute_black": "Dilute Black",
                    "liver": "Brown",
                    "dilute_liver": "Dilute Brown",
                    "cocoa": "Cocoa",
                    "dilute_cocoa": "Dilute Cocoa"}

    liver = mod_liver(variantCalls,dogIndex,variantTable)
    dilute = mod_dilute(variantCalls,dogIndex,variantTable)
    dom_black = mod_dom_black(variantCalls,dogIndex,variantTable)
    rec_black = mod_rec_black(variantCalls,dogIndex,variantTable)
    rec_red = mod_rec_red(variantCalls,dogIndex,variantTable)
    cocoa = mod_cocoa(variantCalls,dogIndex,variantTable)

    predictPhenoID = "uncalled"
    predictPhenoSTR = "error"
    predictPhenoVAL = ""
    predictPhenoHEX = "#ff0000"
    predictPhenoIMG = ""

    if cocoa == "Y" and dilute == "N":
        phenoID = "cocoa"
        if rec_red == "Y":
            predictPhenoID = "cocoa_leathers"
            predictPhenoSTR = "Cocoa (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "cocoa_solid"
            predictPhenoSTR = "Cocoa (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "cocoa"
            predictPhenoSTR = "Cocoa"
            predictPhenoHEX = eumelaninHex[phenoID]
    elif cocoa == "Y" and dilute == "Y":
        phenoID = "dilute_cocoa"
        if rec_red == "Y":
            predictPhenoID = "dilute_cocoa_leathers"
            predictPhenoSTR = "Dilute Cocoa (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "dilute_cocoa_solid"
            predictPhenoSTR = "Dilute Cocoa (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "dilute_cocoa"
            predictPhenoSTR = "Dilute Cocoa"
            predictPhenoHEX = eumelaninHex[phenoID]
    elif liver == "N" and dilute == "N":
        phenoID = "black"
        if rec_red == "Y":
            predictPhenoID = "black_leathers"
            predictPhenoSTR = "Black (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "black_solid"
            predictPhenoSTR = "Black (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "black"
            predictPhenoSTR = "Black"
            predictPhenoHEX = eumelaninHex[phenoID]
    elif liver == "N" and dilute == "Y":
        phenoID = "dilute_black"
        if rec_red == "Y":
            predictPhenoID = "dilute_black_leathers"
            predictPhenoSTR = "Dilute Black (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "dilute_black_solid"
            predictPhenoSTR = "Dilute Black (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "dilute_black"
            predictPhenoSTR = "Dilute Black"
            predictPhenoHEX = eumelaninHex[phenoID]
    elif liver == "Y" and dilute == "N":
        phenoID = "liver"
        if rec_red == "Y":
            predictPhenoID = "liver_leathers"
            predictPhenoSTR = "Brown (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "liver_solid"
            predictPhenoSTR = "Brown (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "liver"
            predictPhenoSTR = "Brown"
            predictPhenoHEX = eumelaninHex[phenoID]
    elif liver == "Y" and dilute == "Y":
        phenoID = "dilute_liver"
        if rec_red == "Y":
            predictPhenoID = "dilute_liver_leathers"
            predictPhenoSTR = "Dilute Brown (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "dilute_liver_solid"
            predictPhenoSTR = "Dilute Brown (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "dilute_liver"
            predictPhenoSTR = "Dilute Brown"
            predictPhenoHEX = eumelaninHex[phenoID]
    else:
        phenoID = "black"
        if rec_red == "Y":
            predictPhenoID = "black_leathers"
            predictPhenoSTR = "Black (nose and paws)"
            predictPhenoHEX = eumelaninHex[phenoID]
        elif dom_black == "Y" or rec_black == "Y":
            predictPhenoID = "black_solid"
            predictPhenoSTR = "Black (solid coat)"
            predictPhenoHEX = eumelaninHex[phenoID]
        else:
            predictPhenoID = "black"
            predictPhenoSTR = "Black"
            predictPhenoHEX = eumelaninHex[phenoID]

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_red(variantCalls,dogIndex,variantTable):
    "The trait reported as 'red' describes the way pheomelanin is expressed in the base coat color of a dog. This takes into consideration the modules: mod_red_intensity()."

    pheomelaninHex = {"red": "#7E341B",
                      "tan": "#A86B39",
                      "yellow": "#CE9E63",
                      "cream": "#D8C4AF"}

    pheomelaninStr = {"red": "Red",
                     "tan": "Tan",
                     "yellow": "Yellow",
                     "cream": "Cream"}

    #red_intensity = mod_red_intensity(variantCalls,dogIndex,variantTable)

    # get the red intensity alleles: i1,i2,i3,i4,i5:
    for varIndex in variantTable.index[variantTable['module'] == "mod_red_intensity"].tolist():
        recCount = np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    modelRedIntensity = pd.DataFrame({"variant": ["i1","i1","i2","i2","i3","i3","i4","i4","i5","i5"],
                                      "genotype": ["het","hom","het","hom","het","hom","het","hom","het","hom"],
                                      "coeff": [0.472,1.068,0.057,0.208,0.234,0.208,0.700,1.232,0.199,0.222]})
    counts = []
    for i in modelRedIntensity['variant'].unique():
        varIndex = variantTable.index[variantTable['variant'] == i][0]
        thisCount = np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])
        if thisCount == 1:
            counts.append(1)
        else:
            counts.append(0)

        if thisCount == 2:
            counts.append(1)
        else:
            counts.append(0)
    modelRedIntensity['value'] = np.array(counts)

    # model takes in 0/1 for if het and if hom:
    redIntensity = -1.504 + np.sum(modelRedIntensity['value']*modelRedIntensity['coeff'])

    # 1.46 is max red
    # -1.503 is min red
    # .75 increments

    if redIntensity < -0.75:
        predictPhenoID = "cream"
        predictPhenoSTR = pheomelaninStr[predictPhenoID]
        predictPhenoVAL = ""
        predictPhenoHEX = pheomelaninHex[predictPhenoID]
        predictPhenoIMG = ""
    elif redIntensity < 0:
        predictPhenoID = "yellow"
        predictPhenoSTR = pheomelaninStr[predictPhenoID]
        predictPhenoVAL = ""
        predictPhenoHEX = pheomelaninHex[predictPhenoID]
        predictPhenoIMG = ""
    elif redIntensity < 0.75:
        predictPhenoID = "tan"
        predictPhenoSTR = pheomelaninStr[predictPhenoID]
        predictPhenoVAL = ""
        predictPhenoHEX = pheomelaninHex[predictPhenoID]
        predictPhenoIMG = ""
    elif redIntensity >= 0.75:
        predictPhenoID = "red"
        predictPhenoSTR = pheomelaninStr[predictPhenoID]
        predictPhenoVAL = ""
        predictPhenoHEX = pheomelaninHex[predictPhenoID]
        predictPhenoIMG = ""

    # until a prediction model is finalized, default to tan:
    #predictPhenoID = "tan"
    #predictPhenoSTR = pheomelaninStr[predictPhenoID]
    #predictPhenoVAL = ""
    #predictPhenoHEX = pheomelaninHex[predictPhenoID]
    #predictPhenoIMG = ""

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_agouti(variantCalls,dogIndex,variantTable):
    "The trait reported as 'agouti' describes the way eumelanin is expressed against pheomelanin in a dog's coat. This takes into consideration modules: mod_sable(), mod_tan_points(), mod_rec_black(), mod_dom_black(), mod_rec_red()."

    sable = mod_sable(variantCalls,dogIndex,variantTable)
    tan_points = mod_tan_points(variantCalls,dogIndex,variantTable)
    rec_black = mod_rec_black(variantCalls,dogIndex,variantTable)
    dom_black = mod_dom_black(variantCalls,dogIndex,variantTable)
    rec_red = mod_rec_red(variantCalls,dogIndex,variantTable)

    # "coat-pattern.svg", "sable.svg", "Sable", "masking.svg", "Masking"

    predictPhenoID = ""
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    if sable == "Y":
        if (rec_red == "Y") or (dom_black == "Y"):
            predictPhenoID = "sable_hidden"
            predictPhenoSTR = "Sable (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_sable.svg"
        else:
            predictPhenoID = "sable"
            predictPhenoSTR = "Sable"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_sable.svg"
    elif tan_points == "Y":
        if (rec_red == "Y") or (dom_black == "Y"):
            predictPhenoID = "tan_points_hidden"
            predictPhenoSTR = "Tan Points (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_tanPoints.svg"
        else:
            predictPhenoID = "tan_points"
            predictPhenoSTR = "Tan Points"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_tanPoints.svg"
    elif rec_black == "Y":
        if (rec_red == "Y") or (dom_black == "Y"):
            predictPhenoID = "rec_black_hidden"
            predictPhenoSTR = "No sable, agouti, or tan points"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = ""
        else:
            predictPhenoID = "rec_black"
            predictPhenoSTR = "No sable, agouti, or tan points"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = ""
    else:
        if (rec_red == "Y") or (dom_black == "Y"):
            predictPhenoID = "agouti_hidden"
            predictPhenoSTR = "Agouti (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_agouti.svg"
        else:
            predictPhenoID = "agouti"
            predictPhenoSTR = "Agouti"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusA_agouti.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_extension(variantCalls,dogIndex,variantTable):
    "The trait reported as 'extension' describes the way eumelanin is expressed against pheomelanin in a dog's coat. This takes into consideration modules: mod_mask(), mod_grizzle(), mod_domino(), mod_rec_black(), mod_dom_black(), mod_rec_red()."

    mask = mod_mask(variantCalls,dogIndex,variantTable)
    grizzle = mod_grizzle(variantCalls,dogIndex,variantTable)
    domino = mod_domino(variantCalls,dogIndex,variantTable)
    rec_black = mod_rec_black(variantCalls,dogIndex,variantTable)
    dom_black = mod_dom_black(variantCalls,dogIndex,variantTable)
    rec_red = mod_rec_red(variantCalls,dogIndex,variantTable)

    predictPhenoID = "normal_extension"
    predictPhenoSTR = "No mask, grizzle, or domino patterns"
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    mask_n = 0
    for varIndex in variantTable.index[variantTable['module'] == "mod_mask"].tolist():
        mask_n = mask_n + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    grizzle_n = 0
    for varIndex in variantTable.index[variantTable['module'] == "mod_grizzle"].tolist():
        grizzle_n = grizzle_n + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    domino_n = 0
    for varIndex in variantTable.index[variantTable['module'] == "mod_domino"].tolist():
        domino_n = domino_n + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    rec_red_n = 0
    for varIndex in variantTable.index[variantTable['module'] == "mod_rec_red"].tolist():
        rec_red_n = rec_red_n + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if mask_n >= 1:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "mask_hidden"
            predictPhenoSTR = "Facial Mask (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_mask.svg"
        else:
            predictPhenoID = "mask"
            predictPhenoSTR = "Facial Mask"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_mask.svg"
    elif grizzle_n == 2:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "grizzle_hidden"
            predictPhenoSTR = "Sighthound Grizzle (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
        else:
            predictPhenoID = "grizzle"
            predictPhenoSTR = "Sighthound Grizzle"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
    elif domino_n == 2:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "domino_hidden"
            predictPhenoSTR = "Northern Domino (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_domino.svg"
        else:
            predictPhenoID = "domino"
            predictPhenoSTR = "Northern Domino"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    elif grizzle_n == 1 and domino_n == 1:
        predictPheoID = "grizzle_domino"
        predictPhenoSTR = "Grizzle or Domino"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    elif rec_red_n == 2:
        predictPhenoID = "rec_red"
        predictPhenoSTR = "No mask, grizzle, or domino patterns"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""
    elif rec_red_n == 1 and grizzle_n == 1:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "grizzle_hidden"
            predictPhenoSTR = "Sighthound Grizzle (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
        else:
            predictPhenoID = "grizzle"
            predictPhenoSTR = "Sighthound Grizzle"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
    elif rec_red_n == 1 and domino_n == 1:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "domino_hidden"
            predictPhenoSTR = "Northern Domino (hidden)"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_domino.svg"
        else:
            predictPhenoID = "domino"
            predictPhenoSTR = "Northern Domino"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    else:
        if (dom_black == "Y") or (rec_black == "Y"):
            predictPhenoID = "normal_extension_hidden"
            predictPhenoSTR = "No mask, grizzle, or domino patterns"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = ""
        else:
            predictPhenoID = "normal_extension"
            predictPhenoSTR = "No mask, grizzle, or domino patterns"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = ""

    # if mask == "Y":
    #     if (dom_black == "Y") or (rec_black == "Y"):
    #         predictPhenoID = "mask_hidden"
    #         predictPhenoSTR = "Facial Mask (hidden)"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_mask.svg"
    #     else:
    #         predictPhenoID = "mask"
    #         predictPhenoSTR = "Facial Mask"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_mask.svg"
    # elif grizzle == "Y":
    #     if (dom_black == "Y") or (rec_black == "Y"):
    #         predictPhenoID = "grizzle_hidden"
    #         predictPhenoSTR = "Sighthound Grizzle (hidden)"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
    #     else:
    #         predictPhenoID = "grizzle"
    #         predictPhenoSTR = "Sighthound Grizzle"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_grizzle.svg"
    # elif domino == "Y":
    #     if (dom_black == "Y") or (rec_black == "Y"):
    #         predictPhenoID = "domino_hidden"
    #         predictPhenoSTR = "Northern Domino (hidden)"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    #     else:
    #         predictPhenoID = "domino"
    #         predictPhenoSTR = "Northern Domino"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    # elif (grizzle == "C" and domino == "C") or (rec_red == "C" and grizzle == "C") or (rec_red == "C" and domino == "C"):
    #     predictPheoID = "grizzle_domino"
    #     predictPhenoSTR = "Grizzle or Domino"
    #     predictPhenoHEX = ""
    #     predictPhenoVAL = ""
    #     predictPhenoIMG = "CoatPattern_locusE_domino.svg"
    # elif rec_red == "Y":
    #     predictPhenoID = "rec_red"
    #     predictPhenoSTR = "No mask, grizzle, or domino patterns"
    #     predictPhenoHEX = ""
    #     predictPhenoVAL = ""
    #     predictPhenoIMG = ""
    # else:
    #     if (dom_black == "Y") or (rec_black == "Y"):
    #         predictPhenoID = "normal_extension_hidden"
    #         predictPhenoSTR = "No mask, grizzle, or domino patterns"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = ""
    #     else:
    #         predictPhenoID = "normal_extension"
    #         predictPhenoSTR = "No mask, grizzle, or domino patterns"
    #         predictPhenoHEX = ""
    #         predictPhenoVAL = ""
    #         predictPhenoIMG = ""

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_ticking(variantCalls,dogIndex,variantTable):
    "The trait reported as 'ticking' describes flecks of color that appear among otherwise white fur, sometimes in a dense, mottled manner, known as roaning. It takes into consideration the module mod_ticking()."

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    ticking = mod_ticking(variantCalls,dogIndex,variantTable)

    if ticking == "Y":
        predictPhenoID = "ticking"
        predictPhenoSTR = "Ticking or Roaning"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_ticking.svg"
    elif ticking == "N":
        predictPhenoID = "no_ticking"
        predictPhenoSTR = ""
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_brindle(variantCalls,dogIndex,variantTable):
    "The trait reported as 'brindle' describes a pattern of coat color alternating in eumelanin and pheomelanin. It takes into consideration the modules mod_brindle() and mod_dominant_black()."

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    brindle = mod_brindle(variantCalls,dogIndex,variantTable)
    recRed = mod_rec_red(variantCalls,dogIndex,variantTable)
    domBlack = mod_dom_black(variantCalls,dogIndex,variantTable)

    if brindle == "Y" and domBlack != "Y" and recRed != "Y":
        predictPhenoID = "brindle"
        predictPhenoSTR = "Brindle"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_brindle.svg"
    elif brindle == "Y" and recRed == "Y":
        predictPhenoID = "brindle_hidden"
        predictPhenoSTR = "Brindle (hidden)"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_brindle.svg"
    elif brindle == "Y" and domBlack == "Y":
        predictPhenoID = "not_brindle"
        predictPhenoSTR = ""
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""
    elif brindle == "N":
        predictPhenoID = "not_brindle"
        predictPhenoSTR = ""
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_merle(variantCalls,dogIndex,variantTable):
    "add desc"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    merle = mod_merle(variantCalls,dogIndex,variantTable)
    harlequin = mod_harlequin(variantCalls,dogIndex,variantTable)

    if merle == "double" and harlequin == "N":
        predictPhenoID = "double_merle"
        predictPhenoSTR = "Double Merle"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_merle_double.svg"
    elif merle == "single" and harlequin == "N":
        predictPhenoID = "merle"
        predictPhenoSTR = "Merle"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_merle.svg"
    elif merle == "double" and harlequin == "Y":
        predictPhenoID = "harlequin_double_merle"
        predictPhenoSTR = "Harlequin-type Merle (Double Merle carrier)"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_harlequin.svg"
    elif merle == "single" and harlequin == "Y":
        predictPhenoID = "harlequin"
        predictPhenoSTR = "Harlequin-type Merle"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatPattern_harlequin.svg"
    elif merle == "N" and harlequin == "N":
        predictPhenoID = "not_merle"
        predictPhenoSTR = ""
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""
    elif harlequin == "EL":
        predictPhenoID = "embryonic_lethal"
        predictPhenoSTR = ""
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = ""

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_shedding_propensity(variantCalls,dogIndex,variantTable):
    "add desc"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    shedding_propensity = mod_shedding_propensity(variantCalls,dogIndex,variantTable)

    if shedding_propensity == "Y":
        predictPhenoID = "low_shedding"
        predictPhenoSTR = "Low shedding"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatType_Shedding_low.svg"
    elif shedding_propensity == "N":
        predictPhenoID = "normal_shedding"
        predictPhenoSTR = "Normal shedding"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatType_Shedding_normal.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_leg_length(variantCalls,dogIndex,variantTable):
    "The trait reported as 'leg length' describes shortened legs. It takes into consideration the modules mod_short_legs() and mod_long_legs()"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    short_legs = mod_short_legs(variantCalls,dogIndex,variantTable)

    if short_legs == "Y":
        predictPhenoID = "short_legs"
        predictPhenoSTR = "Shortened leg length"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "SpecialFeatures_Limbs_short.svg"
    elif short_legs == "M":
        predictPhenoID = "short_legs_marker"
        predictPhenoSTR = "Shortened leg length"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "SpecialFeatures_Limbs_short.svg"
    else:
        predictPhenoID = "normal_legs"
        predictPhenoSTR = "Normal leg length"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "SpecialFeatures_Limbs_normal.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_high_altitude_adaptation(variantCalls,dogIndex,variantTable):
    "add desc"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    high_altitude_adaptation = mod_high_altitude_adaptation(variantCalls,dogIndex,variantTable)

    if high_altitude_adaptation == "Y":
        predictPhenoID = "hypoxia_adapted"
        predictPhenoSTR = "Adaptable to high altitudes"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "SpecialFeatures_Altitudes_adapted.svg"
    elif high_altitude_adaptation == "N":
        predictPhenoID = "not_hypoxia_adapted"
        predictPhenoSTR = "No adaptation to high altitudes"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "SpecialFeatures_Altitudes_not_adapted.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_coat_length(variantCalls,dogIndex,variantTable):
    "The trait reported as 'coat length' describes how long the dog's fur grows. It takes into consideration the module 'mod_coat_length'."
    thisModule = "mod_coat_length"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            predictPhenoID = "long_coat"
            predictPhenoSTR = "Long coat"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatType_Length_long.svg"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) <= 1:
            predictPhenoID = "short_coat"
            predictPhenoSTR = "Short coat"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatType_Length_short.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_coat_texture(variantCalls,dogIndex,variantTable):
    "The trait reported as 'coat texture' describes whether a dog's coat is straight, wavy, or curly. It considers the module 'mod_coat_length'."
    thisModule = "mod_curly_coat"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    recCount = 0
    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        recCount = recCount + np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele'])

    if recCount == 2:
        predictPhenoID = "curly_coat"
        predictPhenoSTR = "Curly coat"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatType_Curl_curly.svg"
    elif recCount == 1:
        predictPhenoID = "wavy_coat"
        predictPhenoSTR = "Wavy coat"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatType_Curl_wavy.svg"
    elif recCount == 0:
        predictPhenoID = "straight_coat"
        predictPhenoSTR = "Straight coat"
        predictPhenoHEX = ""
        predictPhenoVAL = ""
        predictPhenoIMG = "CoatType_Curl_straight.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

def trait_furnishings(variantCalls,dogIndex,variantTable):
    "The trait reported as 'furnishings' describes whether a dog has long hair on its brows and muzzle. It considers the module 'mod_furnishings'."
    thisModule = "mod_furnishings"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) >= 1:
            predictPhenoID = "furnishings"
            predictPhenoSTR = "Eyebrow and muzzle furnishings"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatType_Furnishings_furnished.svg"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 0:
            predictPhenoID = "no_furnishings"
            predictPhenoSTR = "No eyebrow and muzzle furnishings"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "CoatType_Furnishings_unfurnished.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG


def trait_natural_bob_tail(variantCalls,dogIndex,variantTable):
    "add desc"
    thisModule = "BobTail"

    predictPhenoID = "uncalled"
    predictPhenoSTR = ""
    predictPhenoHEX = ""
    predictPhenoVAL = ""
    predictPhenoIMG = ""

    for varIndex in variantTable.index[variantTable['module'] == thisModule].tolist():
        if np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 1:
            predictPhenoID = "natural_bob_tail"
            predictPhenoSTR = "Natural bob tail"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "SpecialFeatures_Tail_bob_tail.svg"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 2:
            predictPhenoID = "embryonic_lethal"
            predictPhenoSTR = "Natural bob tail"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "SpecialFeatures_Tail_bob_tail.svg"
        elif np.sum(variantCalls['genotypes_alleles'][varIndex,dogIndex,:] == variantTable.iloc[varIndex]['variantAllele']) == 0:
            predictPhenoID = "normal_tail"
            predictPhenoSTR = "Normal length tail"
            predictPhenoHEX = ""
            predictPhenoVAL = ""
            predictPhenoIMG = "SpecialFeatures_Tail_normal_tail.svg"

    return predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG

#def BobTail(genetics):
#    if genetics[36] + genetics[37] == 0:
#        tail = "full length tail"
#    elif genetics[36] + genetics[37] == 1:
#        tail = "natural bobtail"
#    else: tail = "lethal genotype"
#    return tail

def trailblazerFlair(variantCalls,dogIndex,variantTable):
    "Generates the text about coat color and pattern predictions for Darwin's Ark Trailblazer booklets."

    return flair

#================================================================
# Generate Genotypes
#================================================================
def generateGenotypes(variantCalls):
    "loops through dogs and generates genotypes"

    gli = []
    for dogIndex, dog in enumerate(variantCalls['samples']):
        genotypeTableDog = pd.DataFrame(data = {'sample': dog,
                                            'CHR': variantCalls['CHR'],
                                            'POS': variantCalls['POS'],
                                            'tab': variantCalls['tab'],
                                            'trait': variantCalls['trait'],
                                            'locus': variantCalls['locus'],
                                            'variant': variantCalls['variant'],
                                            'A1': variantCalls['genotypes_alleles'][:,dogIndex,0,0],
                                            'A2': variantCalls['genotypes_alleles'][:,dogIndex,0,1],
                                            'CONF': variantCalls['probs'][:,dogIndex]})
        gli.append(genotypeTableDog)
    genotypeTable = pd.concat(gli, axis=0, ignore_index=True)

    return genotypeTable

#================================================================
# Predict Phenotypes
#================================================================
def predictPhenotypes(variantCalls,variantTable):
    "loops through dogs and predicts phenotypes"

    pli = []
    gli = []
    for dogIndex, dog in enumerate(variantCalls['samples']):
        print("Processing phenotypes for dog number " + str(dogIndex) + " with sample ID: " + str(dog))

        tabDog = []
        traitDog = []
        predictPhenoIDDog = []
        predictPhenoSTRDog = []
        predictPhenoVALDog = []
        predictPhenoHEXDog = []
        predictPhenoIMGDog = []

        # Coat Color #

        # trait: black
        tabDog.append("CoatColor")
        traitDog.append("black")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_black(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: red
        tabDog.append("CoatColor")
        traitDog.append("red")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_red(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # Coat Pattern #

        # trait: agouti
        tabDog.append("CoatPattern")
        traitDog.append("agouti")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_agouti(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: extension
        tabDog.append("CoatPattern")
        traitDog.append("extension")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_extension(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: ticking
        tabDog.append("CoatPattern")
        traitDog.append("ticking")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_ticking(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: brindle
        tabDog.append("CoatPattern")
        traitDog.append("brindle")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_brindle(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: merle
        tabDog.append("CoatPattern")
        traitDog.append("merle")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_merle(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # Coat Type #

        # trait: coat_texture
        tabDog.append("CoatType")
        traitDog.append("Coat Texture")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_coat_texture(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: coat_length
        tabDog.append("CoatType")
        traitDog.append("Coat Length")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_coat_length(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: furnishings
        tabDog.append("CoatType")
        traitDog.append("Coat Furnishings")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_furnishings(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: shedding_propensity
        tabDog.append("CoatType")
        traitDog.append("Shedding Propensity")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_shedding_propensity(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # Special Features #

        # trait: natural_bob_tail
        tabDog.append("SpecialFeatures")
        traitDog.append("Skeletal - Tail Length")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_natural_bob_tail(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: leg_length
        tabDog.append("SpecialFeatures")
        traitDog.append("Skeletal - Leg Length")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_leg_length(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        # trait: high_altitude_adaptation
        tabDog.append("SpecialFeatures")
        traitDog.append("High Altitude Adaptation")
        predictPhenoID,predictPhenoSTR,predictPhenoVAL,predictPhenoHEX,predictPhenoIMG = trait_high_altitude_adaptation(variantCalls,dogIndex,variantTable)
        predictPhenoIDDog.append(predictPhenoID)
        predictPhenoSTRDog.append(predictPhenoSTR)
        predictPhenoVALDog.append(predictPhenoVAL)
        predictPhenoHEXDog.append(predictPhenoHEX)
        predictPhenoIMGDog.append(predictPhenoIMG)

        phenotypeTableDog = pd.DataFrame(data = {'sample': dog,
                                                 'tab': np.array(tabDog),
                                                 'trait': np.array(traitDog),
                                                 'result': np.array(predictPhenoIDDog),
                                                 'string': np.array(predictPhenoSTRDog),
                                                 'value': np.array(predictPhenoVALDog),
                                                 'color': np.array(predictPhenoHEXDog),
                                                 'image': np.array(predictPhenoIMGDog)})

        pli.append(phenotypeTableDog)
    phenotypeTable = pd.concat(pli, axis=0, ignore_index=True)

    return phenotypeTable

#================================================================
# Main Function
#================================================================
def main(args):
    "Extracts genotypes for variants of interest and calls interpretation modules to output phenotype predictions."

    # load variants of interest
    variantTable = pd.read_csv(args.variants, header = 0, na_values = "NA")

    # select imputation and reference
    variantTable.drop(variantTable[variantTable['imputation'] != args.imputation].index, inplace = True)
    variantTable.drop(variantTable[variantTable['reference'] != args.reference].index, inplace = True)
    variantTable.reset_index(inplace = True)

    # set CHR and POS to integers
    variantTable['CHR'] = variantTable['CHR'].astype(int)
    variantTable['POS'] = variantTable['POS'].astype(int)

    # load data from variant call file for variants of interest
    variantCalls = extractGenotypes(args.vcf, variantTable)

    # drop unobserved variants and re-index variantTable
    variantTable = variantTable[variantCalls['USE']]
    variantTable.reset_index(inplace = True)

    # trailblazers
    # generate ColorWheel

    # generate genotype table
    genotypeTable = generateGenotypes(variantCalls)
    genotypeTable.to_csv(args.output + "_genotypeTable.csv", sep = ",", na_rep = "NA", index = False)

    # generate phenotype table
    phenotypeTable = predictPhenotypes(variantCalls,variantTable)
    phenotypeTable.to_csv(args.output + "_phenotypeTable.csv", sep = ",", na_rep = "NA", index = False)

    # Darwin's Ark: Trailblazers
    trailblazerTable = genotypeTable.merge(variantTable[['tab','trait','locus','variant','gene','title','desc','normalAllele','variantAllele']], on = ['tab','trait','locus','variant'])[['sample','tab','trait','title','gene','normalAllele','variantAllele','A1','A2','desc']].rename(columns = {'sample': 'Sample ID', 'title': 'Trait', 'gene': 'Gene', 'normalAllele': 'Normal Version', 'variantAllele': 'Variant Version', 'A1': 'First Copy', 'A2': 'Second Copy', 'desc':'Description'})
    trailblazerTable.to_csv(args.output + "_trailblazerGenotypeTable.csv", sep = ",", na_rep = "NA", index = False)
    phenotypeTable[['sample','tab','result','string']].to_csv(args.output + "_trailblazerPhenotypeTable.tsv", sep = "\t", na_rep = "NA", index = False)

    # Dog Aging Project: Genomic Reports

    # sample list
    sampleList = pd.unique(genotypeTable['sample'])

    # tabs list
    tabsList = pd.unique(variantTable['tab'])

    # generate json per dog per tab
    jsonTable = genotypeTable.merge(variantTable[['tab','trait','locus','variant','gene','title','desc','normalAllele','variantAllele']], on = ['tab','trait','locus','variant'])[['sample','tab','title','gene','normalAllele','variantAllele','A1','A2']].dropna().rename(columns = {'title': 'name', 'A1': 'firstCopy', 'A2': 'secondCopy'})
    jsonTable['possibleAlleles'] = jsonTable['normalAllele'].astype(str) + " & " + jsonTable['variantAllele']
    jsonTable = jsonTable[jsonTable['firstCopy'] != "NA"]
    jsonTable = jsonTable[jsonTable['secondCopy'] != "NA"]
    jsonTable = jsonTable.drop(columns = ['normalAllele','variantAllele'])
    jsonTable['effect'] = 0
    jsonTable.to_csv(args.output + "_jsonTable.csv", sep = ",", na_rep = "null", index = False)
    #for s in sampleList:
    #    for tab in tabList:
    #        jsonTable[jsonTable['sample']==s].drop(columns = ['sample']).to_json("StudyID-{}".format(s) + ".genotypesTable." + tab + ".json", orient='records')

if __name__ == '__main__':
    main(args)

#!/usr/bin/python3
#================================================================
# ReportsToJSON.py - populate JSON for Dog Aging Project
#================================================================
# DESCRIPTION
#   This script will take a template JSON and results for a dog
#   and populate fields for that dog.
#================================================================
# Initialization
#================================================================
import sys, argparse
import os
import numpy as np
import pandas as pd
import json
workdir=os.getcwd()

parser = argparse.ArgumentParser(description="This script will take a template JSON and results for a query dog, and populate the JSON")

parser.add_argument("-J", "--json",
                    help="input template JSON",
                    default="ref/GenomicReport_TemplateJSON.json")

parser.add_argument("-S", "--studyid",
                    help="input study ID",
                    default="3")

parser.add_argument("-A", "--ancestry",
                    help="input ancestry table or json with fields `breed` and `percent`",
                    default="StudyID-3.GlobalAncestry.json")

parser.add_argument("-I", "--inbreeding",
                    help="input coefficient of inbreeding value or file of COIs",
                    default="COI.csv")

parser.add_argument("-C", "--coifile",
                    help="",
                    default="true")

parser.add_argument("-SP", "--sizephenofile",
                    help="input phenotypes table for complex trait prediction of body size",
                    default="BodySize.phenotypes.csv")

parser.add_argument("-SG", "--sizegenofile",
                    help="input genotypes table for complex trait prediction of body size",
                    default="BodySize.genotypes.csv")

parser.add_argument("-WP", "--whitephenofile",
                    help="input phenotypes table for complex trait prediction of white spotting",
                    default="WhiteSpotting.phenotypes.csv")

parser.add_argument("-WG", "--whitegenofile",
                    help="input genotypes table for complex trait prediction of white spotting",
                    default="WhiteSpotting.genotypes.csv")

parser.add_argument("-P", "--phenotypes",
                    help="input phenotypes table for simple trait predictions, collated by `tab` and `trait`",
                    default="phenotypeTable.csv")

parser.add_argument("-G", "--genotypes",
                    help="input genotypes table for simple trait predictions, collated by `tab` and `trait`",
                    default="jsonTable.csv")

parser.add_argument("-H", "--html",
                    help="fill in HTML elements for content? `noHTML` if not",
                    default="noHTML")

args = parser.parse_args()

#================================================================
# LOAD JSON
#================================================================
tmpJSON = json.load(open(args.json))

dogJSON = tmpJSON
dogJSON['id'] = str(args.studyid)

#================================================================
# BREED ANCESTRY
#================================================================
dogADM = json.load(open(args.ancestry))

dogJSON['data'] = dogADM

#================================================================
# COEFFICIENT OF INBREEDING
#================================================================
if args.coifile == "true":
    coifile = pd.read_csv(args.inbreeding)
    dogCOI = coifile[coifile['id'].astype('str') == str(args.studyid)]['coi'].values[0]
else:
    dogCOI = args.inbreeding

dogJSON['inbreeding'] = dogCOI

#================================================================
# TRAIT PREDICTIONS: Predictive Models
#================================================================

# HTML 'content'
if args.html == "noHTML":
    # BodySize
    dogJSON['panels'][0]['top']['content'] = "null"
    dogJSON['panels'][0]['results']['content'] = "null"

    # WhiteSpotting
    dogJSON['panels'][3]['top']['content'] = "null"
    dogJSON['panels'][3]['results']['content'] = "null"
else:
    # BodySize
    dogJSON['panels'][0]['top']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.BodySize.top.content.html").read().replace('\n', "")
    dogJSON['panels'][0]['results']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.BodySize.results.content.html").read().replace('\n', "")

    # WhiteSpotting
    dogJSON['panels'][3]['top']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.WhiteSpotting.top.content.html").read().replace('\n', "")
    dogJSON['panels'][3]['results']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.WhiteSpotting.results.content.html").read().replace('\n', "")

# top: 'genericValues'

# BodySize
phenoBodySize = pd.read_csv(args.sizephenofile)
dogPheno_BodySize = phenoBodySize[phenoBodySize['id'].astype('str') == str(args.studyid)]['prediction'].values[0]

# rescale from [0,4] to [0.25,3.75] on [0,4] scale
oldVal = dogPheno_BodySize
rmin=0.25
rmax=3.75
tmin=0.00
tmax=4.00
newVal = (((oldVal-rmin)/(rmax-rmin))*(tmax-tmin))+tmin
print(oldVal)
print(newVal)
dogPheno_BodySize = newVal
dogJSON['panels'][0]['top']['genericValues'] = np.array([dogPheno_BodySize])

# WhiteSpotting
phenoWhiteSpotting = pd.read_csv(args.whitephenofile)
dogPheno_WhiteSpotting = phenoWhiteSpotting[phenoWhiteSpotting['id'].astype('str') == str(args.studyid)]['prediction'].values[0]
dogJSON['panels'][3]['top']['genericValues'] = np.array([dogPheno_WhiteSpotting])

# results: 'traits'

# BodySize
genoBodySize = pd.read_csv(args.sizegenofile)
genoBodySize['effect'] = genoBodySize['effect'].round(0)
dogGeno_BodySize = genoBodySize[genoBodySize['id'].astype('str') == str(args.studyid)].drop(columns = ['id']).to_dict(orient='records')
dogJSON['panels'][0]['results']['traits'] = dogGeno_BodySize

# WhiteSpotting
genoWhiteSpotting = pd.read_csv(args.whitegenofile)
genoWhiteSpotting['effect'] = genoWhiteSpotting['effect'].round(0)
dogGeno_WhiteSpotting = genoWhiteSpotting[genoWhiteSpotting['id'].astype('str') == str(args.studyid)].drop(columns = ['id']).to_dict(orient='records')
dogJSON['panels'][3]['results']['traits'] = dogGeno_WhiteSpotting

#================================================================
# TRAIT PREDICTIONS: SIMPLE MODELS
#================================================================

# load tables
genoTraits = pd.read_csv(args.genotypes)
dogGenoTraits = genoTraits[genoTraits['sample'].astype('str') == str(args.studyid)]
phenoTraits = pd.read_csv(args.phenotypes)
dogPhenoTraits = phenoTraits[phenoTraits['sample'].astype('str') == str(args.studyid)]

# HTML 'content'
if args.html == "noHTML":
    # CoatColor
    dogJSON['panels'][1]['top']['content'] = "null"
    dogJSON['panels'][1]['results']['content'] = "null"

    # CoatPattern
    dogJSON['panels'][2]['top']['content'] = "null"
    dogJSON['panels'][2]['results']['content'] = "null"

    # CoatType
    dogJSON['panels'][4]['top']['content'] = "null"
    dogJSON['panels'][4]['results']['content'] = "null"

    # SpecialFeatures
    dogJSON['panels'][5]['top']['content'] = "null"
    dogJSON['panels'][5]['results']['content'] = "null"
else:
    # CoatColor
    dogJSON['panels'][1]['top']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.CoatColor.top.content.html").read().replace('\n', "")
    dogJSON['panels'][1]['results']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.CoatColor.results.content.html").read().replace('\n', "")

    # CoatPattern
    dogJSON['panels'][2]['top']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.CoatPattern.top.content.html").read().replace('\n', "")
    dogJSON['panels'][2]['results']['content'] = open("/seq/vgb/dap/json/html/GenomicReports.CoatPattern.results.content.html").read().replace('\n', "")

# top: 'genericValues'

# CoatColor
genericValues = []
for index, row in dogPhenoTraits[dogPhenoTraits['tab']=="CoatColor"].iterrows():
    genericValues.append(row['string'])
    genericValues.append(row['color'])
dogJSON['panels'][1]['top']['genericValues'] = np.array(genericValues)

# CoatPattern
genericValues = []
for index, row in dogPhenoTraits[dogPhenoTraits['tab']=="CoatPattern"][dogPhenoTraits['image'].notnull()].iterrows():
    genericValues.append(row['image'])
    genericValues.append(row['string'])
dogJSON['panels'][2]['top']['genericValues'] = np.array(genericValues)

# CoatType
genericValues = []
for index, row in dogPhenoTraits[dogPhenoTraits['tab']=="CoatType"].iterrows():
    genericValues.append(row['trait'])
    genericValues.append(row['image'])
    genericValues.append(row['string'])
dogJSON['panels'][4]['top']['genericValues'] = np.array(genericValues)

# SpecialFeatures
genericValues = []
for index, row in dogPhenoTraits[dogPhenoTraits['tab']=="SpecialFeatures"].iterrows():
    genericValues.append(row['trait'])
    genericValues.append(row['image'])
    genericValues.append(row['string'])
dogJSON['panels'][5]['top']['genericValues'] = np.array(genericValues)

# results: 'traits'

# CoatColor
dogJSON['panels'][1]['results']['traits'] = dogGenoTraits[dogGenoTraits['tab']=="CoatColor"].drop(columns = ['sample','tab']).to_dict(orient = 'records')

# CoatPattern
dogJSON['panels'][2]['results']['traits'] = dogGenoTraits[dogGenoTraits['tab']=="CoatPattern"].drop(columns = ['sample','tab']).to_dict(orient = 'records')

# CoatType
dogJSON['panels'][4]['results']['traits'] = dogGenoTraits[dogGenoTraits['tab']=="CoatType"].drop(columns = ['sample','tab']).to_dict(orient = 'records')

# SpecialFeatures
dogJSON['panels'][5]['results']['traits'] = dogGenoTraits[dogGenoTraits['tab']=="SpecialFeatures"].drop(columns = ['sample','tab']).to_dict(orient = 'records')

#================================================================
# OUTPUT
#================================================================
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

with open("StudyID-" + str(args.studyid) + "_GenomicReport.json", 'w') as outfile:
    json.dump(obj = dogJSON, fp = outfile, sort_keys=False, indent=4, cls=NpEncoder)

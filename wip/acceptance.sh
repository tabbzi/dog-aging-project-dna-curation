#!/bin/bash

# for file in testing/acceptance/*.vcf.gz ; do
for file in "testing/acceptance/black.vcf.gz" ; do
    basename=${file##*/}
    basename=${basename%%.*}
    echo ${basename}

    mkdir -p testing/acceptance/${basename}/{old,new}

    # run new version
    time /tigress/tcomi/.conda/dap/bin/python genotype_to_phenotype.py \
      --vcf ${file} \
      --output testing/acceptance/${basename}/new/test \
      --variants /home/tcomi/projects/dog-aging-project-dna-curation/ref/DogAgingProject_VariantsOfInterest.csv \
      --reference CanFam3.1 \
      --imputation impute-v2

    # run old version
    time /tigress/tcomi/.conda/dap/bin/python \
        -W ignore \
        GenoToPhenoNew.py \
        -F $file \
        -O testing/acceptance/${basename}/old/test \
        -V ../ref/DogAgingProject_VariantsOfInterest.csv

    # compare
    for f in testing/acceptance/${basename}/old/* ; do
        echo ${f##*/}
        diff \
            <(sort testing/acceptance/${basename}/old/${f##*/}) \
            <(sort testing/acceptance/${basename}/new/${f##*/} | sed 's/\.0$//')
    done

done

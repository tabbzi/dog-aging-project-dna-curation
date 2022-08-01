#!/usr/bin/python

import sys, argparse, os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="This script will ")

parser.add_argument("-S",
                    "--samples",
                    help="query sample IDs",
                    default="query.txt")

parser.add_argument("-Q",
                    "--admQ",
                    help="ADMIXTURE .Q file with labels",
                    default="query.Q.labeled")

parser.add_argument("-P",
                    "--pops",
                    help="RedCap population key",
                    default="/seq/vgb/dap/bin/DAPBreedKey.tsv")

args = parser.parse_args()

populations = pd.read_csv(args.pops, sep = " ")

samples = np.array(open(args.samples).read().splitlines())

adm = pd.read_csv(args.admQ, header=None, delim_whitespace=True)

# define populations in proper order
admDog = adm[adm[1].isin(samples)]
admRef = adm[~adm[1].isin(samples)]
refPop = admRef[0].unique()

admDog.columns = np.append(['FID','IID'],refPop)


admDogLong = pd.melt(admDog, id_vars = ['FID','IID'], value_vars = refPop, var_name = 'breed', value_name = 'percent')

admDogLong = admDogLong[admDogLong['percent']>0.01].reset_index()
# admDogLongOut = admDogLong
admDogLongOut = admDogLong[admDogLong['FID']=="0"]
admDogLongOut['percent'] = admDogLongOut['percent'].round(decimals=4)
for dog, data in admDogLongOut.groupby('IID'):
    data[['breed','percent']].sort_values(by = 'percent').to_csv("StudyID-{}.fullResults.GlobalAncestry.csv".format(dog), index = False)

# remap populations not in RedCap that have close populations:
# admDogLong['breed'] = admDogLong['breed'].replace({'yorkshire_terrier': 'oopsie'})
admDogLong['breed'] = admDogLong['breed'].replace({'landseer': 'newfoundland'})

admDogLong['breed'] = admDogLong['breed'].replace({'white_swiss_shepherd_dog': 'german_shepherd_dog'})

admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_australia': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_china': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_india': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_namibia': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_nigeria': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_portugal': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_qatar': 'village_dog'})
admDogLong['breed'] = admDogLong['breed'].replace({'village_dog_vietnam': 'village_dog'})

# admDogLong['breed'] = admDogLong['breed'].replace({'coyote': 'wolf'})
admDogLong['breed'] = admDogLong['breed'].replace({'wolf_north_american': 'wolf'})
admDogLong['breed'] = admDogLong['breed'].replace({'wolf_eurasian': 'wolf'})

# combine same breed by summing percent
admDogLong = admDogLong.groupby(['FID','IID','breed'])['percent'].sum().reset_index()

# map RedCap population index:
admDogMap = admDogLong.merge(populations, how = 'inner', on = 'breed')

# round percents:
admDogMap['percent'] = admDogMap['percent'].round(decimals=4)

# drop refs
admDogMap = admDogMap[admDogMap['FID']=="0"]

# save each IID results as id, breed, percent
for dog, data in admDogMap.groupby('IID'):
    data[['id','breed','percent']].sort_values(by = 'percent').to_csv("StudyID-{}.GlobalAncestry.csv".format(dog), index = False)
    data[['id','percent']].rename(columns = {'id': 'breed'}).to_json("StudyID-{}.GlobalAncestry.json".format(dog), orient='records')

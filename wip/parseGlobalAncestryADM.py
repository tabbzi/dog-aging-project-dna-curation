#!/usr/bin/python

import argparse
import os
import sys

import numpy as np
import pandas as pd


def parse(args=None):
    parser = argparse.ArgumentParser(
        description=("This script will parse admixture outputs "
                     "to csv and json files for DAP"),
    )

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
                        )

    parser.add_argument("-O",
                        "--basename",
                        help="Base directory to output files",
                        default="./")

    return parser.parse_args(args)


def read_data(args):
    samples = np.loadtxt(args.samples, dtype=str)

    admixture = pd.read_csv(args.admQ, header=None, delim_whitespace=True)

    # define populations in proper order
    refPop = admixture[~admixture[1].isin(samples)][0].unique()
    admixture = admixture[admixture[1].isin(samples)]

    admixture.columns = np.append(['FID', 'IID'], refPop)

    admixture = pd.melt(
        admixture,
        id_vars=['FID', 'IID'],
        value_vars=refPop,
        var_name='breed',
        value_name='percent',
    )

    return admixture[admixture['percent']>0.01].reset_index()


def save_full_results(data, basename):
    output = data[data['FID']=="0"]
    output['percent'] = output['percent'].round(decimals=4)

    for dog, data in output.groupby('IID'):
        data.sort_values(by = 'percent')\
            .to_csv(
            f"{basename}StudyID-{dog}.fullResults.GlobalAncestry.csv",
            index = False,
            columns=['breed', 'percent'],
        )


def map_to_redcap(data, population_file):
    populations = pd.read_csv(population_file, sep=" ")

    # remap populations not in RedCap that have close populations:
    data['breed'] = data['breed'].replace(
        regex={
            r'^landseer$': 'newfoundland',
            r'^white_swiss_shepherd_dog$': 'german_shepherd_dog',
            r'^village_dog_.*$': 'village_dog',
            r'^wolf.*$': 'wolf',
        }
    )

    # combine same breed by summing percent
    data = data.groupby(
        ['FID', 'IID', 'breed'],
    )['percent'].sum().reset_index()

    # map RedCap population index:
    result = data.merge(populations, how='inner', on='breed')

    # round percents:
    result['percent'] = result['percent'].round(decimals=4)

    # drop refs
    return result[result['FID']=="0"]


def save_recap_results(redcap_data, basename):
    # save each IID result as id, breed, percent
    for dog, data in redcap_data.groupby('IID'):
        data.sort_values(by='percent')\
            .to_csv(
                f"{basename}StudyID-{dog}.GlobalAncestry.csv",
                index=False,
                columns=['id', 'breed', 'percent'],
            )
        data[['id', 'percent']]\
            .rename(columns = {'id': 'breed'})\
            .to_json(
                f"{basename}StudyID-{dog}.GlobalAncestry.json",
                orient='records',
            )


if __name__ == "__main__":
    args = parse()

    admixture_data = read_data(args)

    save_full_results(admixture_data, args.basename)

    redcap_data = map_to_redcap(admixture_data, args.pops)

    save_recap_results(redcap_data, args.basename)

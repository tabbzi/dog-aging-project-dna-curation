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
import argparse
import numpy as np
import pandas as pd
import json
from contextlib import contextmanager
import requests

def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=("This script will take a template JSON and results for a "
        "query dog, and populate an output JSON"))

    parser.add_argument("-S", "--study-ids",
                        help="File with study ids to analyze",
                        )

    parser.add_argument("-A", "--ancestry",
                        help="input ancestry table or json with fields `breed` and `percent`, use {id} as id placeholder",
                        default="StudyID-{id}.GlobalAncestry.json")

    parser.add_argument("-I", "--inbreeding",
                        help="input coefficient of inbreeding value or file of COIs",
                        default="COI.csv")

    parser.add_argument("-SP", "--body-size-phenotypes",
                        help="input phenotypes table for complex trait prediction of body size",
                        default="BodySize.phenotypes.csv")

    parser.add_argument("-SG", "--body-size-genotypes",
                        help="input genotypes table for complex trait prediction of body size",
                        default="BodySize.genotypes.csv")

    parser.add_argument("-WP", "--white-spotting-phenotypes",
                        help="input phenotypes table for complex trait prediction of white spotting",
                        default="WhiteSpotting.phenotypes.csv")

    parser.add_argument("-WG", "--white-spotting-genotypes",
                        help="input genotypes table for complex trait prediction of white spotting",
                        default="WhiteSpotting.genotypes.csv")

    parser.add_argument("-P", "--phenotypes",
                        help="input phenotypes table for simple trait predictions, collated by `tab` and `trait`",
                        default="phenotypeTable.csv")

    parser.add_argument("-G", "--genotypes",
                        help="input genotypes table for simple trait predictions, collated by `tab` and `trait`",
                        default="jsonTable.csv")

    parser.add_argument("-O", "--output",
                        help="Output file location with {id} as the id placeholder",
                        default="StudyID-{id}_GenomicReport.json")

    parser.add_argument("-H", "--webhook",
                        help="Webhook to signal, if unset will not make any requests",
                        default=None)


    return parser.parse_args(args)


class JsonBuilder:
    def __init__(self, args):
        self.coi = pd.read_csv(args.inbreeding, dtype={'id': str})
        self.ancestry = pd.read_json(args.ancestry, orient="records", dtype=False)

        self.size_pheno = pd.read_csv(args.body_size_phenotypes,
                                      dtype={'id': str})
        # rescale from [0,4] to [0.25,3.75] on [0,4] scale
        self.size_pheno['prediction'] = (3.5 * self.size_pheno['prediction'] + 1) / 4

        self.size_geno = pd.read_csv(args.body_size_genotypes,
                                      dtype={'id': str})
        self.size_geno['effect'] = self.size_geno['effect'].round(0)

        self.pheno_traits = pd.read_csv(args.phenotypes,
                                        dtype={'sample': str})
        self.geno_traits = pd.read_csv(args.genotypes,
                                        dtype={'sample': str})

        self.white_pheno = pd.read_csv(args.white_spotting_phenotypes,
                                       dtype={'id': str})
        self.white_geno = pd.read_csv(args.white_spotting_genotypes,
                                       dtype={'id': str})
        self.white_geno['effect'] = self.white_geno['effect'].round(0)

    def build_json(self, id):
        geno = self.geno_traits[self.geno_traits['sample'] == id]
        pheno = self.pheno_traits[self.pheno_traits['sample'] == id]
        return {
            'id': id,
            'inbreeding': self.get_coi(id),
            'data': self.get_ancestry(id),
            'panels': [
                self.body_size(id),
                self.coat_color(geno, pheno),
                self.coat_pattern(geno, pheno),
                self.white_spotting(id),
                self.coat_type(geno, pheno),
                self.special_features(geno, pheno),
            ],
        }

    def get_coi(self, id):
        # get first matching id
        return self.coi.loc[self.coi['id'] == id, 'coi'].iloc[0]

    def get_ancestry(self, id):
        return json.loads(self.ancestry.loc[
                self.ancestry.IID == id, ['breed', 'percent']
            ].to_json(orient="records"))

    def body_size(self, id):
        top_values = [self.size_pheno.loc[self.size_pheno.id == id, 'prediction'].iloc[0]]

        result_traits = self.size_geno[self.size_geno['id'] == id]
        result_traits = result_traits.drop(columns=['id'])

        return self._panel_base("body-size", top_values, result_traits)

    def coat_color(self, geno, pheno):
        return self._trait_panel("coat-color", geno, pheno, columns=['string', 'color'], name='Colors')

    def coat_pattern(self, geno, pheno):
        pheno = pheno[pheno.image.notnull()]
        return self._trait_panel("coat-pattern", geno, pheno, columns=['image', 'string'])

    def white_spotting(self, id):
        top_values = [self.white_pheno.loc[self.white_pheno.id == id, 'prediction'].iloc[0]]

        result_traits = self.white_geno[self.white_geno['id'] == id]
        result_traits = result_traits.drop(columns=['id'])

        return self._panel_base("white-spotting", top_values, result_traits)

    def coat_type(self, geno, pheno):
        return self._trait_panel("coat-type", geno, pheno)

    def special_features(self, geno, pheno):
        return self._trait_panel("special-features", geno, pheno)

    def _trait_panel(self, id, geno, pheno, columns=['trait', 'image', 'string'], name=None):
        return self._panel_base(
            id,
            *self._filter_frames(id, geno, pheno, columns=columns),
            name=name,
        )

    @staticmethod
    def _filter_frames(id, geno, pheno, columns):
        tab = id.title().replace('-', '')
        return (
            pheno.loc[pheno.tab == tab, columns],
            geno[geno.tab == tab].drop(columns=['sample', 'tab'])
        )

    @staticmethod
    def _panel_base(id, top_values, result_traits, name=None):
        if not isinstance(top_values, list):
            top_values = top_values.values.flatten().tolist()

        if name is None:
            # "coat-type" to "Coat Type"
            name = id.replace('-', ' ').title()

        return {
            'name': name,
            'id': id,
            'top': {
                'title': None,
                'content': "null",
                'genericValues': top_values,
            },
            'results': {
                'title': None,
                'content': "null",
                'traits': result_traits.to_dict(orient = 'records'),
                'genericValues': None,
            },
        }


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        return super(NpEncoder, self).default(obj)

def main():
    args = parse_args()
    builder = JsonBuilder(args)

    for id in open(args.study_ids):
        id = id.strip()
        if id.startswith('Friend') or ':' in id:
            print(f'Excluding {id} from analysis.')
            continue
        with open(args.output.format(id=id), 'w') as outfile:
            data = builder.build_json(id)
            json.dump(data,
                      fp=outfile,
                      sort_keys=False,
                      indent=4,
                      cls=NpEncoder,
                      )
            if args.webhook:
                request = requests.pos(args.webhook, json=data)
                if request.status_code != 200:
                    print(f"Failed to post sample '{id}'. "
                          f"Response code {request.status_code}")
                else:
                    print(f"Posted sample '{id}'.")


if __name__ == "__main__":
    main()

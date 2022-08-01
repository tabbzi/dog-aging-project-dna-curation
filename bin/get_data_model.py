#!/usr/bin/python

import sys, argparse, os
import firecloud.api as fapi
import numpy as np
import pandas as pd
import json
import datetime

parser = argparse.ArgumentParser(description="This script will export the data model tables for a given workspace and namespace.")

parser.add_argument("-W",
                    "--workspace",
                    help="workspace",
                    default="Dog Aging Project - Sequencing Data")

parser.add_argument("-P",
                    "--project",
                    help="billing project aka namespace",
                    default="dap-user-return-reports-test")

args = parser.parse_args()

platformTable = fapi.get_entities(workspace = args.workspace, namespace = args.project, etype = 'platform').json()
platformTable = pd.json_normalize(platformTable)
platformTable.columns = platformTable.columns.str.replace("attributes.","")
platformTable = platformTable.rename(columns = {"name": "entity:platform_id"}).drop(columns = {'entityType'})

platformTable.to_csv(path_or_buf = "platformTable.tsv", sep = '\t', na_rep = '', header = True, index = False)

sampleTable = fapi.get_entities(workspace = args.workspace, namespace = args.project, etype = 'sample').json()
sampleTable = pd.json_normalize(sampleTable)
sampleTable.columns = sampleTable.columns.str.replace("attributes.","")
sampleTable = sampleTable.rename(columns = {"name": "entity:sample_id"}).drop(columns = {'entityType'})

sampleTable.to_csv(path_or_buf = "sampleTable.tsv", sep = '\t', na_rep = '', header = True, index = False)

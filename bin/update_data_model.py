#!/usr/bin/python

import sys, argparse, os
import firecloud.api as fapi
import numpy as np
import pandas as pd
import json
import datetime

parser = argparse.ArgumentParser(description="This script will use Firecloud API to update the data models and generate JSONs for status updates.")

parser.add_argument("-W",
                    "--workspace",
                    help="workspace",
                    default="Dog Aging Project - Sequencing Data")

parser.add_argument("-P",
                    "--project",
                    help="billing project or namespace",
                    default="dap-user-return-reports-test")

parser.add_argument("-A",
                    "--participant",
                    help="participant data model table",
                    default="participant.tsv")

parser.add_argument("-B",
                    "--sample",
                    help="sample data model table",
                    default="sample.tsv")

parser.add_argument("-C",
                    "--platform",
                    help="platform data model table",
                    default="platform.tsv")

args = parser.parse_args()

# testing:
# workspace="Dog Aging Project - Sequencing Data"
# project="dap-user-return-reports-test"
# test = fapi.get_entities(workspace = workspace, namespace = project, etype = 'platform').json()

# Update data models:
fapi.upload_entities_tsv(model = 'flexible', workspace = args.workspace, namespace = args.project, entities_tsv = args.participant)
fapi.upload_entities_tsv(model = 'flexible', workspace = args.workspace, namespace = args.project, entities_tsv = args.sample)
fapi.upload_entities_tsv(model = 'flexible', workspace = args.workspace, namespace = args.project, entities_tsv = args.platform)

# Read new runs from platform data model table:
platform = pd.read_csv(args.platform, sep = '\t', dtype = {'entity:platform_id': str, 'sample': str, 'client': str, 'participant': str})

# For those with participant != "", generate JSONs
platform = platform[platform['participant']!='']

# make JSON for each NEW sample status
for sample in platform['sample']:
    print("Making JSON for " + sample)
    platform[platform['sample']==sample][['sample','participant','entity:platform_id','datetime','status']].rename(columns = {'sample': 'sample_id', 'participant': 'study_id', 'entity:platform_id': 'gencove_id', 'datetime': 'date_time', 'status': 'sample_status'}).to_json(path_or_buf=str(sample) + ".JSON", orient='records')

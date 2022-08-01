#!/usr/bin/python

import sys, argparse, os
import numpy as np
import pandas as pd
import json

parser = argparse.ArgumentParser(description="This script will generate a sample status JSON for a sample")

parser.add_argument("-B",
                    "--sampleid",
                    help="Sample ID, or barcode",
                    default="test")

parser.add_argument("-P",
                    "--participantid",
                    help="Participant ID, or Study ID",
                    default="test")

parser.add_argument("-G",
                    "--gencoveid",
                    help="Gencove ID",
                    default="test")

parser.add_argument("-D",
                    "--datetime",
                    help="Date and time that status was updated",
                    default="test")

parser.add_argument("-S",
                    "--status",
                    help="Status of sequencing data from sample",
                    default="test")

args = parser.parse_args()

sample = pd.DataFrame({'sample_id': [str(args.sampleid)],'study_id': [str(args.participantid)],'gencove_id': [str(args.gencoveid)],'date_time': [str(args.datetime)],'sample_status': [str(args.status)]})
print(sample)
sample.to_json(path_or_buf=str(args.sampleid) + ".JSON", orient='records')

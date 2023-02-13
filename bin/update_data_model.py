import argparse
import firecloud.api as fapi
import pandas as pd
import json
from io import StringIO
import requests


def parse(args=None):
    parser = argparse.ArgumentParser(description="This script will use Firecloud API to update the data models and generate JSONs for status updates.")

    parser.add_argument("-W",
                        "--workspace",
                        help="workspace",
                        default="Dog Aging Project - Sequencing Data")

    parser.add_argument("-P",
                        "--project",
                        help="billing project or namespace",
                        default="dap-user-return-reports-test")

    parser.add_argument("-S",
                        "--status",
                        help="status from gencove",
                        default="statusTable.tsv")

    parser.add_argument("-K",
                        "--kittable",
                        help="kit table from terra",
                        default="kitTable.tsv")

    parser.add_argument("-H",
                        "--webhook",
                        help="webhook url",
                        default="http://httpbin.org/post")

    return parser.parse_args(args)

def get_entities(args):
    """Get the platform and sample tables as dataframes."""
    def get_table(etype: str):
        result = pd.json_normalize(
                fapi.get_entities(
                    workspace=args.workspace,
                    namespace=args.project,
                    etype=etype).json()
                )
        result.columns = result.columns.str.replace("attributes.", "", regex=False)
        return result.drop(columns = {'entityType'}).rename(
                columns={"name": f"entity:{etype}_id"})

    return get_table('platform'), get_table('sample')


def get_status_table(status):
    result = pd.read_csv(
            status,
            sep='\t',
            names=["datetime", "entity:platform_id", "client", "status", "availability"],
            dtype=str,
            )

    result = result[~result["client"].str.contains('-')]

    return result


def get_kit_table(kit):
    return pd.read_csv(
            kit,
            sep='\t',
            dtype=str,
            ).rename(columns={
                'dog_id': 'participant',
                'sample_id': 'entity:sample_id',
                })


def get_participants(kit, output):
    result = kit[['participant', 'cohort']]
    result = result.rename(columns={
        'participant': 'entity:participant_id',
        'cohort': 'cohort_code'})
    result['int_ids'] = result['entity:participant_id'].astype('int')
    result.sort_values(by='int_ids').drop_duplicates().to_csv(
            output,
            sep='\t',
            na_rep='',
            index=False,
            columns=['entity:participant_id', 'cohort_code'],
            )

def get_sample_updates(kit, sample, output):
    """Get samples in the kit not found in sample."""
    kit[~kit['entity:sample_id'].isin(sample['entity:sample_id'])].to_csv(
            output,
            sep='\t',
            na_rep='',
            index=False,
            columns=['entity:sample_id', 'participant',
                     'date_swab_arrival_laboratory', 'sample_type'],
            )


def get_platform_updates(status, platform, samples, output):
    # get status rows not found in platform
    result = status[~status['entity:platform_id'].isin(platform['entity:platform_id'])]
    # join with kit table, subset
    result = result.merge(
        samples[['entity:sample_id', 'participant']],
        left_on='client',
        right_on='entity:sample_id',
    ).drop(
        columns=['entity:sample_id'],
    )
    result['sample'] = result['client']
    result.to_csv(
        output,
        sep='\t',
        na_rep='',
        index=False,
        columns=[
            'entity:platform_id',
            'client',
            'sample',
            'participant',
            'status',
            'availability',
            'datetime',
        ]
    )


def upload_entities(args, table):  # table can be a file or stringio object
    table.seek(0)
    fapi.upload_entities_tsv(
            model='flexible',
            workspace=args.workspace,
            namespace=args.project,
            entities_tsv=table)


def signal_webhook(args, platform):
    # reset position in platform update
    platform.seek(0)
    # Read new runs from platform data model table:
    to_send = pd.read_csv(
        platform,
        sep='\t',
        dtype=str,
    )

    # For those with participant != "", generate JSONs
    to_send = to_send[to_send['participant']!='']

    # only these columns
    to_send = to_send[['sample','participant','entity:platform_id','datetime','status']]

    # change names
    to_send = to_send.rename(columns={
                                 'sample': 'sample_id',
                                 'participant': 'study_id',
                                 'entity:platform_id': 'gencove_id',
                                 'datetime': 'date_time',
                                 'status': 'sample_status',
                             })

    # send request
    for _, sample in to_send.iterrows():
        request = requests.post(args.webhook, json=sample.to_dict())
        if request.status_code != 200:
            print(f"Failed to post sample '{sample['sample_id']}'. "
                  f"Response code {request.status_code}")
        else:
            print(f"Posted sample '{sample['sample_id']}'.")



if __name__ == "__main__":
    args = parse()
    platform, sample = get_entities(args)
    status = get_status_table(open(args.status))
    kit = get_kit_table(open(args.kittable))

    # update participants table
    participants = StringIO()  # in memory file
    get_participants(kit, participants)
    with open('participant.tsv', 'w') as outfile:
        outfile.write(participants.getvalue())

    # update sample table
    samples = StringIO()  # in memory file
    get_sample_updates(kit, sample, samples)
    with open('sample.tsv', 'w') as outfile:
        outfile.write(samples.getvalue())

    # update platform table
    platform_updates = StringIO()  # in memory file
    get_platform_updates(status, platform, samples, platform_updates)
    with open('platform.tsv', 'w') as outfile:
        outfile.write(platform_updates.getvalue())

    # signal webhook about sample status from platform
    signal_webhook(args, platform_updates)

    # upload changes to gcp
    upload_entities(args, participants)
    upload_entities(args, samples)
    upload_entities(args, platform_updates)

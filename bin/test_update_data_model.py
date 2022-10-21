import update_data_model
from io import StringIO
import numpy as np
import pandas as pd


def test_parse():
    defaults = update_data_model.parse([])
    assert defaults.workspace == "Dog Aging Project - Sequencing Data"
    assert defaults.project == "dap-user-return-reports-test"
    assert defaults.status == "statusTable.tsv"
    assert defaults.kittable == "kitTable.tsv"
    assert defaults.webhook == "http://httpbin.org/post"

    non_defaults = update_data_model.parse((
        '--workspace workspace '
        '--project project '
        '--status status '
        '--kittable kit '
        '--webhook web.hook '
        ).split())
    assert non_defaults.workspace == "workspace"
    assert non_defaults.project == "project"
    assert non_defaults.status == "status"
    assert non_defaults.kittable == "kit"
    assert non_defaults.webhook == "web.hook"


def test_get_entities(mocker):
    sample_json = [{'attributes': {'date_swab_arrival_laboratory': '2022-04-18',
         'participant': 'part1',
         'sample_type': 'saliva_DNA_lowcov'},
         'entityType': 'sample',
         'name': 'sample1'},
     {'attributes': {'date_swab_arrival_laboratory': '2022-04-11',
         'participant': 'part2',
         'sample_type': 'saliva_DNA_lowcov'},
         'entityType': 'sample',
         'name': 'sample2'},
     {'attributes': {'date_swab_arrival_laboratory': '2022-04-04',
         'participant': 'part3',
         'sample_type': 'saliva_DNA_lowcov'},
         'entityType': 'sample',
         'name': 'sample3'},
     {'attributes': {'date_swab_arrival_laboratory': '2022-04-11',
         'participant': 'part4',
         'sample_type': 'saliva_DNA_lowcov'},
         'entityType': 'sample',
         'name': 'sample4'}]

    platform_json = [{'attributes': {'assigned': 'TRUE',
         'assigned_2021': 'FALSE',
         'assigned_2022': 'TRUE',
         'availability': 'available',
         'client': '31020061517115',
         'datetime': '2022-02-17T21:12:42Z',
         'participant': '48055',
         'sample': '31020061517115',
         'sex_xratio_bam': 0.559665,
         'status': 'succeeded'},
         'entityType': 'platform',
         'name': 'name1'},
     {'attributes': {'assigned': 'TRUE',
         'assigned_2021': 'TRUE',
         'assigned_2022': 'TRUE',
         'availability': 'archived',
         'bioproject': 'PRJNA800779',
         'biosample': 'SAMN25333232',
         'client': '31020061518490',
         'datetime': '2021-09-30T17:04:20.936472Z',
         'participant': '59642',
         'release': '2021',
         'sample': '31020061518490',
         'sex_xratio_bam': '0.966975',
         'status': 'succeeded'},
         'entityType': 'platform',
         'name': 'name2'}]

    result = {'sample': mocker.MagicMock(), 'platform': mocker.MagicMock()}
    result['sample'].json.return_value = sample_json
    result['platform'].json.return_value = platform_json

    mock_fapi = mocker.patch('update_data_model.fapi.get_entities',
            side_effect=lambda *args, **kwargs: result[kwargs['etype']]
            )

    args = mocker.MagicMock()
    args.workspace = "workspace"
    args.project = "project"
    platform, sample = update_data_model.get_entities(args)

    mock_fapi.assert_has_calls([
        mocker.call(workspace='workspace', namespace='project', etype='platform'),
        mocker.call(workspace='workspace', namespace='project', etype='sample'),
        ])

    # name is renamed
    assert (platform['entity:platform_id'] == ['name1', 'name2']).all()
    # entityType is removed
    assert 'entityType' not in platform.columns
    # attribute is removed from status
    assert (platform['status'] == ['succeeded', 'succeeded']).all()

    # name is renamed
    assert (sample['entity:sample_id'] == [f'sample{i}' for i in range(1,5)]).all()
    # entityType is removed
    assert 'entityType' not in sample.columns


def test_get_status_table():
    data = StringIO(
            '2022-07-09T19:30:07.949485+00:00\tplatform1\t31211050307003\tsucceeded\tavailable\n'
            '2022-07-09T19:30:06.458277+00:00\tplatform2\t31211050304787\tsucceeded\tavailable\n'
            '2022-07-09T19:30:05.184814+00:00\tplatform3\t31211050304685\tsucceeded\tavailable\n'
            '2022-07-09T19:30:03.465062+00:00\tplatform4\t31211050316980\tsucceeded\tavailable\n'
            '2022-07-09T19:30:01.321240+00:00\tplatform5\t31211050304746\tsucceeded\tavailable\n'
            '2022-07-09T19:30:00.120093+00:00\tplatform6\t31211050307026\tsucceeded\tavailable\n'
            # this format happens for some samples from gencove, should be removed
            '2022-07-09T19:29:58.622768+00:00\tplatform7\tabc-31211050316950-def\tsucceeded\tavailable\n'
            )
    status = update_data_model.get_status_table(data)

    assert (status.columns == ["datetime", "entity:platform_id", "client", "status", "availability"]).all()
    assert (status["entity:platform_id"] == [f"platform{i}" for i in range(1, 7)]).all()


def test_get_kit_table():
    data = StringIO(
'sample_id\tdog_id\tcohort\tsample_type\tdate_swab_arrival_laboratory\n'
'31020061514610\t72814\t10\tsaliva_DNA_lowcov\t2022-02-24\n'
'31020061516291\t7580\t9\tsaliva_DNA_lowcov\t2021-03-08\n'
'31020061513376\t107536\t9\tsaliva_DNA_lowcov\t2022-02-14\n'
'31020061518922\t81878\t10\tsaliva_DNA_lowcov\t2022-01-18\n'
            )
    kits = update_data_model.get_kit_table(data)

    assert (kits.columns == 'entity:sample_id\tparticipant\tcohort\tsample_type\tdate_swab_arrival_laboratory'.split()).all()

def test_get_participants():
    output = StringIO()
    kit = pd.DataFrame.from_dict(
            {
                'participant': [f'{i}' for i in range(5)] + ['20', '10', '10'],
                'cohort': [f'{i}' for i in range(5)] + ['11', '21', '21'],
                'unused': ['z' for _ in range(5)] + ['z']*3,
            }
            )
    update_data_model.get_participants(kit, output)

    result = output.getvalue().split('\n')
    # only these columns, renamed
    assert result[0] == 'entity:participant_id\tcohort_code'
    # sorted numerically
    assert result[1] == '0\t0'
    assert result[2] == '1\t1'
    assert result[3] == '2\t2'
    assert result[4] == '3\t3'
    assert result[5] == '4\t4'
    assert result[6] == '10\t21'  # no duplicaates
    assert result[7] == '20\t11'


def test_get_sample_updates():
    output = StringIO()
    kit = pd.DataFrame.from_dict(
            {
                'entity:sample_id': [f"{i}" for i in range(1, 10)],
                'participant': [f"p{i}" for i in range(1, 10)],
                'date_swab_arrival_laboratory': [f"2022-02-0{i}" for i in range(1, 10)],
                'sample_type': [f"sample" for _ in range(1, 10)],
            }
        )

    sample = pd.DataFrame.from_dict(
            {
                'entity:sample_id': [f"{i}" for i in range(1, 10, 2)],
                'other': [f"o{i}" for i in range(1, 10, 2)],
            }
        )


    update_data_model.get_sample_updates(kit, sample, output)

    result = output.getvalue().split('\n')
    assert result[0] == 'entity:sample_id\tparticipant\tdate_swab_arrival_laboratory\tsample_type'
    assert result[1] == '2\tp2\t2022-02-02\tsample'
    assert result[2] == '4\tp4\t2022-02-04\tsample'
    assert result[3] == '6\tp6\t2022-02-06\tsample'
    assert result[4] == '8\tp8\t2022-02-08\tsample'
    assert result[5] == ''
    assert len(result) == 6


def test_get_platform_updates():
    output = StringIO()

    kit = pd.DataFrame.from_dict(
        {
            'entity:sample_id': [f"client{i}" for i in range(1, 10)],
            'participant': [f"p{i}" for i in range(1, 10)],
            'date_swab_arrival_laboratory': [f"2022-02-0{i}" for i in range(1, 10)],
            'sample_type': [f"sample" for _ in range(1, 10)],
        }
    )
    status = pd.DataFrame.from_dict(
        {
            "datetime": [f"date{i}" for i in range(1,5)],
            "entity:platform_id": [f"id{i}" for i in range(1, 5)],
            "client": [f"client{i}" for i in range(1, 5)],
            "status": [f"status{i}" for i in range(1, 5)],
            "availability": [f"avail{i}" for i in range(1, 5)],
        }
    )

    platform = pd.DataFrame.from_dict(
        {
            'entity:platform_id': [f"id{i}" for i in range(2, 6, 2)],
            'status': ["asdf" for _ in range(2, 6, 2)],
        }
    )

    update_data_model.get_platform_updates(status, platform, kit, output)
    result = output.getvalue().split('\n')
    assert result[0] == 'entity:platform_id\tclient\tsample\tparticipant\tstatus\tavailability\tdatetime'
    assert result[1] == 'id1\tclient1\tclient1\tp1\tstatus1\tavail1\tdate1'
    assert result[2] == 'id3\tclient3\tclient3\tp3\tstatus3\tavail3\tdate3'
    assert result[3] == ''
    assert len(result) == 4


def test_signal_webhook(mocker, capsys):
    platform_updates = StringIO(
        'entity:platform_id\tclient\tsample\tparticipant\tstatus\tavailability\tdatetime\n'
        'id1\tclient1\tclient1\tp1\tstatus1\tavail1\tdate1\n'
        'id3\tclient3\tclient3\tp3\tstatus3\tavail3\tdate3\n'
    )

    mock_result = mocker.MagicMock()
    mock_result.status_code = 200
    mock_post = mocker.patch('update_data_model.requests.post',
                             return_value=mock_result)

    args = update_data_model.parse([])
    update_data_model.signal_webhook(args, platform_updates)

    stdout = capsys.readouterr().out.split('\n')
    assert stdout[0] == "Posted sample 'client1'."
    assert stdout[1] == "Posted sample 'client3'."

    mock_post.assert_has_calls(
        [
            mocker.call(args.webhook, json={
                            'date_time': 'date1',
                            'gencove_id': 'id1',
                            'sample_id': 'client1',
                            'sample_status': 'status1',
                            'study_id': 'p1',
                        }),
            mocker.call(args.webhook, json={
                            'date_time': 'date3',
                            'gencove_id': 'id3',
                            'sample_id': 'client3',
                            'sample_status': 'status3',
                            'study_id': 'p3',
                        }),
        ]
    )


def test_signal_webhook_error(mocker, capsys):
    platform_updates = StringIO(
        'entity:platform_id\tclient\tsample\tparticipant\tstatus\tavailability\tdatetime\n'
        'id1\tclient1\tclient1\tp1\tstatus1\tavail1\tdate1\n'
        'id3\tclient3\tclient3\tp3\tstatus3\tavail3\tdate3\n'
    )

    mock_result = mocker.MagicMock()
    mock_result.status_code = 234
    mock_post = mocker.patch('update_data_model.requests.post',
                             return_value=mock_result)

    args = update_data_model.parse([])
    update_data_model.signal_webhook(args, platform_updates)

    stdout = capsys.readouterr().out.split('\n')
    assert stdout[0] == ("Failed to post sample 'client1'. "
                         "Response code 234")
    assert stdout[1] == ("Failed to post sample 'client3'. "
                         "Response code 234")

    mock_post.assert_has_calls(
        [
            mocker.call(args.webhook, json={
                            'date_time': 'date1',
                            'gencove_id': 'id1',
                            'sample_id': 'client1',
                            'sample_status': 'status1',
                            'study_id': 'p1',
                        }),
            mocker.call(args.webhook, json={
                            'date_time': 'date3',
                            'gencove_id': 'id3',
                            'sample_id': 'client3',
                            'sample_status': 'status3',
                            'study_id': 'p3',
                        }),
        ]
    )

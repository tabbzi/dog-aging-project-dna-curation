import pytest
import ReportsToJSON as r2j
import json


@pytest.fixture
def default_args(roh_coi_csv, report_body_pheno, report_body_geno,
                 report_pheno, report_json_table, report_white_pheno,
                 report_white_geno, get_adm):
    args = r2j.parse_args([])
    args.inbreeding = roh_coi_csv
    args.ancestry = get_adm
    args.body_size_phenotypes = report_body_pheno
    args.body_size_genotypes = report_body_geno 
    args.genotypes = report_json_table 
    args.phenotypes = report_pheno
    args.white_spotting_phenotypes = report_white_pheno
    args.white_spotting_genotypes = report_white_geno 
    return args


@pytest.fixture
def builder(default_args):
    return r2j.JsonBuilder(default_args)

def test_builder_get_coi(builder):
    assert builder.get_coi('1') == 0.179865
    assert builder.get_coi('22') == 0.235657
    assert builder.get_coi('333') == 0.0801061

def test_builder_build_json(builder, mocker, json_1, json_22, json_333):
    result = json.load(json_1)

    assert builder.build_json('1') == result

    assert builder.build_json('22') == json.load(json_22)

    assert builder.build_json('333') == json.load(json_333)

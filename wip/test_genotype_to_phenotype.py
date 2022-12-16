import genotype_to_phenotype as g2p
from io import StringIO
import pytest
import numpy as np
import pandas as pd
import pathlib
from pprint import pprint
from itertools import product


@pytest.fixture()
def default_args(mocker, tmp_path, variants_of_interest, sample_vcf):
    args = mocker.MagicMock()
    args.variants = variants_of_interest

    vcf = str(tmp_path / "test.vcf")

    with open(vcf, 'w') as temp:
        temp.write(sample_vcf.getvalue())
    args.vcf = vcf
    args.imputation = "impute-v2"
    args.reference = "CanFam3.1"
    args.output = str(tmp_path / "output")
    return args

def test_read_variant_table(default_args):
    vt = g2p.read_variant_table(default_args)
    assert len(vt) == 44
    assert vt['CHR'].dtype == int
    assert vt['POS'].dtype == int


def test_read_vcf(default_args):
    result = g2p.read_vcf(default_args.vcf)
    assert len(result['variants/POS']) == 151
    assert result['variants/CHR'].dtype == int
    assert result['variants/CHR'][-1] == 38  # strip leading chr


@pytest.fixture
def sample_caller(default_args):
    vcf = g2p.read_vcf(default_args.vcf)
    variants = g2p.read_variant_table(default_args)
    return g2p.TraitCaller(variants, vcf)

def test_caller_initialization(sample_caller):
    assert sample_caller.mask.shape == (44, 151)
    # handle ./. properly; gt == NA
    assert (sample_caller.gts.iloc[-1, :] == (['NA', 'NA'] + ['G']*8)).all()
    # not a variant count
    assert (sample_caller.variant_counts.iloc[-1, :] == 0).all()

def assert_result_equals(results, traits):
    for (dog, result), trait in zip(results.iterrows(), traits):
        assert result['result'] == trait.id, f'failed for id with dog {dog}'
        assert result['string'] == trait.display, f'failed for display with dog {dog}'
        assert result['value'] == trait.value, f'failed for value with dog {dog}'
        assert result['color'] == trait.hex, f'failed for hex with dog {dog}'
        assert result['image'] == trait.image, f'failed for image with dog {dog}'

def test_caller_trait_black(sample_caller):
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_black()
    assert_result_equals(
        result, [
            g2p.TRAITS['black'],
            g2p.TRAITS['dilute_liver_leathers'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['dilute_liver_leathers'],
        ])

def test_decode_black_trait():
    modules = pd.DataFrame(
        list(product(
                 [True, False],
                 [True, False],
                 [True, False],
                 [True, False],
                 [True, False],
                 ['Y', 'C', 'N'],
             )),
        columns=[
                'liver',
                'cocoa',
                'dilute',
                'dom_black',
                'rec_black',
                'rec_red',
        ]
    )
    result = g2p.TraitCaller.decode_black_trait(modules)
    assert_result_equals(
        result, [
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa'],
            g2p.TRAITS['dilute_cocoa'],
            #
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa'],
            g2p.TRAITS['cocoa'],
            #
            g2p.TRAITS['dilute_liver_leathers'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_leathers'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_leathers'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_solid'],
            g2p.TRAITS['dilute_liver_leathers'],
            g2p.TRAITS['dilute_liver'],
            g2p.TRAITS['dilute_liver'],
            #
            g2p.TRAITS['liver_leathers'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_leathers'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_leathers'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_solid'],
            g2p.TRAITS['liver_leathers'],
            g2p.TRAITS['liver'],
            g2p.TRAITS['liver'],
            #
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_solid'],
            g2p.TRAITS['dilute_cocoa_leathers'],
            g2p.TRAITS['dilute_cocoa'],
            g2p.TRAITS['dilute_cocoa'],
            #
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_solid'],
            g2p.TRAITS['cocoa_leathers'],
            g2p.TRAITS['cocoa'],
            g2p.TRAITS['cocoa'],
            #
            g2p.TRAITS['dilute_black_leathers'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_leathers'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_leathers'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_solid'],
            g2p.TRAITS['dilute_black_leathers'],
            g2p.TRAITS['dilute_black'],
            g2p.TRAITS['dilute_black'],
            #
            g2p.TRAITS['black_leathers'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_leathers'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_leathers'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_solid'],
            g2p.TRAITS['black_leathers'],
            g2p.TRAITS['black'],
            g2p.TRAITS['black'],
        ])

def test_trait_red(sample_caller):
    # make up interesting data
    modules = pd.DataFrame(
        list(product(
                 [0, 1, 2],
                 [0, 1, 2],
                 [0, 1, 2],
                 [0, 1, 2],
                 [0, 1, 2],
             )),
    )
    sample_caller.variant_counts = pd.DataFrame(
        0,
        index=range(44),
        columns=pd.Index((str(r) for r in range(len(modules))), name='sample'))
    sample_caller.variant_counts.loc[sample_caller.variants['module'] == 'mod_red_intensity', :] = modules.to_numpy().T
    result = sample_caller.trait_red()
    expected = [g2p.TRAITS['cream']] * len(result)

    # this is ugly but space efficient
    for index in [
        4,   5,   6,   7,   8,  12,  13,  14,  15,  21,  22,  23,  24,
        30,  31,  32,  33,  34,  39,  40,  41,  48,  49,  50,  51,  57,
        58,  59,  60,  66,  67,  68,  75,  76,  77,  84,  85,  86,  91,
        92,  93, 100, 101, 102, 111, 112, 113, 117, 118, 119, 120, 127,
        128, 129, 136, 137, 138, 144, 145, 146, 153, 154, 155, 162, 163,
        164, 171, 172, 180, 181, 182, 189, 190, 191, 198, 207, 216, 217,
        218, 234]:
        expected[index] = g2p.TRAITS['yellow']

    for index in [
        16,  17,  25,  26,  35,  42,  43,  44,  52,  53,  61,  62,  69,
        70,  71,  78,  79,  80,  87,  88,  89,  94,  95,  96,  97,  98,
        103, 104, 105, 106, 107, 114, 115, 116, 121, 122, 123, 124, 125,
        130, 131, 132, 133, 134, 139, 140, 141, 142, 143, 147, 148, 149,
        150, 156, 157, 158, 159, 165, 166, 167, 173, 174, 175, 176, 183,
        184, 185, 192, 193, 194, 199, 200, 201, 208, 209, 210, 211, 219,
        220, 221, 225, 226, 227, 228, 235, 236, 237]:
        expected[index] = g2p.TRAITS['tan']

    for index in [
        151, 152, 160, 161, 168, 169, 170, 177, 178, 179, 186, 187, 188,
        195, 196, 197, 202, 203, 204, 205, 206, 212, 213, 214, 215, 222,
        223, 224, 229, 230, 231, 232, 233, 238, 239, 240, 241, 242]:
        expected[index] = g2p.TRAITS['red']

    assert_result_equals(result, expected)


def test_caller_trait_agouti(sample_caller):
    result = sample_caller.trait_agouti()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['sable'],
            g2p.TRAITS['tan_points_hidden'],
            g2p.TRAITS['tan_points_hidden'],
            g2p.TRAITS['sable_hidden'],
            g2p.TRAITS['tan_points'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_agouti()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['agouti'],
            g2p.TRAITS['sable_hidden'],
            g2p.TRAITS['sable_hidden'],
            g2p.TRAITS['sable_hidden'],
            g2p.TRAITS['sable_hidden'],
        ])

def test_caller_trait_extension(sample_caller):
    result = sample_caller.trait_extension()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['mask'],
            g2p.TRAITS['rec_red'],
            g2p.TRAITS['mask_hidden'],
            g2p.TRAITS['normal_extension_hidden'],
            g2p.TRAITS['mask'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_extension()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_extension'],
            g2p.TRAITS['mask_hidden'],
            g2p.TRAITS['mask_hidden'],
            g2p.TRAITS['grizzle_hidden'],
            g2p.TRAITS['mask_hidden'],
        ])


def test_decode_extension_trait():
    modules = pd.DataFrame(
        [
            [True, True, False, 0, 0, 1],
            [True, False, False, 0, 2, 2],
            [False, False, True, 2, 0, 1],
            [False, False, False, 2, 1, 1],
            [False, True, True, 1, 2, 1],
            [False, False, False, 0, 2, 1],
            [False, False, False, 1, 1, 2],
            [False, True, False, 1, 1, 2],
            [False, True, False, 0, 1, 2],
            [False, True, True, 0, 1, 2],
            [False, False, True, 1, 0, 1],
            [False, False, False, 1, 0, 1],
            [False, True, False, 0, 1, 1],
            [False, False, False, 0, 1, 1],
            [False, True, False, 0, 1, 0],
            [False, False, False, 0, 1, 0],

        ],
        columns=[
                'mask',
                'dom_black',
                'rec_black',
                'grizzle',
                'domino',
                'rec_red',
        ]
    )
    result = g2p.TraitCaller.decode_extension_trait(modules)
    assert_result_equals(
        result,
        [
            g2p.TRAITS['mask_hidden'],
            g2p.TRAITS['mask'],
            g2p.TRAITS['grizzle_hidden'],
            g2p.TRAITS['grizzle'],
            g2p.TRAITS['domino_hidden'],
            g2p.TRAITS['domino'],
            g2p.TRAITS['grizzle_domino'],
            g2p.TRAITS['grizzle_domino'],
            g2p.TRAITS['rec_red'],
            g2p.TRAITS['rec_red'],
            g2p.TRAITS['grizzle_hidden'],
            g2p.TRAITS['grizzle'],
            g2p.TRAITS['domino_hidden'],
            g2p.TRAITS['domino'],
            g2p.TRAITS['normal_extension_hidden'],
            g2p.TRAITS['normal_extension'],
        ])


def test_caller_trait_ticking(sample_caller):
    result = sample_caller.trait_ticking()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_ticking()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['ticking'],
            g2p.TRAITS['no_ticking'],
            g2p.TRAITS['no_ticking'],
        ])


def test_caller_trait_brindle(sample_caller):
    result = sample_caller.trait_brindle()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_brindle()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['brindle_hidden'],
            g2p.TRAITS['brindle_hidden'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['brindle_hidden'],
        ])

def test_decode_brindle_trait():
    modules = pd.DataFrame(
        list(product(
                 [True, False],
                 [True, False],
                 [True, False],
             )),
        columns=[
            'brindle', 'dom_black', 'rec_red',
        ]
    )
    result = g2p.TraitCaller.decode_brindle_trait(modules)
    assert_result_equals(
        result,
        [
            g2p.TRAITS['brindle_hidden'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['brindle_hidden'],
            g2p.TRAITS['brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
            g2p.TRAITS['not_brindle'],
        ])


def test_caller_trait_merle(sample_caller):
    result = sample_caller.trait_merle()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_merle'],
            g2p.TRAITS['not_merle'],
            g2p.TRAITS['not_merle'],
            g2p.TRAITS['not_merle'],
            g2p.TRAITS['not_merle'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_merle()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_merle'],
            g2p.TRAITS['uncalled'],
            g2p.TRAITS['embryonic_lethal'],
            g2p.TRAITS['uncalled'],
            g2p.TRAITS['embryonic_lethal'],
        ])


def test_caller_trait_coat_length(sample_caller):
    result = sample_caller.trait_coat_length()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['short_coat'],
            g2p.TRAITS['short_coat'],
            g2p.TRAITS['long_coat'],
            g2p.TRAITS['short_coat'],
            g2p.TRAITS['short_coat'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_coat_length()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['short_coat'],
            g2p.TRAITS['long_coat'],
            g2p.TRAITS['long_coat'],
            g2p.TRAITS['short_coat'],
            g2p.TRAITS['short_coat'],
        ])

def test_caller_trait_coat_texture(sample_caller):
    result = sample_caller.trait_coat_texture()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['straight_coat'],
            g2p.TRAITS['straight_coat'],
            g2p.TRAITS['wavy_coat'],
            g2p.TRAITS['straight_coat'],
            g2p.TRAITS['straight_coat'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_coat_texture()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['straight_coat'],
            g2p.TRAITS['curly_coat'],
            g2p.TRAITS['curly_coat'],
            g2p.TRAITS['wavy_coat'],
            g2p.TRAITS['curly_coat'],
        ])


def test_caller_trait_furnishings(sample_caller):
    result = sample_caller.trait_furnishings()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['no_furnishings'],
            g2p.TRAITS['no_furnishings'],
            g2p.TRAITS['no_furnishings'],
            g2p.TRAITS['no_furnishings'],
            g2p.TRAITS['no_furnishings'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_furnishings()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['no_furnishings'],
            g2p.TRAITS['furnishings'],
            g2p.TRAITS['furnishings'],
            g2p.TRAITS['furnishings'],
            g2p.TRAITS['furnishings'],
        ])


def test_caller_trait_shedding_propensity(sample_caller):
    result = sample_caller.trait_shedding_propensity()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_shedding'],
            g2p.TRAITS['low_shedding'],
            g2p.TRAITS['low_shedding'],
            g2p.TRAITS['normal_shedding'],
            g2p.TRAITS['low_shedding'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_shedding_propensity()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_shedding'],
            g2p.TRAITS['low_shedding'],
            g2p.TRAITS['low_shedding'],
            g2p.TRAITS['low_shedding'],
            g2p.TRAITS['low_shedding'],
        ])


def test_caller_trait_natural_bob_tail(sample_caller):
    result = sample_caller.trait_natural_bob_tail()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['normal_tail'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_natural_bob_tail()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['natural_bob_tail'],
            g2p.TRAITS['embryonic_lethal'],
            g2p.TRAITS['normal_tail'],
            g2p.TRAITS['normal_tail'],
        ])


def test_caller_trait_leg_length(sample_caller):
    result = sample_caller.trait_leg_length()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_legs'],
            g2p.TRAITS['normal_legs'],
            g2p.TRAITS['normal_legs'],
            g2p.TRAITS['short_legs_marker'],
            g2p.TRAITS['normal_legs'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_leg_length()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['normal_legs'],
            g2p.TRAITS['short_legs_marker'],
            g2p.TRAITS['short_legs_marker'],
            g2p.TRAITS['short_legs_marker'],
            g2p.TRAITS['normal_legs'],
        ])


def test_caller_trait_high_altitude_adaptation(sample_caller):
    result = sample_caller.trait_high_altitude_adaptation()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
        ])
    # make up interesting data
    sample_caller.variant_counts.iloc[:, 0] = 0
    sample_caller.variant_counts.iloc[:, 1] = 1
    sample_caller.variant_counts.iloc[:, 2] = 2
    sample_caller.variant_counts.iloc[::2, 3] = 1
    sample_caller.variant_counts.iloc[::2, 4] = 2
    result = sample_caller.trait_high_altitude_adaptation()
    assert_result_equals(
        result,
        [
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['hypoxia_adapted'],
            g2p.TRAITS['hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
            g2p.TRAITS['not_hypoxia_adapted'],
        ])


def test_predict_phenotypes(sample_caller, old_predict_output):
    old_result = pd.read_csv(old_predict_output, keep_default_na=False, dtype={'sample': str})
    old_result = old_result.sort_values(by=old_result.columns.to_list()).reset_index(drop=True)
    result = sample_caller.predict_phenotypes()
    result = result.sort_values(by=result.columns.to_list()).reset_index(drop=True)

    pd.testing.assert_frame_equal(old_result, result)

def test_main(default_args, old_outputs):
    g2p.main(default_args)

    files = sorted(
        str(file.name)
        for file in pathlib.Path(default_args.output).parent.glob('output*'))

    assert files == [
        'output_genotypeTable.csv',
        'output_jsonTable.csv',
        'output_phenotypeTable.csv',
        'output_trailblazerGenotypeTable.csv',
        'output_trailblazerPhenotypeTable.tsv',
    ]

    for file, value in old_outputs.items():
        # this works if the order doesn't matter
        new_values = sorted(line.rstrip('\n') for line in open(f'{default_args.output}_{file}').readlines())
        old_values = sorted(value.getvalue().split('\n'))
        # for new, old in zip(new_values, old_values):
        #     if old not in new_values:
        #         print(old, 'old missing')
        #     if new not in old_values:
        #         print(new, 'new missing')
        assert new_values == old_values, f"failed on {default_args.output}_{file}"

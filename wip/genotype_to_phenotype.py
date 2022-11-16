#================================================================
# DESCRIPTION
#   This script will take a canine VCF, either single sample or
#   multi-sample, and extract genetic variants of interest,
#   using the consensus genotypes to predict phenotypes for
#   well-studied, simple genetic traits.
#================================================================

import argparse
import numpy as np
import pandas as pd
import allel
import re
from dataclasses import dataclass
import itertools


def parse_args(args=None):
    parser = argparse.ArgumentParser(description="This script extracts and interprets genotype calls from a canine VCF into phenotype predictions for well-studied physical traits.")

    parser.add_argument("-F", "--vcf",
                        help="input variant call file",
                        default="input.vcf.gz")
    parser.add_argument("-V", "--variants",
                        help="input variants of interest, comma-delimited file",
                        default="bin/DogAgingProject_VariantsOfInterest.csv")
    parser.add_argument("-I", "--imputation",
                        help="If imputed, then imputation panel version? OPTS: impute-v1, impute-v2",
                        default="impute-v2")
    parser.add_argument("-R", "--reference",
                        help="Which reference genome? OPTS: CanFam3.1, UU_Cfam_GSD_1.0, Dog10K_Boxer_Tasha",
                        default="CanFam3.1")
    parser.add_argument("-O", "--output",
                        help="output filepath and prefix",
                        default="output")

    args = parser.parse_args(args)

    print(f"Variant call file (VCF): {args.vcf}")
    print(f"Variants of interest: {args.variants}")
    print(f"Asserting that imputation panel is {args.imputation}")
    print(f"Asserting that reference panel is {args.reference}")

    return args


@dataclass
class Trait():
    id: str
    display: str = ""
    hex: str = ""
    value: str = ""
    image: str = ""

    def to_dict(self):
        return {
            'result': self.id,
            'string': self.display,
            'value': self.value,
            'color': self.hex,
            'image': self.image,
        }

# this simplifies testing
TRAITS = {
    'dilute_cocoa_leathers': Trait(id='dilute_cocoa_leathers', display='Dilute Cocoa (nose and paws)', hex='#947D80'),
    'dilute_cocoa_solid': Trait(id='dilute_cocoa_solid', display='Dilute Cocoa (solid coat)', hex='#947D80'),
    'dilute_cocoa': Trait(id='dilute_cocoa', display='Dilute Cocoa', hex='#947D80'),
    'cocoa_leathers': Trait(id='cocoa_leathers', display='Cocoa (nose and paws)', hex='#21100e'),
    'cocoa_solid': Trait(id='cocoa_solid', display='Cocoa (solid coat)', hex='#21100e'),
    'cocoa': Trait(id='cocoa', display='Cocoa', hex='#21100e'),
    'dilute_liver_leathers': Trait(id='dilute_liver_leathers', display='Dilute Brown (nose and paws)', hex='#947D80'),
    'dilute_liver_solid': Trait(id='dilute_liver_solid', display='Dilute Brown (solid coat)', hex='#947D80'),
    'dilute_liver': Trait(id='dilute_liver', display='Dilute Brown', hex='#947D80'),
    'liver_leathers': Trait(id='liver_leathers', display='Brown (nose and paws)', hex='#442825'),
    'liver_solid': Trait(id='liver_solid', display='Brown (solid coat)', hex='#442825'),
    'liver': Trait(id='liver', display='Brown', hex='#442825'),
    'dilute_black_leathers': Trait(id='dilute_black_leathers', display='Dilute Black (nose and paws)', hex='#746F76'),
    'dilute_black_solid': Trait(id='dilute_black_solid', display='Dilute Black (solid coat)', hex='#746F76'),
    'dilute_black': Trait(id='dilute_black', display='Dilute Black', hex='#746F76'),
    'black_leathers': Trait(id='black_leathers', display='Black (nose and paws)', hex='#17181D'),
    'black_solid': Trait(id='black_solid', display='Black (solid coat)', hex='#17181D'),
    'black': Trait(id='black', display='Black', hex='#17181D'),

    'cream': Trait(id='cream', display='Cream', hex='#D8C4AF'),
    'yellow': Trait(id='yellow', display='Yellow', hex='#CE9E63'),
    'tan': Trait(id='tan', display='Tan', hex='#A86B39'),
    'red': Trait(id='red', display='Red', hex='#7E341B'),

    'agouti': Trait(id='agouti', display='Agouti', image='CoatPattern_locusA_agouti.svg'),
    'agouti_hidden': Trait(id='agouti_hidden', display='Agouti (hidden)', image='CoatPattern_locusA_agouti.svg'),
    'rec_black': Trait(id='rec_black', display='No sable, agouti, or tan points'),
    'rec_black_hidden': Trait(id='rec_black_hidden', display='No sable, agouti, or tan points'),
    'sable': Trait(id='sable', display='Sable', image='CoatPattern_locusA_sable.svg'),
    'sable_hidden': Trait(id='sable_hidden', display='Sable (hidden)', image='CoatPattern_locusA_sable.svg'),
    'tan_points': Trait(id='tan_points', display='Tan Points', image='CoatPattern_locusA_tanPoints.svg'),
    'tan_points_hidden': Trait(id='tan_points_hidden', display='Tan Points (hidden)', image='CoatPattern_locusA_tanPoints.svg'),

    'domino': Trait(id='domino', display='Northern Domino', image='CoatPattern_locusE_domino.svg'),
    'domino_hidden': Trait(id='domino_hidden', display='Northern Domino (hidden)', image='CoatPattern_locusE_domino.svg'),
    'grizzle': Trait(id='grizzle', display='Sighthound Grizzle', image='CoatPattern_locusE_grizzle.svg'),
    'grizzle_domino': Trait(id='grizzle_domino', display='Grizzle or Domino', image='CoatPattern_locusE_domino.svg'),
    'grizzle_hidden': Trait(id='grizzle_hidden', display='Sighthound Grizzle (hidden)', image='CoatPattern_locusE_grizzle.svg'),
    'mask': Trait(id='mask', display='Facial Mask', image='CoatPattern_locusE_mask.svg'),
    'mask_hidden': Trait(id='mask_hidden', display='Facial Mask (hidden)', image='CoatPattern_locusE_mask.svg'),
    'normal_extension': Trait(id='normal_extension', display='No mask, grizzle, or domino patterns'),
    'normal_extension_hidden': Trait(id='normal_extension_hidden', display='No mask, grizzle, or domino patterns'),
    'rec_red': Trait(id='rec_red', display='No mask, grizzle, or domino patterns'),
    # these are really non-hidden but including them simplifies decoding
    'rec_red_hidden': Trait(id='rec_red', display='No mask, grizzle, or domino patterns'),
    'grizzle_domino_hidden': Trait(id='grizzle_domino', display='Grizzle or Domino', image='CoatPattern_locusE_domino.svg'),

    'ticking': Trait(id='ticking', display='Ticking or Roaning', image='CoatPattern_ticking.svg'),
    'no_ticking': Trait(id='no_ticking'),

    'not_brindle': Trait(id='not_brindle'),
    'brindle': Trait(id='brindle', display='Brindle', image='CoatPattern_brindle.svg'),
    'brindle_hidden': Trait(id='brindle_hidden', display='Brindle (hidden)', image='CoatPattern_brindle.svg'),

    'not_merle': Trait(id='not_merle'),
    'embryonic_lethal': Trait(id='embryonic_lethal'),
    'uncalled': Trait(id='uncalled'),

    'long_coat': Trait(id='long_coat', display="Long coat", image='CoatType_Length_long.svg'),
    'short_coat': Trait(id='short_coat', display="Short coat", image='CoatType_Length_short.svg'),

    'straight_coat': Trait(id='straight_coat', display="Straight coat", image='CoatType_Curl_straight.svg'),
    'wavy_coat': Trait(id='wavy_coat', display="Wavy coat", image='CoatType_Curl_wavy.svg'),
    'curly_coat': Trait(id='curly_coat', display="Curly coat", image='CoatType_Curl_curly.svg'),

    'furnishings': Trait(id='furnishings', display="Eyebrow and muzzle furnishings", image='CoatType_Furnishings_furnished.svg'),
    'no_furnishings': Trait(id='no_furnishings', display="No eyebrow and muzzle furnishings", image='CoatType_Furnishings_unfurnished.svg'),

    'low_shedding': Trait(id='low_shedding', display="Low shedding", image='CoatType_Shedding_low.svg'),
    'normal_shedding': Trait(id='normal_shedding', display="Normal shedding", image='CoatType_Shedding_normal.svg'),

    'natural_bob_tail': Trait(id='natural_bob_tail', display="Natural bob tail", image='SpecialFeatures_Tail_bob_tail.svg'),
    'normal_tail': Trait(id='normal_tail', display="Normal length tail", image='SpecialFeatures_Tail_normal_tail.svg'),

    'short_legs': Trait(id='short_legs', display="Shortened leg length", image='SpecialFeatures_Limbs_short.svg'),
    'short_legs_marker': Trait(id='short_legs_marker', display="Shortened leg length", image='SpecialFeatures_Limbs_short.svg'),
    'normal_legs': Trait(id='normal_legs', display="Normal leg length", image='SpecialFeatures_Limbs_normal.svg'),

    'hypoxia_adapted': Trait(id='hypoxia_adapted', display="Adaptable to high altitudes", image='SpecialFeatures_Altitudes_adapted.svg'),
    'not_hypoxia_adapted': Trait(id='not_hypoxia_adapted', display="No adaptation to high altitudes", image='SpecialFeatures_Altitudes_not_adapted.svg'),
}

class TraitCaller():
    def __init__(self, variants, vcf):
        self.variants = variants
        self.vcf = vcf

        self.mask = np.logical_and(
                  variants['CHR'].to_numpy()[:, None] == vcf['variants/CHR'][None, :],
                  variants['POS'].to_numpy()[:, None] == vcf['variants/POS'][None, :])

        self.variant_counts = pd.DataFrame(
            data=np.nan,
            columns=pd.Index(vcf['samples'], name="sample"),
            index=variants.index,
            dtype=float,
        )
        self.gts = pd.DataFrame(
            data='NC',
            columns=pd.MultiIndex.from_product(
                (vcf['samples'], ['A1', 'A2']),
                names=("sample", "GT")),
            index=variants.index,
            dtype=str,
        )
        self.probs = pd.DataFrame(
            data=np.nan,
            columns=pd.Index(vcf['samples'], name="sample"),
            index=variants.index,
            dtype=float,
        )

        # for each variant of interest, need to count the number of matching
        # alleles to variantAllele in vcf.  Keep highest matching number if
        # there are multiple
        for i in self.variant_counts.index:
            if np.sum(self.mask[i, :]) == 0:
                continue
            gts = np.where(
                self.vcf['calldata/GT'][self.mask[i, :]] == 0,
                self.vcf['variants/REF'][self.mask[i, :]][:, None, None],
                self.vcf['variants/ALT'][self.mask[i, :]][:, 0][:, None, None],
            )
            matches = np.sum(gts == self.variants.loc[i, 'variantAllele'], axis=2)
            # taking the max uses the site with the most variant alleles
            self.variant_counts.loc[i, :] = np.amax(matches, axis=0)
            # record the alleles of that site
            self.gts.loc[i, :] = gts[
                np.argmax(matches, axis=0),  # site with most variants
                np.arange(gts.shape[1]),  # each element in dim 1
                np.arange(gts.shape[2])[:, None]  # broadcast to trigger advanced indexing
            ].flatten('F')  # flatten in fortran order to map properly

            # record probabilities
            probs = self.vcf['probability_averages'][self.mask[i, :]]
            self.probs.loc[i, :] = probs[
                np.argmax(matches, axis=0),  # site with most variants
                np.arange(probs.shape[1]),  # each element in dim 1
            ]

    def trait_black(self):
        modules = pd.DataFrame(
            {
                'liver': self._count('mod_liver') > 1,
                'cocoa': self._count('mod_cocoa') > 1,
                'dilute': self._count('mod_dilute') > 1,
                'dom_black': self._count('mod_dom_black') > 0,
                'rec_black': self._count('mod_rec_black') > 0,
                'rec_red': self._stepped(default='Y', steps={1: 'C', 0: 'N'}, module='mod_rec_red'),
            },
            index=self.variant_counts.columns,
        )

        return self.decode_black_trait(modules)

    @classmethod
    def decode_black_trait(cls, modules):
        dilute = np.where(modules['dilute'], 'dilute_', '')

        color = np.where(modules['liver'], 'liver', 'black')
        color = np.where(modules['cocoa'], 'cocoa', color)

        variant = np.where(modules['dom_black'] | modules['rec_black'], '_solid', '')
        variant = np.where(modules['rec_red'] == 'Y', '_leathers', variant)
        return pd.DataFrame.from_records(
            (TRAITS[''.join(id)].to_dict() for id in zip(dilute, color, variant)),
            index=modules.index,
        )


    def trait_red(self):
        # get red intensity modules
        red_variants = self.variant_counts[self.variants['module'] == 'mod_red_intensity'].copy(deep=True)
        # add in variant
        red_variants['variant'] = self.variants.loc[
                self.variants['module'] == 'mod_red_intensity', 'variant']

        # columns are genotype, index is the variant, values are the coefficients
        red_intensity_coefficients = pd.DataFrame(
            {
                0: 0,
                1: [0.472, 0.057, 0.234, 0.700, 0.199],
                2: [1.068, 0.208, 0.208, 1.232, 0.222],
            },
            index=red_variants['variant'],
        )

        # convert to long form for indexing into coefficients
        red_variants = pd.wide_to_long(red_variants, '', 'variant', 'dog')
        # replace nans with 0
        red_variants[red_variants.isna()] = 0
        red_variants = red_variants.astype(int)

        # this looks up the coefficient value for each dog/variant
        idx, cols = pd.factorize(red_variants.index.get_level_values('variant'))
        red_variants['coeff'] = red_intensity_coefficients.to_numpy()[idx, red_variants['']]

        intensity = red_variants['coeff'].groupby('dog').sum() - 1.504

        # convert to trait, this is the opposite order of the if/else statements
        color = np.where(intensity < 0.75, 'tan', 'red')
        color = np.where(intensity < 0, 'yellow', color)
        color = np.where(intensity < -0.75, 'cream', color)
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict() for id in color),
            index=self.variant_counts.columns,
        )

    def trait_agouti(self):
        # num modules x num dogs
        modules = pd.DataFrame(
            {
                'sable': self._count('mod_sable') > 0,
                'tan_points': self._count('mod_tan_points') > 1,
                'rec_black': self._count('mod_rec_black') > 1,
                'dom_black': self._count('mod_dom_black') > 0,
                'rec_red': self._stepped(default='Y', steps={1: 'C', 0: 'N'}, module='mod_rec_red'),
            },
            index=self.variant_counts.columns,
        )

        return self.decode_agouti_trait(modules)

    @classmethod
    def decode_agouti_trait(cls, modules):
        variant = np.where(modules['rec_black'], 'rec_black', 'agouti')
        variant = np.where(modules['tan_points'], 'tan_points', variant)
        variant = np.where(modules['sable'], 'sable', variant)

        hidden =np.where(
            (modules['rec_red'] == 'Y') | modules['dom_black'],
            '_hidden', '')

        return pd.DataFrame.from_records(
            (TRAITS[''.join(id)].to_dict() for id in zip(variant, hidden)),
            index=modules.index,
        )


    def trait_extension(self):
        # num modules x num dogs
        modules = pd.DataFrame(
            {
                'rec_black': self._count('mod_rec_black') > 1,
                'dom_black': self._count('mod_dom_black') > 0,
                'mask': self._count('mod_mask') > 0,
                'grizzle': self._count('mod_grizzle'),
                'domino': self._count('mod_domino'),
                'rec_red': self._count('mod_rec_red'),
            },
            index=self.variant_counts.columns,
        )

        return self.decode_extension_trait(modules)

    @classmethod
    def decode_extension_trait(cls, modules):

        mask = np.where((modules['domino'] == 1) & (modules['rec_red'] == 1),
                        'domino', 'normal_extension')
        mask = np.where((modules['grizzle'] == 1) & (modules['rec_red'] == 1),
                        'grizzle', mask)
        mask = np.where(modules['rec_red'] == 2, 'rec_red', mask)
        mask = np.where(modules['domino'] == 2, 'domino', mask)
        mask = np.where(modules['grizzle'] == 2, 'grizzle', mask)
        mask = np.where((modules['grizzle'] == 1) & (modules['domino'] == 1),
                        'grizzle_domino', mask) 
        mask = np.where(modules['mask'], 'mask', mask)

        hidden = np.where(modules['rec_black'] | modules['dom_black']
                          , '_hidden', '')

        return pd.DataFrame.from_records(
            (TRAITS[''.join(id)].to_dict() for id in zip(mask, hidden)),
            index=modules.index,
        )

    def trait_ticking(self):
        # num modules x num dogs
        ticking = self._count('mod_ticking') > 1

        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in np.where(ticking, 'ticking', 'no_ticking')),
            index=self.variant_counts.columns,
        )

    def trait_brindle(self):
        modules = pd.DataFrame(
            {
                'brindle': self._count('mod_brindle') > 1,
                'dom_black': self._count('mod_dom_black') > 0,
                'rec_red': self._count('mod_rec_red') > 1,
            },
            index=self.variant_counts.columns,
        )
        return self.decode_brindle_trait(modules)

    @classmethod
    def decode_brindle_trait(cls, modules):
        brindle = np.where(modules['brindle'] & modules['rec_red'],
                           'brindle_hidden',
                           'not_brindle')
        brindle = np.where(modules['brindle'] & ~modules['rec_red'] & ~modules['dom_black'],
                           'brindle',
                           brindle)
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in brindle),
            index=modules.index,
        )

    def trait_merle(self):
        modules = pd.DataFrame(
            {
                'merle': False,  # not in variants of interest
                'harlequin': self._stepped(default="N", steps={1: "Y", 2: "EL"}, module='mod_harlequin'),
            },
            index=self.variant_counts.columns,
        )
        return self.decode_merle_trait(modules)

    @classmethod
    def decode_merle_trait(cls, modules):

        merle = np.where(modules['harlequin'] == 'EL',
                           'embryonic_lethal',
                           'uncalled')
        merle = np.where(~modules['merle'] & (modules['harlequin'] == 'N'),
                           'not_merle',
                           merle)
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in merle),
            index=modules.index,
        )

    def trait_coat_length(self):
        modules = pd.DataFrame(
            {
                'coat_length': self._count('mod_coat_length') > 0,
            },
            index=self.variant_counts.columns,
        )
        length = np.where(modules['coat_length'], 'long_coat', 'short_coat')
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in length),
            index=modules.index,
        )

    def trait_coat_texture(self):
        modules = pd.DataFrame(
            {
                'coat_texture': self._stepped(default='curly_coat', 
                                              steps={
                                                  0: 'straight_coat',
                                                  1: 'wavy_coat',
                                              },
                                              module='mod_curly_coat'),
            },
            index=self.variant_counts.columns,
        )
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in modules['coat_texture']),
            index=modules.index,
        )

    def trait_furnishings(self):
        modules = pd.DataFrame(
            {
                'furnishings': self._count('mod_furnishings') > 0,
            },
            index=self.variant_counts.columns,
        )
        furnishings = np.where(modules['furnishings'], 'furnishings', 'no_furnishings')

        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in furnishings),
            index=modules.index,
        )

    def trait_shedding_propensity(self):
        modules = pd.DataFrame(
            {
                'shedding': self._count('Shedding') > 0,
                'curly': self._count('mod_curly_coat') > 0,
            },
            index=self.variant_counts.columns,
        )
        shedding = np.where(modules['shedding'] | modules['curly'], 'low_shedding', 'normal_shedding')
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in shedding),
            index=modules.index,
        )

    def trait_natural_bob_tail(self):
        modules = pd.DataFrame(
            {
                'bobtail': self._stepped(default='uncalled', 
                                              steps={
                                                  0: 'normal_tail',
                                                  1: 'natural_bob_tail',
                                                  2: 'embryonic_lethal',
                                              },
                                              module='BobTail'),
            },
            index=self.variant_counts.columns,
        )

        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in modules['bobtail']),
            index=modules.index,
        )

    def trait_leg_length(self):
        modules = pd.DataFrame(
            {
                'insertion': self._count("fgf4c18", column='locus') > 2,
                'marker': self._count("fgf4c18gwas", column='locus') > 0,
            },
            index=self.variant_counts.columns,
        )

        length = np.where(modules['marker'], 'short_legs_marker', 'normal_legs')
        length = np.where(modules['insertion'], 'short_legs', length)
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in length),
            index=modules.index,
        )

    def trait_high_altitude_adaptation(self):
        # this is specific to this trait as it counts the number of variant *sites*
        # instead of just the number of variants
        counts = np.sum(self.variant_counts[self.variants['module'] == 'AltitudeAdaptation'] > 0, axis=0)
        # counts == 4 indicates all sites with altitude module have a variant
        adapted = np.where(counts == 4, 'hypoxia_adapted', 'not_hypoxia_adapted')
        return pd.DataFrame.from_records(
            (TRAITS[id].to_dict()
                for id in adapted),
            index=self.variant_counts.columns,
        )

    def _count(self, module, column="module"):
        '''count alternate variants
        '''
        return np.sum(self.variant_counts[self.variants[column] == module], axis=0)

    def _stepped(self, default, steps, module):
        '''stepped phenotype
            Will return default unless the sum of recessive counts matches a
            key in the steps dict[int, str]
        '''

        counts = np.sum(self.variant_counts[self.variants['module'] == module], axis=0)
        result = np.where(counts == -1, '', default)
        for count, value in steps.items():
            result = np.where(counts == count, value, result)
        return result

    def predict_phenotypes(self):
        phenotypes = {
            "CoatColor": {
                "black": self.trait_black(),
                "red": self.trait_red(),
            },

            "CoatPattern": {
                "agouti": self.trait_agouti(),
                "extension": self.trait_extension(),
                "ticking": self.trait_ticking(),
                "brindle": self.trait_brindle(),
                "merle": self.trait_merle(),
            },


            "CoatType": {
                "Coat Texture": self.trait_coat_texture(),
                "Coat Length": self.trait_coat_length(),
                "Coat Furnishings": self.trait_furnishings(),
                "Shedding Propensity": self.trait_shedding_propensity(),
            },


            "SpecialFeatures": {
                "Skeletal - Tail Length": self.trait_natural_bob_tail(),
                "Skeletal - Leg Length": self.trait_leg_length(),
                "High Altitude Adaptation": self.trait_high_altitude_adaptation(),
            },

        }
        result = pd.concat(
            (
                results.reset_index().assign(trait=trait, tab=tab)
                for tab, traits in phenotypes.items()
                for trait, results in traits.items()
            ),
            ignore_index=True)
        result = result['sample,tab,trait,result,string,value,color,image'.split(',')]
        return result

    def genotypes(self):
        result = self.variants[['CHR', 'POS', 'tab', 'trait', 'locus', 'variant']].copy()
        # convert wide table to long
        probs = pd.melt(self.probs, ignore_index=False).rename(columns={'value': 'CONF'})
        # convert wide table to long, harder with multiindex columns
        gts = self.gts.unstack().unstack(level=1).reset_index(level=0, drop=False)
        # merge the 3 tables
        result = result.join(gts.merge(probs, on=['variant_index', 'sample']))
        return result['sample,CHR,POS,tab,trait,locus,variant,A1,A2,CONF'.split(',')]


def read_variant_table(args):
    # load variants of interest
    variant_table = pd.read_csv(
        args.variants,
        header=0,
        na_values="NA",
        dtype={'CHR': int, 'POS': int}
    )

    # select imputation and reference
    variant_table = variant_table[
        (variant_table['imputation'] == args.imputation) &
        (variant_table['reference'] == args.reference)
    ]
    variant_table.reset_index(inplace = True)
    variant_table.index.name = "variant_index"

    return variant_table

def read_vcf(vcf):
    try:
        result = allel.read_vcf(
            vcf,
            rename_fields = {'variants/CHROM': 'variants/CHR'},
        )
    except RuntimeError:
        raise ValueError("Update the index file!")

    # remove chr from start of chromosome
    match = re.compile(r'^chr')
    result['variants/CHR'] = np.fromiter(
        (match.sub('', str(chr)) for chr in result['variants/CHR']),
        dtype=int,
    )

    # if no GP data found in vcf, replace with all ones
    if 'calldata/GP' not in result:
        result['calldata/GP'] = np.ones(
            (result['variants/CHR'].size, result['samples'].size),
            dtype=float,
        )
        result['probability_averages'] = np.ones_like(result['calldata/GP'])
    else:
        raise ValueError("need to handle when GP is present")
        result['probability_averages'] = np.mean(np.amax(result['calldata/GP'], axis=2), axis=0)

    return result

def main(args):
    "Extracts genotypes for variants of interest and calls interpretation modules to output phenotype predictions."

    variant_table = read_variant_table(args)

    caller = TraitCaller(
        variants=variant_table,
        vcf=read_vcf(args.vcf),
    )

    phenotypes = caller.predict_phenotypes()
    phenotypes.to_csv(
        f'{args.output}_phenotypeTable.csv',
        sep=",", na_rep="NA", index=False)
    phenotypes.to_csv(
        f'{args.output}_trailblazerPhenotypeTable.tsv',
        columns=['sample','tab','result','string'],
        sep="\t", na_rep="NA", index=False)

    genotypes = caller.genotypes()
    genotypes.to_csv(
        f'{args.output}_genotypeTable.csv',
        sep=",", na_rep="NA", index=False)

    # Darwin's Ark: Trailblazers
    trailblazer_table = genotypes.merge(
        variant_table,
        on=['tab', 'trait', 'locus', 'variant'])[
        ['sample', 'tab', 'trait', 'title', 'gene', 'normalAllele', 'variantAllele', 'A1', 'A2', 'desc']
    ]
    trailblazer_table = trailblazer_table.rename(columns={
                 'sample': 'Sample ID',
                 'title': 'Trait',
                 'gene': 'Gene',
                 'normalAllele': 'Normal Version',
                 'variantAllele': 'Variant Version',
                 'A1': 'First Copy',
                 'A2': 'Second Copy',
                 'desc':'Description',
             })
    trailblazer_table.to_csv(
        f'{args.output}_trailblazerGenotypeTable.csv',
        sep = ",", na_rep = "NA", index = False)

    jsonTable = genotypes.merge(
        variant_table,
        on=['tab', 'trait', 'locus', 'variant'])[
        ['sample', 'tab', 'title', 'gene', 'normalAllele', 'variantAllele', 'A1', 'A2']
    ].dropna().rename(columns={'title': 'name', 'A1': 'firstCopy', 'A2': 'secondCopy'})
    jsonTable['possibleAlleles'] = jsonTable['normalAllele'].astype(str) + " & " + jsonTable['variantAllele']
    jsonTable = jsonTable[(jsonTable['firstCopy'] != "NA") & (jsonTable['secondCopy'] != "NA")]
    jsonTable = jsonTable.drop(columns=['normalAllele','variantAllele'])
    jsonTable['effect'] = 0
    jsonTable.to_csv(
        f'{args.output}_jsonTable.csv',
        sep=",", na_rep="null", index=False)

if __name__ == '__main__':
    args = parse_args()
    main(args)

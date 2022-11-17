from collections import defaultdict
from pprint import pprint
import gzip
import itertools

def main():
    # map trait to dependent modules
    trait_mapping = {
        "black": ["mod_liver", "mod_dilute", "mod_dom_black", "mod_rec_red", "mod_rec_black", "mod_cocoa"],
        "red": ["mod_red_intensity"],
        "agouti": ["mod_sable", "mod_tan_points", "mod_rec_black", "mod_dom_black", "mod_rec_red"],
        "extension": ["mod_mask", "mod_grizzle", "mod_domino", "mod_rec_black", "mod_dom_black", "mod_rec_red"],
        "ticking": ["mod_ticking"],
        "brindle": ["mod_brindle", "mod_rec_red", "mod_dom_black"],
        "merle": ["mod_harlequin"],
        "texture": ["mod_curly_coat"],
        "length": ["mod_coat_length"],
        "furnishings": ["mod_furnishings"],
        "shedding": ["Shedding"],
        "tail": ["BobTail"],
        "legs": ["mod_short_legs"],
        "altitude": ["AltitudeAdaptation"],
    }

    # get chrom/pos of each module
    variants = defaultdict(list)
    for line in open('vars.csv'):
        if line.startswith('#'):
            continue
        tokens = line.strip().split(',')
        variants[tokens[0]].append(tuple(tokens[1:]))

    # read in vcf, removing samples
    preamble = []
    header = ''
    base_vcf = []
    for line in gzip.open('testing/DogAgingProject_2022-08-15_gp-0.7_trait-predictions.vcf.gz', 'rt'):
        line = line.strip()
        if line.startswith('##'):
            preamble.append(line)
            continue
        if line.startswith('#'):
            header = '\t'.join(line.split()[:9])
            continue
        base_vcf.append(line.split()[:9])

    for trait, modules in trait_mapping.items():
        print(trait)
        locations = {
            variant
            for module in modules
            for variant in variants[module]
        }
        # find how many locations in vcf are present
        sites = sum(tuple(line[:2]) in locations for line in base_vcf)
        individuals = 3 ** sites + 1  # each combination of 0/0, 0/1, 1/1 with one all ./.
        blank_line = '\t'.join('0/0' for _ in range(individuals))

        # note this mapping is odd so it isn't always alphabetical and is a string
        # R -> 0/0 (ref)
        # H -> 0/1 (het)
        # A -> 1/1 (alt)
        indiv_names = '\t'.join(''.join(name) for name in itertools.product('RHA', repeat=sites))
        indiv_names += '\t01001'  # ./.

        genotypes = zip(
            *itertools.product(['0/0', '0/1', '1/1'], repeat=sites), 
            itertools.repeat('./.'),
        )

        with gzip.open(f'testing/{trait}.vcf.gz', 'wt') as outfile:
            # write preamble
            outfile.write('\n'.join(preamble))

            # write header
            outfile.write(f'\n{header}\t{indiv_names}\n')

            for line in base_vcf:
                if tuple(line[:2]) in locations:
                    outfile.write("\t".join(line) + '\t' + "\t".join(next(genotypes)))
                else:
                    outfile.write("\t".join(line) + '\t' + blank_line)
                outfile.write("\n")



if __name__ == "__main__":
    main()
    print('done')

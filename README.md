## Repository

### dog-aging-project-dna-curation
Curation of genome sequencing data for the Dog Aging Project

---

## Project

### Dog Aging Project

The [Dog Aging Project](https://dogagingproject.org/) is a nationwide, longitudinal effort to study all dynamics of aging and health span in tens of thousands of companion dogs. As an open, community science initiative, dog-owning participants are an integral part of the project, and all the data that we collect will be shared with scientists around the world. As part of this initiative, we will perform low-coverage DNA sequencing in 10,000 dogs and launch analyses that delineate the genetic and environmental components of aging, longevity, health disorders, and age-related cognitive decline.

### Project 2: Genetics of aging and longevity-related traits in the domesticated dog
Longevity is a complex phenotype governed by genetic and environmental factors. Although considerable progress has been made in understanding the molecular mechanisms of aging, primarily in model organisms, novel interdisciplinary approaches are needed to define the genetic architecture of aging and aging-related diseases, the interaction between genetics and environment, and how these factors contribute to longevity. The domesticated dog is a genetically complex natural model that shares much of its environment with humans.  The close relationships between pet dogs and their owners has enabled the Dog Aging Project to assemble large, well-phenotyped cohorts through community science, and these cohorts continue to grow. The unusual population structure of pet dogs in the United States, with both purebred (ancestry from just one breed) and mixed breed dogs, makes it possible to investigate the interplay of genetics and environment in cohorts considerably smaller than required in human populations.

We propose to leverage state-of-the-art genomics and statistical tools to comprehensively dissect the genetic and environmental determinants of longevity and age-related phenotypes in dogs. By integrating environmental and molecular phenotypes with the largest publicly available genomic resource ever assembled for any natural mammalian model, we will investigate the mechanisms shaping aging-related phenotypes, and how these mechanisms affect longevity and healthspan. Specifically, we will:

(1) Generate high-resolution genetic data that will transform dogs into a natural model for research into the genomics of aging and aging-related diseases. We will produce and process sequencing data for 2000 dogs per year, for a total of 20,000 dogs by 2028, each with detailed phenotype data for aging-related traits and diseases. We will rapidly share genomic data, aging-related phenotype data, and analysis pipelines as an open data resource through the cloud-native Terra platform.

(2) Comprehensively investigate the genetic architecture of aging and aging-related diseases. We will estimate the heritability of aging-related diseases, mobility, frailty, and healthspan, and, through genome-wide association, identify the genomic loci and cellular pathways involved. We will also investigate the genetic architecture of exceptional longevity and identify protective variants associated with long life in our cohort of canine “centenarians”. We will also identify gene-by-environment interactions using novel statistical methods.

(3) Investigate the function of aging-related variants through intermediary quantitative trait loci. We will map DNA methylation QTLs, plasma metabolome QTLs and microbiome QTLs in dogs, intersect these with genome-wide association results, and to explore the functional effects of associated variants and the tissue types in which they might be acting.

---

## Overview

This repository contains and maintains data, code, and workflows needed to curate the low-pass, whole-genome sequencing data generated from saliva samples of dogs enrolled in the scientific cohorts of the Dog Aging Project.

### Sampling

Following cohort enrollment, a saliva sample kit is assigned to the dog and mailed to the participant's address by GBF, Inc., High Point, NC. The participant samples the dog using a DNA Genotek Performagene swab, and returns the sample kit to the mail.

### Library preparation and sequencing

The sample kit arrives at Neogen Corporation - GeneSeek Operations (Lincoln, NE) for library preparation. Prepared libraries undergo short read, low-coverage Illumina DNA sequencing at the Neogen platform.

### Alignment and imputation

Raw sequencing read data is imported by Illumina BaseSpace to the Gencove Platform for alignment and statistical genotype imputation using software *loimpute* by Gencove, Inc. (New York, NY). Not all sequencing runs reach the quality assurance and control metrics needed for genotype imputation. Sequencing reads are aligned to the CanFam3.1 reference genome assembly (NCBI accession GCF_000000145.2). Genotype imputation is performed using  *loimpute* and the dog low-pass v2.0 [0.1x-6x] panel of reference haplotypes, consisting of 34,191,821 single nucleotide polymorphisms and 11,943,064 insertions / deletions and representing 540 dogs of known breed ancestry distributed among ~133 breeds, 28 dogs of mixed breed ancestry, 12 dogs of unknown ancestry, 62 worldwide indigenous or village dogs, 33 wolves, and 1 coyote.

### Curation and assignment

Successful data deliverables migrate to a Dog Aging Project workspace on Terra.bio for curation and assignment to study IDs. For sequencing runs that do not reach quality needed for genotype imputation, only the raw sequencing reads (FASTQ files) are migrated. In the case of multiple samples per study ID, only genotyping data from the sequencing run with the higher effective genome coverage are assigned to the dog.

---

## Glossary

**REDCap** - ???

**Portal** - the user interface for a participant's dog on the Dog Aging Project website

**Study ID** - unique identifier for a dog in REDCap

**Swab ID** - unique identifier for a saliva sample aka barcode


---

## Workspace

### Data models

The Terra data model tables uniquely identify and provide information on each entity of the genome sequencing process, from dog to saliva sample to sequencing run and resulting data.

While many columns are set by outputs of workflows, many have been maintained manually (indicated in italics).

#### Dogs / Study IDs: `participant` table

Columns:

- `entity:participant_id`: dog's Study ID
- `cohort`: _name of scientific cohort_
- `cohort_code`: REDCap ID of scientific cohort
- `report_flag_issue`: _logical, indicates if Genomic Report should not be pushed to dog's Portal_
- `report_flag_pushed**`: _logical, indicates if Genomic Report JSON received on REDCap_

#### Samples / Swab IDs: `sample` table

Columns:

- `entity:sample_id`: dog's Swab ID

#### Runs / Platform IDs: `platform` table

Columns:

- `entity:platform_id`: unique ID for sequencing run

---

## Release

### Raw sequencing data

Raw sequencing read data (FASTQ files) are deposited to the NCBI Sequence Read Archive under BioProject PRJNA800779. Sample IDs correspond to swab IDs (`swab_id`), not study IDs (`dog_id`).

### Genotyping data

Genotyping data (imputed variant calls) are released as a *PLINK1* bfile set. Sample IDs correspond to study IDs (`dog_id`).

### Genomic reports

The genomic reports are data derived and inferred from genetic data that are neither raw sequencing nor genotyping data. Sample IDs correspond to study IDs (`dog_id`).

---

## Methods

### Sex confirmation

We inferred the sex of each dog from aligned sequencing read data (BAM files) using *SAMtools* to measure the ratio of X chromosome coverage to autosomal coverage. Ratios over or equal to 0.7 were inferred as female, under 0.7 inferred as male. We compared genetically-inferred sex to owner-reported sex given in the Health and Life Experience Survey (HLES).

### Filtering

Each variant call file (VCF) was filtered for calls of genotype probability above 70% (max(GP) > 0.7) using *BCFtools*. Each dog had on average XXX ± SD:XXX single nucleotide polymorphisms (SNPs) and XXX ± SD:XXX insertions/deletions called at this GP threshold. Sample IDs were converted from swab IDs to study IDs and genotyping data from multiple samples were merged in batches using *BCFtools*. Each merged VCF was converted to a *PLINK2* pfile data set (.pgen / .pvar / .psam) containing all variant calls (insertions, deletions, single nucleotide polymorphisms, multiallelic variants). Then, only biallelic SNPs were selected (XXX SNPs total) and converted to a *PLINK1* bfile data set (.bed / .bim / .fam). After filtering for a minimum minor allele frequency of 1% and minimum genotype rate of 95%, the final data set included XXX SNPs. Confirmed sexes were encoded into this final *PLINK1* data set.

### Coefficients of inbreeding

We estimated coefficients of inbreeding as autozygosity, or the proportion of genome covered by runs of homozygosity (ROH). We scanned for ROH using *PLINK1* across the XXX unfiltered biallelic SNP genotypes of genotype probability >70% with the following settings: minimum run length of 500kb (`--homozyg-kb 500`) and minimum SNP count of 100 SNPs (`--homozyg-snp 100`), at a density of 1kb per SNP (`--homozyg-density 1`), with no two SNPs more than 500kb apart (`--homozyg-gap 500`), and only 1 heterozygous genotype tolerated per window (`--homozyg-window-het 1`), performing scans without LD-based pruning on chromosomes 1-38 (`--chr 1-38`). All other settings were *PLINK1* defaults. We then calculated the autosomal ROH-estimated coefficient of inbreeding (F<sub>ROH</sub>, or CoI) from the total ROH segment length divided by the total SNP-covered autosomal length (2,203,765,000 bases) used for ROH detection.

### Genomic size prediction

We applied a random forest model developed on the Darwin's Ark (darwinsark.org) data set of 1,730 dogs with surveyed height phenotypes, coded as 0 = tiny, 1 = small, 2 = medium, 3 = large, or 4 = giant as assessed by dog owners, and 2,733 SNPs associated with canine height to predict body sizes (mean square error = 0.3).

### Global ancestry

We selected publicly available genotype data from 109 modern breeds with at least 4 dogs per breed, 3 regional village dog populations (4 Nigerian village dogs, 5 Vietnamese village dogs, 55 Chinese village dogs), and 2 wolf populations (19 North American wolves and 25 Eurasian wolves) (see *Populations* for full list of population labels and counts). We used *PLINK2* to identify ancestry-informative markers. We selected 4,267,732 biallelic single nucleotide polymorphisms with <10% missing genotypes, and calculated the Wright’s F-statistics using Hudson method for (1) each dog breed versus all other purebred dogs; (2) all village dogs versus all other purebred dogs; (3) each regional village dog population; (4) all wolves versus all other dogs; (5) North American wolves versus Eurasian wolves. We selected 1,569,037 SNPs with F<sub>ST</sub> > 0.5 across all comparisons, and performed LD-based pruning in 250kb windows for r<sup>2</sup> > 0.2 to extract 115,427 markers for global ancestry inference.

We merged genotype data for these biallelic SNPs from query samples with genotype data from reference samples, then performed global ancestry inference using *ADMIXTURE* in supervised mode (random seed: 43) to infer ancestry from these populations. We report only admixture proportions over 1% for each dog.

### Genetic relationship matrix

We calculated the variance-standardized relationship matrix generated using *PLINK2* (`--make-rel`) across all autosomes.

### Pairwise kinship

We applied the *PLINK2* implementation of the KING-robust estimator to measure kinship *k* between pairs of dogs from mixed populations. Relationships were either unrelated, and removed from the resulting table, or labeled as related (*k* > 0), second-degree (*k* >= 0.125), or first-degree (*k* >= 0.25). We apply a cutoff *k* >= 0.35 to identify duplicate samples. No data from duplicate samples are included in this release.

---
### Software

The following versions of software were used when referenced:

- *SAMtools*: SAMtools v1.11
- *BCFtools*: BCFtools v1.8
- *PLINK1*: PLINK v.1.9 Stable (19 Oct 2020)
- *PLINK2*: PLINK v.2.00 Development (02 Mar 2021)
- *ADMIXTURE*: ADMIXTURE v.1.3.0

### Populations

The full list of populations and sample sizes used for ancestry inference are as follows:

**Dog Breeds**

- `afghan_hound` (13)
- `airedale_terrier` (17)
- `akita` (14)
- `alaskan_malamute` (16)
- `american_cocker_spaniel` (21)
- `american_pit_bull_terrier` (30)
- `australian_cattle_dog` (19)
- `australian_shepherd` (32)
- `basenji` (15)
- `basset_hound` (17)
- `beagle` (30)
- `bearded_collie` (13)
- `belgian_groenendael` (7)
- `belgian_malinois` (15)
- `belgian_tervuren` (23)
- `berger_picard` (5)
- `bernese_mountain_dog` (24)
- `bichon_frise` (15)
- `bloodhound` (16)
- `border_collie` (31)
- `border_terrier` (17)
- `borzoi` (13)
- `boston_terrier` (18)
- `bouvier_des_flandres` (4)
- `boxer` (29)
- `brittany` (16)
- `bull_terrier` (13)
- `bullmastiff` (23)
- `cairn_terrier` (14)
- `cavalier_king_charles_spaniel` (20)
- `chesapeake_bay_retriever` (16)
- `chihuahua` (20)
- `chinese_crested` (17)
- `chinook` (14)
- `chow_chow` (15)
- `collie` (15)
- `dachshund` (32)
- `dalmatian` (16)
- `doberman_pinscher` (22)
- `english_bulldog` (20)
- `english_cocker_spaniel` (14)
- `english_setter` (15)
- `english_shepherd` (16)
- `english_springer_spaniel` (16)
- `entlebucher` (13)
- `finnish_spitz` (13)
- `french_bulldog` (17)
- `german_shepherd_dog` (33)
- `german_shorthaired_pointer` (17)
- `golden_retriever` (88)
- `gordon_setter` (14)
- `great_dane` (15)
- `great_pyrenees` (14)
- `greater_swiss_mountain_dog` (17)
- `greenland_sled_dog` (11)
- `greyhound` (30)
- `havanese` (18)
- `irish_setter` (15)
- `irish_wolfhound` (21)
- `italian_greyhound` (13)
- `jack_russell_terrier` (22)
- `labrador_retriever` (65)
- `lagotto_romagnolo` (4)
- `leonberger` (51)
- `lhasa_apso` (12)
- `maltese` (20)
- `mastiff` (15)
- `miniature_pinscher` (16)
- `miniature_schnauzer` (30)
- `newfoundland` (17)
- `norfolk_terrier` (12)
- `norwegian_elkhound` (14)
- `norwich_terrier` (4)
- `nova_scotia_duck_tolling_retriever` (15)
- `old_english_sheepdog` (13)
- `papillon` (17)
- `pekingese` (13)
- `pembroke_welsh_corgi` (22)
- `pomeranian` (21)
- `poodle` (29)
- `portuguese_water_dog` (15)
- `pug` (23)
- `rhodesian_ridgeback` (12)
- `rottweiler` (28)
- `saint_bernard` (17)
- `saluki` (12)
- `samoyed` (14)
- `schipperke` (13)
- `scottish_terrier` (15)
- `shar_pei` (12)
- `shetland_sheepdog` (17)
- `shiba_inu` (17)
- `shih_tzu` (18)
- `siberian_husky` (21)
- `sloughi` (4)
- `soft_coated_wheaten_terrier` (15)
- `staffordshire_bull_terrier` (15)
- `standard_schnauzer` (4)
- `tibetan_mastiff` (12)
- `tibetan_spaniel` (12)
- `tibetan_terrier` (17)
- `toy_poodle` (28)
- `vizsla` (12)
- `weimaraner` (15)
- `west_highland_white_terrier` (21)
- `whippet` (14)
- `wire_fox_terrier` (16)
- `wirehaired_pointing_griffon` (13)
- `yorkshire_terrier` (71)

**Regional Village Dogs**

- `village_dog_china` (55)
- `village_dog_nigeria` (4)
- `village_dog_vietnam` (5)

**Wolf Populations**

- `wolf_eurasian` (25)
- `wolf_north_american` (19)

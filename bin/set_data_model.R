library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-g", "--statustable"),
              type="character",
              default= "statusTable.tsv",
              help="status table file name",
              metavar="character"),
  make_option(c("-p", "--platformtable"),
              type="character",
              default= "platformTable.tsv",
              help="platform table file name",
              metavar="character"),
  make_option(c("-k", "--kittable"),
              type="character",
              default= "kitTable.tsv",
              help="kit table file name",
              metavar="character"),
  make_option(c("-s", "--sampletable"),
              type="character",
              default= "sampleTable.tsv",
              help="sample table file name",
              metavar="character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# This script will prepare the tables for updating the data model

# Load statusTable (sequencing run statuses from platform)
statusTable = read_tsv(opt$statustable,
                       col_names = c("datetime",
                                     "entity:platform_id",
                                     "client",
                                     "status",
                                     "availability"))
statusTable$client = as.character(statusTable$client)

# Load platformTable (data model for sequencing runs)
platformTable = read_tsv(opt$platformtable)
platformTable$sample = as.character(platformTable$sample)
platformTable$client = as.character(platformTable$client)
platformTable$participant = as.character(platformTable$participant)

# Load kitTable (sample kit assignments)
kitTable = read_tsv(opt$kittable)
# kitTable$`entity:sample_id` = as.character(kitTable$`entity:sample_id`)
kitTable$`entity:sample_id` = as.character(kitTable$`sample_id`)
kitTable$dog_id = as.character(kitTable$dog_id)
kitTable = kitTable %>% rename(participant = dog_id)

# Load sampleTable (data model for sample kits)
sampleTable = read_tsv(opt$sampletable)
sampleTable$`entity:sample_id` = as.character(sampleTable$`entity:sample_id`)
sampleTable$participant = as.character(sampleTable$participant)

# entity:participant_id
# columns to update: entity:participant_id, cohort_code
participantDataModel = kitTable %>%
  select(`entity:participant_id` = participant,
         cohort_code = cohort) %>%
  arrange(as.numeric(`entity:participant_id`)) %>%
  unique()

write_tsv(x = participantDataModel, file = "participant.tsv", na = '', col_names = TRUE)

# entity:sample_id
# columns to update: entity:sample_id, participant, date_swab_arrival_laboratory, sample_type
sampleDataModel = kitTable %>%
  filter(!`entity:sample_id` %in% sampleTable$`entity:sample_id`) %>%
  select(`entity:sample_id`,
         participant,
         date_swab_arrival_laboratory,
         sample_type)

write_tsv(x = sampleDataModel, file = "sample.tsv", na = '', col_names = TRUE)

# entity:platform_id
# columns to update: entity:platform_id, sample, client, participant, status, availability, datetime
newRuns = anti_join((statusTable %>% select(`entity:platform_id`)),
                    (platformTable %>% select(`entity:platform_id`)),
                    by = 'entity:platform_id') %>% pull(`entity:platform_id`)
colnames(platformTable)

platformDataModel = statusTable %>%
  filter(`entity:platform_id` %in% newRuns) %>%
  merge((kitTable %>% select(`entity:sample_id`,participant)),
        by.x = "client",
        by.y = "entity:sample_id",
        all.x = T) %>%
  mutate(sample = client) %>%
  select(`entity:platform_id`,
         client,
         sample,
         participant,
         status,
         availability,
         datetime)

write_tsv(x = platformDataModel, file = "platform.tsv", na = '', col_names = TRUE)

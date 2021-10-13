library(dplyr)
library(readr)
library(tidyr)

# read in gtdb taxonomy csv.
# format accession to filename, and output. 
# required input for charcoal decontamination

gtdb_taxonomy <- read_csv(snakemake@input[["gtdb_taxonomy"]]) %>%
  select(-is_representative, -signame) %>%
  mutate(ident = paste0(ident, "_genomic.fna.gz"))
write_csv(gtdb_taxonomy, snakemake@output[["gtdb_lineages"]], 
          col_names = F)

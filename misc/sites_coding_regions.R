setwd("/local/workdir/yc2644/CV_NYC_misc/coding_regions")

library(tidyverse)

gff <- read_tsv("genomic.gff.tsv", col_names = F) %>% dplyr::select(-X6,-X8,-X9) %>% 
  dplyr::rename(chr=X1,source=X2,type=X3,start=X4,end=X5, strand=X7) %>% 
  dplyr::filter(chr != "NC_007175.2") # exclude mtDNA

gff_cds <- gff %>% filter(type=="CDS")

gff_cds_bed <- gff_cds %>% dplyr::select(chr, start, end)

write_tsv(gff_cds_bed, "gff_cds.bed", col_names = F)

# proportion of the nuclear genome that are coding regions --------


# Calculate the total length of the genome
total_genome_length <- 684723884 # excluding mtDNA

# Calculate the total length of CDS regions
gff_cds.sorted.merged <- read_tsv("gff_cds.sorted.merged.bed", col_names = F) # excluding mtDNA
coding_length <- sum(gff_cds.sorted.merged$X3 - gff_cds.sorted.merged$X2 + 1)
coding_length

# Calculate the proportion of coding regions
coding_proportion <- coding_length / total_genome_length
coding_proportion

# proportion of the haplotigs-masked nuclear genome that are coding regions --------

# Calculate the total length of haplotigs
haplotigs <- read_tsv("haplotigs.bed", col_names = F) %>% 
  dplyr::select(X1,X2,X3)
haplotigs_length <- sum(haplotigs$X3 - haplotigs$X2 + 1)
haplotigs_length

# Calculate the total length of the genome
total_genome_length <- 684723884 - 106492916 # excluding mtDNA and haplotigs

# Calculate the total length of CDS regions not overlapping with haplotigs
gff_cds.sorted.merged.nohaplotigs <- read_tsv("gff_cds.sorted.merged.nohaplotigs.bed", col_names = F) # excluding mtDNA
coding_length <- sum(gff_cds.sorted.merged.nohaplotigs$X3 - gff_cds.sorted.merged.nohaplotigs$X2 + 1)
coding_length

# Calculate the proportion of coding regions
coding_proportion <- coding_length / total_genome_length
coding_proportion

# proportion of the mappable, haplotigs-masked nuclear genome that are coding regions --------
# for lcWGS
# Calculate the total length of CDS regions that are mappable and not overlapping with haplotigs
CDS_mapped_lcWGS <- read_tsv("CDS_mapped_lcWGS.bed", col_names = F) # excluding mtDNA
CDS_mapped_lcWGS_length <- sum(CDS_mapped_lcWGS$X3 - CDS_mapped_lcWGS$X2 + 1)
CDS_mapped_lcWGS_length 

mapped_lcWGS <- read_tsv("mapped_nonzero_chr_pos.merged.lcWGS.bed", col_names = F)
mapped_lcWGS_length <- sum(mapped_lcWGS$X3 - mapped_lcWGS$X2 + 1)
mapped_lcWGS_length 

CDS_mapped_lcWGS_length/mapped_lcWGS_length

# for ddRAD
# Calculate the total length of CDS regions that are mappable and not overlapping with haplotigs
CDS_mapped_ddRAD <- read_tsv("CDS_mapped_ddRAD.bed", col_names = F) # excluding mtDNA
CDS_mapped_ddRAD_length <- sum(CDS_mapped_ddRAD$X3 - CDS_mapped_ddRAD$X2 + 1)
CDS_mapped_ddRAD_length 

mapped_ddRAD <- read_tsv("mapped_nonzero_chr_pos.merged.ddRAD.bed", col_names = F)
mapped_ddRAD_length <- sum(mapped_ddRAD$X3 - mapped_ddRAD$X2 + 1)
mapped_ddRAD_length 

CDS_mapped_ddRAD_length/mapped_ddRAD_length 


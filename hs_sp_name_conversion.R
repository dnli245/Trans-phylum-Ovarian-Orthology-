library(dplyr)
library(stringr)
library(purrr)
# load reciprocal blast tables
sp_to_hm <- read.table("example_data/maps/hmsp/sp_to_hm.txt")
hm_to_sp <- read.table("example_data/maps/hmsp/hm_to_sp.txt")

# read in human gtf 
library(data.table)
hm_gtf <- fread(
  "human_v47_annotations.gtf",
  sep = "\t",
  header = FALSE,
  quote = "",
  stringsAsFactors = FALSE,
  col.names = c("seqname","source","feature","start","end",
                "score","strand","frame","attributes")
)

# extract gtf info
header <- hm_gtf$"attributes"

gene_id <- gsub('.*gene_id "([^"]+)".*', '\\1', header)
gene_name <- gsub('.*gene_name "([^"]+)".*', '\\1', header)

# edit blast tables
sp_to_hm$V2 <- sapply(
  strsplit(sp_to_hm$V2, "\\|"),
  function(x) if (length(x) >= 2) x[2] else NA
)

hm_to_sp$V1 <- sapply(
  strsplit(hm_to_sp$V1, "\\|"),
  function(x) if (length(x) >= 2) x[2] else NA
)

hm_conversion_table <- data.frame(GeneID = gene_id, GeneName = gene_name, stringsAsFactors = FALSE)
hm_conversion_table <- read.csv("hm_conversion_table.txt", sep="\t", header = FALSE)
sp_conversion_table <- read.csv("sp_conversion_table.txt", sep="\t", header = FALSE)

# find alternate names given reference name

sp_map <- setNames(sp_conversion_table$V3,
                   sp_conversion_table$V2)
hm_map <- setNames(hm_conversion_table$V2,
                   hm_conversion_table$V1)

sp_to_hm_new <- sp_to_hm
sp_to_hm_new$V1 <- sp_map[ sp_to_hm$V1 ]
sp_to_hm_new$V2 <- hm_map[ sp_to_hm$V2 ]
hm_to_sp_new <- hm_to_sp
hm_to_sp_new$V1 <- hm_map[ hm_to_sp$V1 ]
hm_to_sp_new$V2 <- sp_map[ hm_to_sp$V2 ]

write.table(sp_to_hm_new, "sp_to_hm_alt.txt",sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(hm_to_sp_new, "hm_to_sp_alt.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# use alt tables for future SAMap analysis
write.table(hm_conversion_table, "hm_conversion_table.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

hm_conversion_table <- read.table("hm_conversion_table.txt")

# Example 1
# Please read everything!
# This example uses the bugphyzz::makeSignatures function, if it doesn't work for you, try
# example number 2
# Additionally, this example uses cMD data on the species level, which Bugphyzz has "less
# of." Take a look at example number 2, which uses cMD on a genus level

### Load required packages for this project
# Make sure to install these packages
pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage")

### Create GMT File ###

# Step 1: Make signatures
sigs <- bugphyzz::physiologies(c("gram stain", "aerophilicity", "butyrate producing"))
sigs <- lapply(sigs, bugphyzz::makeSignatures) %>% 
      unlist(recursive = FALSE) 

# Write the signatures to a GMT data structure and load it into R environment
EnrichmentBrowser::writeGMT(sigs, gmt.file = "my.gmt")
gmt <- qusage::read.gmt("my.gmt")

### Load and format Data from cMD ###

# Notes: 
# Head to https://github.com/waldronlab/curatedMetagenomicDataCuration/tree/master/inst/curated
# Choose a study that has the right parameters for you (disease, sample size, etc.)
# Remember/save the name of the study, for this example "QinN_2014"

qin <- curatedMetagenomicData("QinN_2014.relative_abundance", dryrun = F, counts = T, rownames = "NCBI")
qin.se <- qin[[1]]

# Differentiate control and experimental groups
# Here, we assign the number 1 to cirrhosis, 0 to the controls
grp1 = ifelse(qin.se$study_condition == "cirrhosis" , 1, 0)
qin.se$GROUP = grp1

### Data analysis portion ###
# Please research these analyses if you are unfamiliar

# Perform DESeq2 on the SummarizedExperiment
qin.se.Ana <- deAna(qin.se, de.method = "DESeq2", filter.by.expr = F)

# Overrepresentation analysis (ORA)
ora_results <- sbea("ora", qin.se.Ana, gmt, perm = 0)
ora_results <- gsRanking(ora_results, signif.only = F) %>%
      as.data.frame()
colnames(ora_results) <- sub("GENE", "MICROBE", colnames(ora_results))

# Gene set enrichment analysis (GSEA)
gsea_results <- sbea("gsea", qin.se.Ana, gmt, perm = 1000)
gsea_results <- gsRanking(gsea_results, signif.only = F) %>%
      as.data.frame()
colnames(gsea_results) <- sub("GENE", "MICROBE", colnames(gsea_results))

gdata::keep(qin.se, ora_results, gsea_results, sure = TRUE)

# Try making your own Volcano Plot for fun

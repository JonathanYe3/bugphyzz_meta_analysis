# Example 3
# Please read everything!
# This example uses the old make_signatures function, which should work for everyone
# IMPORTANT: This example uses all species level data through downstream propagation. 

### Load required packages for this project
# Make sure to install these packages
pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage", "taxize")
source("R/old_make_signatures.R")

### Create GMT File ###

# Step 1: Make signatures and propagate downstream
sigs <- bugphyzz::physiologies(c("gram stain", "aerophilicity")) 
sigs <- lapply(sigs, make_signatures, taxids = "NCBI_ID") %>% 
      unlist(recursive = F) %>% 
      lapply(unique)
children <- sigs %>% 
      lapply(taxizedb::children, db = "ncbi") %>% 
      lapply(data.table::rbindlist) %>% 
      lapply(dplyr::pull, id)
sigs <- list(children, sigs) %>% 
      unlist(recursive=F) 
sigs <- tapply(unlist(sigs, use.names = FALSE), rep(names(sigs), lengths(sigs)), FUN = c)

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

#You don't have to run this next line, it clears the clutter but you may want to keep more
gdata::keep(qin.se, ora_results, gsea_results, sure = TRUE)

# Try making your own Volcano Plot for fun

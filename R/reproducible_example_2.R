# Example 2
# Please read everything!
# This example uses the old make_signatures function, which should work for everyone
# Additionally, this example uses cMD data on the genus level, which Bugphyzz has "more of"

### Load required packages for this project
# Make sure to install these packages
pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage")
source("R/old_make_signatures.R")

### Create GMT File ###

# Step 1: Make signatures
# Feel free to customize the physiologies you want to study
sigs <- bugphyzz::physiologies(c("gram stain", "aerophilicity", "butyrate producing"))
sigs <- lapply(sigs, make_signatures) %>% 
      unlist(recursive = FALSE) 

# Write the signatures to a GMT data structure and load it into R environment
EnrichmentBrowser::writeGMT(sigs, gmt.file = "my.gmt")
gmt <- qusage::read.gmt("my.gmt")

### Load and format Data from cMD ###

# Head to https://github.com/waldronlab/curatedMetagenomicDataCuration/tree/master/inst/curated
# Choose a study that has the right parameters for you (disease, sample size, etc.)
# Remember/save the name of the study, for this example "YeZ_2018"

ye <- curatedMetagenomicData("YeZ_2018.relative_abundance", dryrun = F, counts = T, rownames = "long")
ye.se <- ye[[1]]

# Differentiate control and experimental groups
# Here, we assign the number 1 to BD, 0 to the controls
grp1 = ifelse(ye.se$study_condition == "BD" , 1, 0)
ye.se$GROUP = grp1

### Data analysis portion ###
# Please research these analyses if you are unfamiliar

# Perform DESeq2 on the SummarizedExperiment
ye.se.Ana <- deAna(ye.se, de.method = "DESeq2", filter.by.expr = F)

# Format bacteria names to genus only - they're currently in MetaPhlAn format
rownames(ye.se.Ana) <- vapply(rownames(ye.se.Ana), get_Genus, character(1), USE.NAMES = FALSE)

# Overrepresentation analysis (ORA)
ora_results <- sbea("ora", ye.se.Ana, gmt, perm = 0)
ora_results <- gsRanking(ora_results, signif.only = F) %>%
      as.data.frame()
colnames(ora_results) <- sub("GENE", "MICROBE", colnames(ora_results))

# Gene set enrichment analysis (GSEA)
gsea_results <- sbea("gsea", ye.se.Ana, gmt, perm = 1000)
gsea_results <- gsRanking(gsea_results, signif.only = F) %>%
      as.data.frame()
colnames(gsea_results) <- sub("GENE", "MICROBE", colnames(gsea_results))

# This line of code gets rid of all items in global env besides those listed
gdata::keep(ye.se, ora_results, gsea_results, sure = TRUE)

# View your results
View(gsea_results)
View(ora_results)

# Try making your own Volcano Plot for fun

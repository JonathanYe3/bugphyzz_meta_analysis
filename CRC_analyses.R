# Load packages and functions
pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage", "taxize")
source("R/old_make_signatures.R")
source("R/deAnaPhyzz.R")

# Make GMT - If you've pulled the gmt, no need to run this section
phyzz_dfs <- bugphyzz::physiologies("all") %>% 
      within(rm("habitat", "extreme environment", "plant pathogenicity", 
                "isolation site", "growth medium", "length", "width"
                ))
phyzz_dfs <- lapply(phyzz_dfs, mutate, Attribute_value = str_squish(Attribute_value)) %>% 
      lapply(mutate, Attribute = str_squish(Attribute))

sigs <- phyzz_dfs %>%  
      lapply(make_signatures, taxids = "NCBI_ID") %>% 
      unlist(recursive = F) %>% 
      lapply(unique) %>% 
      purrr::keep( ~ !is.null(.) )

gmt <- makeSpeciesPhyzzGMT(sigs)

# Load GMT (if you haven't already)
gmt <- qusage::read.gmt("my.gmt")

# Differential Analyses
studies <- c("ThomasAM_2019_c", "WirbelJ_2018")
thomas.se.ana <- deAnaPhyzz(study = studies[1], condition = "CRC")
wirbel.se.ana <- deAnaPhyzz(study = studies[2], condition = "CRC")

# ORA and GSEA
results <- lapply(list(thomas.se.ana, wirbel.se.ana), oraGsea)
names(results) <- studies
results <- unlist(results, recursive = F)

#Write to workbook
library(openxlsx)
write.xlsx(results, file = "~bugphyzz_meta_analysis/Data/results.xlsx")

my_list <- list(temp, temp2)
test_lapply <- lapply(my_list, taxizedb::children, db = "ncbi") %>% 
      lapply(data.table::rbindlist) %>% 
      lapply(dplyr::pull, id)


final <- lapply(test_lapply, data.table::rbindlist)
taxid <- lapply(final, dplyr::pull, id)

#samuel moment
taxid <- temp
new_taxa <- taxizedb::children(taxid, db = 'ncbi')[[1]] |>
      dplyr::rename(NCBI_ID = id, Taxon_name = name, Rank = rank) |>
      dplyr::filter(Rank == 'species')
      z <- j |>
      dplyr::select(-Parent_NCBI_ID, -Parent_name, -Parent_rank) |>
      dplyr::rename(
      Parent_NCBI_ID = NCBI_ID, Parent_name = Taxon_name, Parent_rank = Rank
      )
      
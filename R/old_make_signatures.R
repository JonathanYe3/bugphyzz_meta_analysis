make_signatures <- function(df, taxids = "Taxon_name", rank = "all"){
      attribute_names <- unique(df[["Attribute"]])
      dat <- lapply(attribute_names, fetch_bugs, dat = df, taxids = taxids)
      names(dat) <- attribute_names
      return(dat)
}

fetch_bugs <- function(attribute_name, dat, taxids){
      bugs <- dat[dat[["Attribute"]] == attribute_name,] %>% 
            filter(Attribute_value == TRUE | Attribute_value == "TRUE") %>% 
            c()
      
      bugs <- bugs[[taxids]]
      if(length(bugs)==0){
            return (NULL)
      }
      else return(bugs[!is.na(bugs)])
}

get_Genus <- function(MetaPhlAn){
      temp <- sub(".*g__", "", MetaPhlAn)
      temp <- gsub("\\|s__.*", "", temp)
      return(temp)
      #gsub("_noname", "", temp)
}

.getLast <- function(n)
{
      spl <- unlist(strsplit(n, "\\|"))
      spl[length(spl)]
}
return_dgidb <- function(input, type = "genes"){
  require(stringr)
  require(jsonlite)
  require(httr)
  require(stringr)
  
  url = "http://dgidb.genome.wustl.edu"
  path = paste0("/api/v2/interactions.json",
                "?", type,"=",input %>% str_replace_all(., " ", ""))
  
  response <- GET(url = url, path = path)
  
  response$content %>%
    rawToChar() %>% 
    fromJSON() -> temp
  
  do.call(what = "rbind",
          args = lapply(temp, as.data.frame)) -> temp
  
  if(nrow(temp) > 0 & !("suggestions" %in% colnames(temp))){
    as.data.frame(temp$interactions) %>%
      as_tibble() %>%
      select(contains("Name"),  source, interactionType, contains("Id")) %>%
      mutate(input = input,
             alt_input = NA,
             input_type = type) %>%
      return()
  } else if(("suggestions" %in% colnames(temp)) & (temp$suggestions %>% unlist %>% is.null() != TRUE)){
    input.alt <- temp$suggestions %>% unlist %>% .[1]
    cat(paste0(c("corrected ", input, "to", input.alt, "\n")))
    path.alt = paste0("/api/v2/interactions.json",
                      "?", type,"=",input.alt %>% str_replace_all(., " ", ""))
    
    GET(url = url, path = path.alt) %>%
      .$content%>%
      rawToChar() %>% 
      fromJSON() -> temp.alt
    
    do.call(what = "rbind",
            args = lapply(temp.alt, as.data.frame)) -> temp.alt
    
    if(nrow(temp.alt) > 0 & !("suggestions" %in% colnames(temp.alt))){
      as.data.frame(temp.alt$interactions) %>%
        as_tibble() %>%
        select(contains("Name"),  source, interactionType, contains("Id")) %>%
        mutate(input = input,
               alt_input = input.alt,
               input_type = type) %>%
        return()
    }
  } else 
    cat(paste0(c("could not find ", input, "\n")))
}


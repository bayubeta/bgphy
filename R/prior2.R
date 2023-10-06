




setPriors.bgphy_model <- function(model){
  # get parameter names
  parnames <- getParamNames(model)

  parnames_list <- sapply(parnames, function(x) NULL)

  class(parnames_list) <- "bgphy_prior"
  attr(parnames_list, "bgphy_model") <- model

  return(parnames_list)
}

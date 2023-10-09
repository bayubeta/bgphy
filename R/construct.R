
# create a function to define a model for Global/Mixed models
#' @export
setModel <- function(tree, modeltypes, startNodes = NULL){

  # check if elements of modeltypes is in {"OU", "BM"}
  stopifnot("Model name not found! Choose between BM or OU." = all(modeltypes %in% c("BM", "OU")))
  # check if startnodes are not null if the there are regimes

  # check startNodes values (numeric & in the tree)



  # create a PCM object given the tree and specified model and regimes
  r <- length(modeltypes) # number of regimes
  if (r == 1){
    # Global model
    PCMmodel <- PCMBase::PCM(modeltypes)
  }else{
    # Mixed model
    # retrieve the names of regimes
    rnames <- names(modeltypes)
    # code names in PCMBase
    modelStrings <- ifelse(modeltypes == "BM",
                           "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                           "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
    # Create the PCM object
    PCMmodel <- PCMBase::MixedGaussian(k = 1, modelTypes = modelStrings, mapping = stats::setNames(1:r, rnames))

    # update the tree with the regimes information
    tree <- PCMBase::PCMTreeSetPartRegimes(PCMBase::PCMTree(tree),
                                           part.regime = setNames(1:r, startNodes),
                                           setPartition = TRUE, inplace = FALSE)
  }

  model <- structure(list(model = PCMmodel, tree = tree), class = "bgphy_model")
  attr(model, "modeltypes") <- modeltypes

  # set default priors for the parameters
  model$priors <- defaultPriors(model)


  return(model)
}


# print object of class bgphy (the constructed model)
#' @export
print.bgphy_model <- function(model){
  # retrieve regime names
  rnames <- attr(model$model, "regimes") # get regime names
  r <- length(rnames) # get number of regimes

  # retrieve model types
  modeltypes <- attr(model, "modeltypes")

  if (r == 1){
    cat(paste0(rnames, ": "))
    modelprint(modeltypes)
  }else{
    for (i in 1:r){
      cat(paste0(rnames[i], ": "))
      modelprint(modeltypes[i], i)
      cat("\n")
    }
  }
}

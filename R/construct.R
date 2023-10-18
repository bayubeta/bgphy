
# create a function to define a model for Global/Mixed models
#' @export
setModel <- function(tree, modeltypes, startNodes = NULL){

  # check if tree is of class phylo
  stopifnot("Tree must be of class phylo." = class(tree) != "phylo")
  # check if elements of modeltypes is in {"OU", "BM"}
  stopifnot("Model name not found. Choose between BM or OU." = all(modeltypes %in% c("BM", "OU")))
  # check if startnodes are not null if the there are regimes
  stopifnot("Starting nodes for regimes are not specified." = (length(modeltypes) >= 1) & is.null(startNodes))
  # check startNodes values (numeric & in the tree)
  stopifnot("startNodes must be numeric values." = !is.numeric(startNodes))
  stopifnot("startNodes are not in the tree." = any(is.na(match(startNodes, tree$edge))))


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
                                           part.regime = setNames(rnames, startNodes),
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


  # set spacing size
  nchars <- sapply(rnames, nchar) # number of characters
  max_char <- max(c(7, max(nchars)))
  spaces <- max_char - nchars # number of additional spaces

  nspaces <- 4
  col_space <- strrep(" ", max_char - 7 + nspaces)
  def_space <- strrep(" ", nspaces - 1)


  # retrieve model types
  modeltypes <- attr(model, "modeltypes")


  cat(paste0("Regime", col_space,"Model\n"))
  cat(paste0("------", col_space,"-----\n"))
  if (r == 1){
    cat(paste0(rnames, strrep(" ", spaces), def_space))
    modelprint(modeltypes)
    cat("\n")
  }else{
    for (i in 1:r){
      cat(paste0(rnames[i], strrep(" ", spaces[i]), def_space))
      modelprint(modeltypes[i], i)
      cat("\n")
    }
  }

  cat("\n")
  cat("Priors: \n")
  print(model$priors)
}

#' Set a Gaussian phylogenetic model
#'
#' Defines a Gaussian phylogenetic model.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @param modeltypes A named vector of model types. The elements define the model types and names define the regime names.
#'    Currently only `"BM"` and `"OU"` are available as model types.
#' @param startNodes A named vector of nodes. Each element determines the start of a regime on the tree.
#'
#' @returns An object of class `bgphy_model`, which contains:
#' * `model`: The PCM model as defined by [PCMBase].
#' * `tree`: The phylogenetic tree. If there are multiple regimes, the tree is converted into `PCMTree` class.
#' * `priors`: A list of default priors. Each of them is of class `priorpdf`.
#'
#' @examples
#' \dontrun{
#' modeltypes <- setNames(c("OU", "OU"), c("ancestral", "new"))
#' startNodes <- setNames(c(101, 135), c("ancestral", "new"))
#' OU <- setModel(tree = lizardTree, modeltypes = ("OU"))
#' }
#'
#' @export
setModel <- function(tree, regime_names, modeltypes, startNodes = NULL){

  r <- length(regime_names) # number of regimes

  # # check if tree is of class phylo
  # stopifnot("Tree must be of class phylo." = class(tree) == "phylo")
  # # check if elements of modeltypes is in {"OU", "BM"}
  # stopifnot("Model name not found. Choose between BM or OU." = all(modeltypes %in% c("BM", "OU")))
  # # check if startnodes are not null if the there are regimes
  # if (length(modeltypes) > 1){
  #   stopifnot("Starting nodes for regimes are not specified." = !is.null(startNodes))
  #   # check startNodes values (numeric & in the tree)
  #   stopifnot("startNodes must be numeric values." = is.numeric(startNodes))
  #   stopifnot("startNodes are not in the tree." = !any(is.na(match(startNodes, tree$edge))))
  # }

  # create a PCM object given the tree and specified model and regimes
  r <- length(modeltypes) # number of regimes
  if (r == 1){
    # Global model
    PCMmodel <- PCMBase::PCM(model = modeltypes, regimes = regime_names)
  }else{
    # Mixed model
    # retrieve the names of regimes
    rnames <- names(modeltypes)
    # code names in PCMBase
    modelStrings <- ifelse(modeltypes == "BM",
                           "BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                           "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x")
    # Create the PCM object
    PCMmodel <- PCMBase::MixedGaussian(k = 1, modelTypes = modelStrings, mapping = stats::setNames(1:r, regime_names))

    # update the tree with the regimes information
    tree <- PCMBase::PCMTreeSetPartRegimes(PCMBase::PCMTree(tree),
                                           part.regime = unlist(sapply(1:r, function(i){stats::setNames(rep(names(startNodes)[i], length(startNodes[[i]])), startNodes[[i]])})),
                                           setPartition = TRUE, inplace = FALSE)
  }

  model <- structure(list(model = PCMmodel, tree = tree), class = "bgphy_model")
  attr(model, "modeltypes") <- modeltypes

  # set default priors for the parameters
  model$priors <- defaultPriors(model)

  return(model)
}


# print object of class bgphy (the constructed model)
#' Print model information
#'
#' This function prints the information of a model object of class `bgphy_model` on the console.
#'
#' @param post An object of class \code{bgphy_model}.
#'
#' @examples
#' \dontrun{
#' modeltypes <- setNames(c("OU", "OU"), c("ancestral", "new"))
#' startNodes <- setNames(c(101, 135), c("ancestral", "new"))
#' OU <- setModel(tree = lizardTree, modeltypes = ("OU"))
#' print(OU)
#' }
#'
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



#' @export
plot.PCMTree <- function(tree, cols = NULL, ...){

  if (is.null(cols)){
    cols <- c("red", "blue", "green1", "skyblue", "purple1")
  }

  # unique regime names
  ur <- unique(tree$part.regime)

  # selected colors
  col <- c()

  for (i in 1:length(tree$part.regime)){
    col[i] <- cols[which(ur == tree$part.regime[i])]
  }


  edge.col <- tree$edge.part
  for (i in 1:length(tree$part.regime)){
    edge.col[edge.col == names(tree$part.regime)[i]] <- col[i]
  }

  plot.phylo(tree, edge.col = edge.col, ...)
  # add legend
}



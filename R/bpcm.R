# new_bpcm <- function(X = double(), tree, model, prior = NULL){
#
# }


bpcm <- function(X, tree, model, prior = NULL){
  structure(list(X = X,
                 tree = tree,
                 model = model,
                 prior = prior),
            class = "bpcm")
}

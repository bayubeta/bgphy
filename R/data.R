#' Greater Antillean \emph{Anolis} tree
#'
#' A phylogeny of class \code{phylo} consisting of 100 different species of genus
#'   \emph{Anolis} inhabiting the Greater Antillean Islands. `lizardTreeMixed` is partitioned
#'   into two evolutionary regimes for simulation purposes.
#'
#' @name lizardTree
#'
#' @format An object of class \code{phylo}.
#'
#' @source Mahler, D. Luke, et al. "Exceptional convergence on the macroevolutionary
#'   landscape in island lizard radiations." \emph{Science} 341.6143 (2013): 292-295. <https://doi.org/10.1126/science.1232392>
"lizardTree"

#' @rdname lizardTree
#' @format `lizardTreeMixed` is of class \code{PCMTree} in addition to \code{phylo}.
"lizardTreeMixed"


#' Simulated datasets
#'
#' Two datasets simulated on the \emph{Anolis} phylogeny under different evolutionary models.
#'
#' @name datasets
#'
#' @format `XOU` is a matrix of size 100 x 100. Each row of was simulated on `lizardTree` under the `modelOU`.
#'
"XOU"

#' @rdname datasets
#' @format `XMixedOU` is a matrix of size 100 x 100. Each row of was simulated on `lizardTree` under the `modelMixedOU`.
"XMixedOU"



#' Single and Mixed Ornstein-Uhlenbeck model
#'
#' Ornstein-Uhlenbeck (OU) models of evolution. `modelOU` is a single (global) model
#'   while `modelMixedOU` is consisted of two OU models.
#'
#' @name models
#'
#' @format An object of class \code{PCM}.
#'
#' @source Mitov, Venelin, et al. "Fast likelihood calculation for multivariate Gaussian phylogenetic models with shifts."
#' \emph{Theoretical Population Biology} 131 (2020): 66-78. <https://doi.org/10.1016/j.tpb.2019.11.005>
"modelOU"

#' @rdname models
#' @format An object of class \code{PCM}.
"modelMixedOU"

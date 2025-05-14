# bgphy
 R package for Bayesian Inference of Mixed Gaussian Phylogenetic Models. 

 [![R-CMD-check](https://github.com/bayubeta/bgphy/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bayubeta/bgphy/actions/workflows/R-CMD-check.yaml)

## Installation
```
devtools::install_github("bayubeta/bgphy", build_vignettes = TRUE)
```

## Quickstart
(see `vignettes/bgphy_tutorial.pdf` for more information)

Loading the package:
```
library(bgphy)
```

Constructing a model with one regime:
```
# Brownian Motion (BM)
BM <- setModel(tree = lizardTree, regime_names = "R1", modeltypes = "BM")

# Ornstein-Uhlenbeck (OU)
OU <- setModel(tree = lizardTree, regime_names = "R1", modeltypes = "OU")
```

Print the defined model:
```
print(BM)
print(OU)
```

Constructing a model with multiple regimes:
```
BMOU <- setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                 modeltypes = c("BM", "OU"), 
                 startNodes = list(Ancestral = "101", New = "135"))
```

Visualizing regimes on the tree:
```
plot_regimes(BMOU$tree, no.margin = TRUE, x.lim = c(0,1.6), cex = 0.7)
```

Specifying priors:
```
BMOU$priors$alpha_2 <- prior_halfnormal(sigma = 2)
BMOU$priors$theta_2 <- prior_normal(mean = 2, sd = 5)
BMOU$priors$sigma_2 <- prior_halft(nu = 1, sigma = 2)
```

Running posterior inference:
```
# define the data as a row matrix with species names
X <- matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU)))

# posterior inference
post_BMOU <- bgphy(model = BMOU, X = X)

# print the results
print(post_BMOU)
```

Running posterior predictive check:
```
post_pred_check(post_BMOU)
```


## Reference
Brahmantio, B., Bartoszek, K., & Yapar, E. (2024). Bayesian inference of mixed Gaussian phylogenetic models. arXiv preprint arXiv:2410.11548.

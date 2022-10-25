library(ggtree)
# set.seed(420)
# tree = TreeSim::sim.bd.taxa.age(n = 50, age = 1, numbsim = 1, lambda=1, mu=0)[[1]]
# plot(tree)
load("data/tree.rda")


# set regimes
tree2 <- PCMTree(tree)
PCMTreeSetLabels(tree2)
PCMTreeSetPartRegimes(tree2, part.regime = c('65'= 2), setPartition = TRUE)

# plot tree
PCMTreePlot(tree2, size = 1, ladderize = TRUE) +
  geom_tiplab(size = 2) +
  geom_nodelab(size = 2, color = "black")

# save Tree
save(tree2, file = "data/treeMixed.rda")


# define the model
modelMixed = MixedGaussian(k = 1,
                           modelTypes = c("OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x",
                                          "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"),
                           mapping = c('1'=1, '2'=2))


# set parameters
# for regime 1
modelMixed$`1`$H[,,1] = 2
modelMixed$`1`$Theta[] = 2
modelMixed$`1`$Sigma_x[,,1] = 1

# for regime 2
modelMixed$`2`$H[,,1] = 5
modelMixed$`2`$Theta[] = 0.5
modelMixed$`2`$Sigma_x[,,1] = 0.5

print(PCMTable(modelMixed))

# save(modelMixed, file = "modelMixed.rda")

# simulate data
Xmixed = PCMSim(tree = tree2, model = modelMixed, X0 = modelMixed$X0[1])
plot(Xmixed[,tree2$tip.label])

Xtips2 = matrix(Xmixed[,tree2$tip.label], nrow = 1)
colnames(Xtips2) = tree2$tip.label
save(Xtips2, file = "data/Xmixed.rda")


#
# modelObject <- PCM("BM", k = 2L, regimes = c("a", "b", "c"))
# vec <- seq_len(PCMParamCount(modelObject))
# PCMParamLoadOrStore(modelObject, vec, offset = 0, load=TRUE)
#
#
#
# v = numeric(PCMParamCount(modelMixed))
# PCMParamLoadOrStore(modelMixed, v, load = FALSE, offset = 0)
# v

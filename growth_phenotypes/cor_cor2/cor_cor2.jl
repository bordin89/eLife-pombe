#=
Try transforming wideform growth phenotypes using cor and cor2
=#

using PombeAgeingGenes
using Distances, Statistics, Plots, DataFrames

df = load(GrowthPhenotypesWideform)

M = Matrix(df[:,Not(:id)])
c1 = cor(M')
c2 = cor(c1)

clusteredheatmap(M, linkage=:complete, c=:RdBu, clims=(-3,3))

clusteredheatmap(c1, metric=CorrDist(), linkage=:average, cor=true)
clusteredheatmap(c1, metric=CorrDist(), linkage=:complete, cor=true)

clusteredheatmap(c2, metric=CorrDist(), linkage=:average, cor=true)
clusteredheatmap(c2, metric=CorrDist(), linkage=:complete, cor=true)

#=
PCA
=#

using MultivariateStats

# Genes
m = fit(PCA, M')
transform(m, M')

m = fit(PCA, c1)
transform(m, c1)

m = fit(PCA, c2)
transform(m, c2)

# Conditions
fit(PCA, M)

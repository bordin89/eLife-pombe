#=
For Seaborn plots use `plot_clustermap.py`.
=#

cd(@__DIR__)

using PombeAgeingGenes, Statistics
using DataFrames, DataFramesMeta

# using Seaborn
# Seaborn.set(font_scale=0.3)
# pygui(true)

df = load(GrowthPhenotypesWideform)
A = Matrix(df[2:end])

# Conditions

C = cor(A)
clusteredheatmap(C, cor=true)

# fig = Seaborn.clustermap(C,
#     cmap="RdBu_r",
#     vmax=1.,
#     vmin=-1.,
#     xticklabels=string.(names(wf)[2:end]),
#     yticklabels=string.(names(wf)[2:end]))
#
# Seaborn.savefig("clustermap_conditions.pdf")

# ids

C = cor(A')
clusteredheatmap(C, cor=true)

# fig = clustermap(C,
#     cmap="RdBu_r",
#     vmax=1.,
#     vmin=-1.,
#     xticklabels=string.(df[:id]),
#     yticklabels=string.(df[:id]))
#
# Seaborn.savefig("clustermap_ids.png")

# ids x conditions

clusteredheatmap(A)

# fig = Seaborn.clustermap(A,
#     cmap="RdBu_r",
#     vmax=2.,
#     vmin=0.,
#     # cmap="viridis",
#     xticklabels=string.(names(df)[2:end]),
#     yticklabels=string.(df[:id]))
#
# Seaborn.savefig("clustermap_ids_conditions.png")

# PCA of genes, then clustermap of PCs

# using MultivariateStats
#
# X, _ = load(ML, Matrix)
#
# M = fit(PCA, X)
# Y = transform(M, X)
# clusteredheatmap(cor, Y)
#
# # PCA of conditions, then clustermap of PCs
#
# M = fit(PCA, X')
# Y = transform(M, X')
# clusteredheatmap(cor, Y)

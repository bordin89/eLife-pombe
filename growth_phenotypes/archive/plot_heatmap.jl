#=
Plot...
=#

cd(@__DIR__)

using PombeAgeingGenes, Plots, DataFrames, Clustering, Statistics, Distances

plotlyjs()

df = load(GrowthPhenotypesNoOutliersWideform)

A = convert(Matrix, df[2:end])
A ./= sum(A, dims=2)

# hc = hclust(pairwise(Euclidean(), C), linkage=:average)
#
# clustermap(A, hc, string.(names(df)[2:end]))
#     # cor=true,
#     # ticks=string.(names(df)[2:end]),
#     # tickfontsize=4)

using Seaborn, PyPlot

Seaborn.set(font_scale=0.5)

# conditions

C = cor(A)

fig = Seaborn.clustermap(C,
    cmap="RdBu_r",
    vmax=1.,
    vmin=-1.,
    xticklabels=string.(names(df)[2:end]),
    yticklabels=string.(names(df)[2:end]),
    figsize=(8, 8))

fig[:fig][:subplots_adjust](right=.75, top=.95, bottom=.25, left=.05)

Seaborn.savefig("clustermap_conditions.pdf")

# ids

C = cor(A')

fig = Seaborn.clustermap(C,
    cmap="RdBu_r",
    vmax=1.,
    vmin=-1.,
    xticklabels=string.(df[:id]),
    yticklabels=string.(df[:id]))

Seaborn.savefig("clustermap_ids.png")

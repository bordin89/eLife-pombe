include("src.jl")

using UMAP

df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
A = permutedims(A)
rescale!(A)

embedding = umap(A, 10)

d = 1
scatter(embedding[d,:], embedding[d+1,:], legend=false)

clusteredheatmap(embedding, clims=(-16,16), c=:RdBu_r)
heatmap(embedding, clims=(-16,16), c=:RdBu_r)

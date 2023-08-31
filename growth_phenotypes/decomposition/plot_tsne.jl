include("src.jl")

using TSne

# TSNE of genes
df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
A = permutedims(A)
rescale!(A)

Y = tsne(A, 2, 50, 1000, 20.0)

plottsne(Y)

savefig("tsne_ids.pdf")

# TSNE of conditions
df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
rescale!(A)

Y = tsne(A', 2, 50, 1000, 20.0)

plottsne(Y)

savefig("tsne_conditions.pdf")

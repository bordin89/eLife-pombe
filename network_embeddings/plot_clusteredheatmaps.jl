#=
Plot clusteredheatmaps of deepNF network embeddings and correlations.

NB for the network fusion thesis chapter, I used the Python script contained in this
directory to produce the figures
=#

cd(@__DIR__)

using PombeAgeingGenes, Statistics, Plots, DataFrames

module sns
    using Seaborn
end

# df, _ = load(ML, growthphenotypes=false, networkembeddings=true)

df = load(NetworkEmbeddings)

A = Matrix(df[:,Not(:id)])
vec(sum(A, dims=1)) .!= 0


A = A[:, vec(sum(A, dims=1)) .!= 0]

idxs = sortperm(vec(sum(A, dims=1)))

clusteredheatmap(
    A,#[:,idxs],
    xlabel="Embedding dimensions",
    ylabel="IDs",
    c=:viridis,
    xticks = false,
    yticks = false,
    # linkage=:complete,
    # cluster_x=false,
)

savefig(joinpath(@__DIR__, "network_embeddings.pdf"))

clusteredheatmap(cor(A'),
    xlabel="IDs",
    ylabel="IDs",
    cor=true,
    )

savefig(joinpath(@__DIR__, "network_embeddings_cor_ids.pdf"))

A_ = cor(A)
A_[isnan.(A_)] .= one(eltype(A_))

clusteredheatmap(A_,
    xlabel="Embedding dimensions",
    ylabel="Embedding dimensions",
    cor=true,
    )

savefig("network_embeddings_cor_embedding_dimensions.pdf")

# Seaborn

fig = sns.clustermap(
    A,
    xticklabels=false,
    yticklabels=false,
    cmap="viridis_r",
)
fig.ax_heatmap.set_xlabel("Embedding dimensions")
fig.ax_heatmap.set_ylabel("IDs")
sns.savefig("seaborn_network_embeddings.png")

fig = sns.clustermap(cor(A'),
    xticklabels=false,
    yticklabels=false,
    cmap="RdBu_r",
    vmax=1.,
    vmin=-1.,
    )
fig.ax_heatmap.set_xlabel("IDs")
fig.ax_heatmap.set_ylabel("IDs")
sns.savefig("seaborn_network_embeddings_cor_ids.png")

A_ = cor(A)
A_[isnan.(A_)] .= one(eltype(A_))

fig = sns.clustermap(A_,
    xticklabels=false,
    yticklabels=false,
    cmap="RdBu_r",
    vmax=1.,
    vmin=-1.,
    )
fig.ax_heatmap.set_xlabel("Embedding dimensions")
fig.ax_heatmap.set_ylabel("Embedding dimensions")
sns.savefig("seaborn_network_embeddings_cor_embedding_dimensions.png")

#=
PCA

https://stats.stackexchange.com/questions/82249/how-to-normalize-poisson-distributed-data-before-pca

https://stackoverflow.com/questions/10119913/pca-first-or-normalization-first
=#

include("src.jl")

# PCA of genes
df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
A = permutedims(A)
rescale!(A)

M = fit(PCA, A)
PCs = transform(M, A)

for i = 1:4
    plotpca(PCs, i)
    savefig("pca_ids_PC$(i).pdf")
end

# Highlight genes with annotations

# GO

for i = 1:length(ageinggoslimterms)
    goterm = ageinggoslimterms[i]
    ids = idswithgoterm(goterm)
    inds = [x ∈ ids for x = df[:id]]
    title = string(goterm)
    plotpca(PCs, 1, inds, title=title, legend=false)
    savefig("pca_genes_with_goterm/$title.pdf")
end

using PombeAgeingGenes.GeneOntology

# FYPO

for (fypoterm, ids) = GeneOntology.AGEING_FYPO_TERMS
    inds = [x ∈ ids for x = df[:id]]
    title = string(fypoterm)
    plotpca(PCs, 1, inds, title=title, legend=false)
    savefig("pca_genes_with_goterm/$title.pdf")
end


# PCA of conditions
df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
rescale!(A)

M = fit(PCA, A)
PCs = transform(M, A)

plotpca(PCs, 4)

savefig("pca_conditions.pdf")

#=
NB old code for a previous version of the data

# The points that are very different in the PCA are the six conditions that are very
# uncorrelated with the other conditions in the clustermap of the condition correlation
# matrix
selection = select(PCs, xlo=50, indices=true)
names(df)[2:end][selection]

#####################

#=
GO enrichment analysis using AnGeLi

http://bahlerweb.cs.ucl.ac.uk/AnGeLi
=#

using Rotations, DelimitedFiles

df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
A = permutedims(A)
rescale!(A)

M = fit(PCA, A)
PCs = transform(M, A)

# Rotate points around the Z axis, so that the clusters can be selected with a box in the
# XY plane, where the vertical edges are parallel to the Y axis and the horizontal edges
# are parallel to the X axis.
xyz = vcat(PCs[1:2,:], ones(size(PCs,2))')
r = RotZ(-0.04) # by trial and error
rxyz = r * xyz
plotpca(rxyz)

# Clusters
clusters = [
    (ylo=1.8, yhi=2.2),
    (ylo=1.3, yhi=1.7),
    (ylo=.8, yhi=1.2),
    (ylo=-2.7, yhi=-2.3),
    (ylo=-5, yhi=-2.7),
    ]

for i = 1:length(clusters)
    ylo, yhi = clusters[i]
    plotpca(select(rxyz, ylo=ylo, yhi=yhi), title="Cluster $i")
    savefig("pca_clusters/cluster$i.pdf")
    idsincluster(select(rxyz, ylo=ylo, yhi=yhi, indices=true), "pca_clusters/cluster$i.txt")
end

=#

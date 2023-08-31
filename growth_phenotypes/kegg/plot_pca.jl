using PombeAgeingGenes, Plots
plotlyjs()

include(joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "decomposition", "src.jl"))

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "kegg")

df = load(GrowthPhenotypesWideform)

A = convert(Matrix, df[2:end])
A = permutedims(A)
rescale!(A)

M = fit(PCA, A)
PCs = transform(M, A)

for pathway = KEGGPATHWAYS
    d = load(KEGGPathway(pathway))
    ids = keys(d)
    inds = [x ∈ ids for x = df[:id]]
    plotpca(PCs, 1, inds, title=pathway, legend=false)
    savefig(joinpath(dir, "pca", pathway*".pdf"))
    # break
end

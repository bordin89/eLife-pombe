using PombeAgeingGenes, Plots
plotlyjs()

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "kegg")

df = load(GrowthPhenotypesWideform)

# clusteredheatmaps

for pathway = KEGGPATHWAYS
    d = load(KEGGPathway(pathway))
    genes = keys(d)
    sub = @in(df, :id, genes)
    size(sub,1) < 2 && continue
    clusteredheatmap(Matrix(sub[2:end]),
        c=:RdBu_r, clims=(0.,2.),
        xlabel="Conditions", ylabel="IDs"
        )
    savefig(joinpath(dir, "size_clusteredheatmaps", pathway*".pdf"))
end

#=
Plot clustermaps of growth phenotype data and strain correlations:
    - ageing_genes_growth_*
    - ageing_genes_correlation_*
=#

using PombeAgeingGenes, Statistics
using PombeAgeingGenes.GeneOntology: AGEING_FYPO_TERMS

module sns
    using Seaborn
end

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "ageing")

df = load(GrowthPhenotypesWideform)

iter = [("FYPO0004344", AGEING_FYPO_TERMS[:FYPO0004344]),
        ("FYPO0001309", AGEING_FYPO_TERMS[:FYPO0001309]),
        ("all", [AGEING_FYPO_TERMS[:FYPO0004344]; AGEING_FYPO_TERMS[:FYPO0001309]])]

function color_mapper(id)
    cp = sns.Seaborn.seaborn.color_palette("colorblind")

    if id in AGEING_FYPO_TERMS[:FYPO0004344]
        if id in AGEING_FYPO_TERMS[:FYPO0001309]
            return cp[1]
        else
            return cp[2]
        end
    elseif id in AGEING_FYPO_TERMS[:FYPO0001309]
        return cp[3]
    else
        return "w"
    end
end

function plot_data(fname, df, ids)
    sub = @in(df, :id, ids)
    A = Matrix(sub[2:end])
    height = size(A,1)/5

    sns.clustermap(A,
        cmap="RdBu_r",
        yticklabels=string.(sub[:id]),
        xticklabels=string.(names(sub)[2:end]),
        # xticklabels=false,
        # yticklabels=false,
        vmax=2.,
        vmin=0.,
        row_colors=color_mapper.(sub[:id]),
        figsize=(12,height))
    sns.savefig(fname, bbox_inches="tight")
    sns.close()
end

for (fname, ids) = iter
    plot_data(joinpath(dir, "ageing_genes_growth_"*fname*"_clustermap.pdf"), df, ids)
end

# Correlation

function plot_correlation(fname, df, ids)
    sub = @in(df, :id, ids)
    A = Matrix(sub[2:end])
    c = cor(A')
    height = size(A,1)/5

    sns.clustermap(c,
        cmap="RdBu_r",
        # yticklabels=string.(sub[:id]),
        # xticklabels=string.(sub[:id]),
        xticklabels=false,
        yticklabels=false,
        vmax=1.,
        vmin=-1.,
        row_colors=color_mapper.(sub[:id]),
        figsize=(height,height))

    sns.savefig(fname, bbox_inches="tight")
    sns.close()
end

for (fname, ids) = iter
    plot_correlation(joinpath(dir, "ageing_genes_correlation_"*fname*"_clustermap.pdf"), df, ids)
end

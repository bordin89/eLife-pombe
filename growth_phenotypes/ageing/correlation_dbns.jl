using PombeAgeingGenes, Statistics, LinearAlgebra, Combinatorics, HypothesisTests, MultipleTesting, Printf, Plots, StatsPlots
using PombeAgeingGenes.GeneOntology: AGEING_FYPO_TERMS

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "ageing")

iter = [("FYPO0004344", AGEING_FYPO_TERMS[:FYPO0004344]),
        ("FYPO0001309", AGEING_FYPO_TERMS[:FYPO0001309]),
        ("all_ageing", unique([AGEING_FYPO_TERMS[:FYPO0004344]; AGEING_FYPO_TERMS[:FYPO0001309]]))]

df = load(GrowthPhenotypesWideform)
ids = df[:id]
A = Matrix(df[2:end])
C = cor(A')
mapping = Dict(ids[i] => i for i = 1:length(ids))

genes = AGEING_FYPO_TERMS[:FYPO0004344]
genes = AGEING_FYPO_TERMS[:FYPO0001309]
genes = unique([AGEING_FYPO_TERMS[:FYPO0004344]; AGEING_FYPO_TERMS[:FYPO0001309]])

# Correlations

function get_ageing_gene_idxs(genes, mapping)
    idxs = Vector{Int}[]
    for (g1,g2) = combinations(genes, 2)
        g1 in ids && g2 in ids || continue
        m1, m2 = mapping[g1], mapping[g2]
        # Ensure indices are in the upper triangle
        m1 > m2 ? push!(idxs, [m2,m1]) : push!(idxs, [m1,m2])
    end
    idxs
end

function plot_correlation_violin(genes, name)
    ageing_idxs = get_ageing_gene_idxs(genes, mapping)
    ageing_cor = map(idx->C[idx...], ageing_idxs)

    mask = triu(trues(size(C)), 1)
    foreach(idx->mask[idx...] = false, ageing_idxs)
    other_cor = C[mask]

    @assert length(other_cor) == ((size(C,1) * (size(C,1)-1)) ÷ 2) - length(ageing_idxs)

    p = HypothesisTests.pvalue(MannWhitneyUTest(ageing_cor, other_cor))

    violin([ageing_cor, other_cor],
        ylabel="Pearson correlation coefficient",
        c=:grey75,
        xticks=(1:2, ["Ageing genes", "Other genes"]),
        legend=false,
        title="Mann-Whitney p-value = "*@sprintf("%.3g", p)
        )

    Plots.savefig(joinpath(dir, "ageing_genes_correlation_$(name)_violin.pdf"))

    # histogram([ageing_cor, other_cor],
    #     xlabel="Pearson correlation coefficient",
    #     ylabel="Frequency",
    #     label=["Ageing genes", "Other genes"],
    #     norm=true,
    #     title="Mann-Whitney p-value = "*@sprintf("%.3g", p)
    #     )
    #
    # Plots.savefig(joinpath(dir, "ageing_genes_correlation_$(name)_histogram.pdf"))
end

for (name, genes) = iter
    plot_correlation_violin(genes, name)
end

# Cumulative distributions

correlation(C, mapping, g1, g2) = C[mapping[g1],mapping[g2]]

all_cor = C[triu(trues(size(C)), 1)]
sort!(all_cor)
all_cor_plot = all_cor[1:1000:length(all_cor)]

plotlyjs()

Plots.plot(all_cor_plot, range(0,1;length=length(all_cor_plot)), xticks=collect(-1:.2:1), label="All genes", xlims=(-1,1))

for (name, genes) = [("Nitrogen starvation ageing", AGEING_FYPO_TERMS[:FYPO0004344]),
                     ("Stationary phase ageing", AGEING_FYPO_TERMS[:FYPO0001309]),
                     ("All ageing", unique([AGEING_FYPO_TERMS[:FYPO0004344];
                                            AGEING_FYPO_TERMS[:FYPO0001309]]))]
    genes = intersect(genes, ids)
    cors = map(gs->correlation(C, mapping, gs...), combinations(genes, 2))
    Plots.plot!(sort(cors), range(0,1;length=length(cors)), label=name*" genes")
end

Plots.plot!(
    xlabel="Pearson correlaton coefficient",
    ylabel="Cumulative relative frequency",
    legend=true,
    left_margin=5Plots.PlotMeasures.mm,
    )

Plots.savefig(joinpath(dir, "ageing_genes_correlation_cumulative_dbn.pdf"))


# Comparison with network embeddings
plotlyjs()

df, _ = load(ML, growthphenotypes=false, networkembeddings=true)
ids = df[:id]
A = Matrix(df[2:end])
C = cor(A')
mapping = Dict(ids[i] => i for i = 1:length(ids))

all_cor = C[triu(trues(size(C)), 1)]
sort!(all_cor)
all_cor_plot = all_cor[1:1000:length(all_cor)]

Plots.plot(all_cor_plot, range(0,1;length=length(all_cor_plot)), xticks=collect(-1:.2:1), label="All genes", xlims=(-1,1))

for (name, genes) = [("Nitrogen starvation ageing", AGEING_FYPO_TERMS[:FYPO0004344]),
                     ("Stationary phase ageing", AGEING_FYPO_TERMS[:FYPO0001309]),
                     ("All ageing", unique([AGEING_FYPO_TERMS[:FYPO0004344];
                                            AGEING_FYPO_TERMS[:FYPO0001309]]))]
    genes = intersect(genes, ids)
    cors = map(gs->correlation(C, mapping, gs...), combinations(genes, 2))
    Plots.plot!(sort(cors), range(0,1;length=length(cors)), label=name*" genes")
end

Plots.plot!(
    xlabel="Pairwise Pearson correlaton coefficient",
    ylabel="Cumulative relative frequency",
    )

Plots.savefig(joinpath(dir,
                       "ageing_genes_networkembeddings_correlation_cumulative_dbn.pdf"))

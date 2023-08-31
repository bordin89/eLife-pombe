#=
Calculate p-values for network statistics:
    - degree
    - betweenness_centrality

Plot p-values against correlation of growth phenotypes across all condtions for each gene
compared to the wt.
=#
using PombeAgeingGenes, Plots, Statistics, DataFrames

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "networks")

df = load(GrowthPhenotypesWideform)

A = Matrix(df[2:end])
c = cor(A, dims=2)

# Correlation of each gene to wt
x = c[:,end]
histogram(x, xlabel="Correlation to wt", ylabel="Frequency")
savefig(joinpath(dir, "wt_cor_histogram.pdf"))

df = DataFrame(id = df[:id], cor = x)

sort!(df, :cor)

using LightGraphs, CodecZlib

################
# Maybe move this section into PombeAgeingGenes.jl

import PombeAgeingGenes: @file
import FileIO: load

@file STRINGNetwork joinpath(ENV["POMBEAGEINGGENES"], "data", "networks",
                             "4896.protein.links.v11.0.txt.gz")

function load(x::STRINGNetworkFile)
    g = Graph()
    mapping = Dict{String,Int}()
    i = 1

    open(GzipDecompressorStream, filepath(x), "r") do io
        for l = eachline(io)
            startswith(l, "protein1") && continue
            p1, p2, s = split(l)

            p1 = p1[6:end-2]
            p2 = p2[6:end-2]

            for p = (p1, p2)
                if !haskey(mapping, p)
                    mapping[p] = i
                    add_vertex!(g)
                    i += 1
                end
            end

            add_edge!(g, mapping[p1], mapping[p2])
        end
    end

    g, mapping
end

################

using DelimitedFiles, GLM

g, mapping = load(STRINGNetwork)

# bc = betweenness_centrality(g)
# writedlm(joinpath(ENV["POMBEAGEINGGENES"], "data", "networks",
#                   "4896.protein.links.v11.0.bc.txt"), bc)

d = degree(g)
histogram(d)

# Degree

function degree_pvalue(gene::AbstractString)
    if haskey(mapping, gene)
        g_idx = mapping[gene]
        PValue(d[g_idx], null).p
    else
        missing
    end
end

d = degree(g)
null = NullDistribution(d)
df[:p] = map(gene->degree_pvalue(gene), df[:id])

let x = dropmissing(df)
    ols = lm(@formula(p ~ cor), df)

    scatter(x[:cor], x[:p],
        xlabel="Correlation to wt", ylabel="Degree p-value",
        legend=false,
        title="r² = $(round(r2(ols), digits=3))"
        )
    savefig(joinpath(dir, "cor_against_degree_pvalue.pdf"))
end

# betweenness_centrality

function bc_pvalue(gene::AbstractString)
    if haskey(mapping, gene)
        g_idx = mapping[gene]
        PValue(bc[g_idx], null).p
    else
        missing
    end
end

bc = readdlm(joinpath(ENV["POMBEAGEINGGENES"], "data", "networks",
             "4896.protein.links.v11.0.bc.txt"))[:,1]
null = NullDistribution(bc)
df[:p] = map(gene->bc_pvalue(gene), df[:id])

let x = dropmissing(df)
    ols = lm(@formula(p ~ cor), df)

    scatter(x[:cor], x[:p],
        xlabel="Correlation to wt", ylabel="Betweenness centrality p-value",
        legend=false,
        title="r² = $(round(r2(ols), digits=3))"
        )
    savefig(joinpath(dir, "cor_against_bc_pvalue.pdf"))
end

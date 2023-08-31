#=
Plots in ./pccs/
=#

using PombeAgeingGenes, Statistics, Combinatorics, LinearAlgebra, JSON, HypothesisTests,
    Printf, Plots, StatsPlots

plotlyjs()

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "growth_phenotypes", "kegg")

# Make data set

df = load(GrowthPhenotypesWideform)

A = Matrix(df[2:end])

c = cor(A, dims=2)

null = NullDistribution(c[triu(trues(size(c)), 1)])

function gethits(c, null, ids; threshold=.01)
    hits = Dict{String,Vector{Tuple{String,String,Float64}}}()
    mapping = Dict(ids[i] => i for i = 1:length(ids))

    for pathway = KEGGPATHWAYS
        @show pathway
        hits[pathway] = Tuple{String,String,Float64}[]
        genes = collect(keys(load(KEGGPathway(pathway))))

        for (g1,g2) = combinations(genes, 2)
            g1 in ids && g2 in ids || continue
            p = PValue(c[mapping[g1], mapping[g2]], null)
            if p.p < threshold
                push!(hits[pathway], (g1,g2,p.p))
            end
        end

        length(hits[pathway]) == 0 && delete!(hits, pathway)
    end

    hits
end

threshold = 0.01
@time hits = gethits(c, null, df[:id], threshold=threshold)
write(joinpath(dir, "significant_hits_threshold_$threshold.json"), JSON.json(hits))

threshold = 1.
hits = gethits(c, null, df[:id], threshold=threshold)
write(joinpath(dir, "significant_hits_threshold_$threshold.json"), JSON.json(hits))

########################

# Plot distributions

df = load(GrowthPhenotypesWideform)
A = Matrix(df[2:end])
c = cor(A, dims=2)
null = c[triu(trues(size(c)), 1)]

threshold = 1.
d = JSON.parsefile(joinpath(dir, "significant_hits_threshold_$threshold.json"))

# violin plots

for pathway = KEGGPATHWAYS
    !(pathway in keys(d)) && continue
    cors = getindex.(d[pathway], 3)
    p = pvalue(MannWhitneyUTest(cors, null))

    violin([cors, null],
        ylabel="Pearson correlation coefficient",
        c=:grey75,
        xticks=(1:2, ["KEGG pathway", "Null distribution"]),
        legend=false,
        title="Mann-Whitney p-value = "*@sprintf("%.3g", p)
    )

    savefig(joinpath(dir, "pccs", pathway*".pdf"))
    # break
end

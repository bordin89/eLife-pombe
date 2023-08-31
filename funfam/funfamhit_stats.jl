#=
Calculate useful stats about FunFamHits.
=#

using PombeAgeingGenes, DataFrames, Statistics, Plots

plotlyjs()

hits = load(FunFamHits, threshold=1e-3)
df = DataFrame(hits)
unique!(df, [:id, :ff])
df.score = -log10.(df.score)


# Number of proteins with significant hits

length(unique(df.id))

# Number of significant hits

size(df, 1)

# Hits per FunFam

hits_per_ff = by(df, :ff, n_hits = :id => length)

sort!(hits_per_ff, :n_hits)

@show(hits_per_ff[end,:])

n_hits = hits_per_ff.n_hits

median(n_hits)
mean(n_hits)
maximum(n_hits)
quantile(n_hits, 0.95)

plot(
    [0; n_hits],
    (0:length(n_hits)),
    xlabel = "# hits per FunFam",
    ylabel = "Cumulative frequency",
    linetype = :steppost,
    xticks = let x = maximum(n_hits)
        ticks = [0:20:x..., x]
        tick_labels = string.(ticks)
        tick_labels[end-1] = ""
        (ticks, tick_labels)
    end,
    yticks = let y = length(n_hits)
        [0:5000:y..., y]
    end,
    legend = false,
    c = :black,
    lw = 2,
    size = (400,300),
    right_margin = 2Plots.PlotMeasures.mm,
)

# histogram(
#     n_hits,
#     xlabel="Hits per FunFam",
#     ylabel="Frequency",
#     yaxis=:log,
#     legend=false,
#     bins=0:5:maximum(n_hits),
#     c=:grey80,
#     size=(400,300),
# )

savefig(joinpath(@__DIR__, "hits_per_funfam_histogram.pdf"))


# Number of FunFams with only one significant hit

count(hits_per_ff.n_hits .== 1)


# FunFams per GO Slim

using PombeAgeingGenes.GeneOntology, OBOParse, Format

go = load(GO)

d = PombeAgeingGenes.map_goterms_to_funfams()

ffs_per_go_term = Int[]
for t = GO_SLIM_TERMS
    haskey(d, t) && push!(ffs_per_go_term, length(unique(d[t])))
end

sort!(ffs_per_go_term)

plot(
    [0; ffs_per_go_term],
    (0:length(ffs_per_go_term)),
    xlabel = "# FunFams associated with GO Slim term",
    ylabel = "Cumulative frequency",
    linetype = :steppost,
    xticks = let ticks = [0:2500:10000..., maximum(ffs_per_go_term)]
        (ticks, format.(ticks, commas=true))
    end,
    yticks = [0:10:50..., length(ffs_per_go_term)],
    xformatter = :plain,
    legend = false,
    c = :black,
    lw = 2,
    size = (400,300),
    right_margin = 2Plots.PlotMeasures.mm,
)

stephist(
    ffs_per_go_term,
    xlabel="# FunFams associated with GO Slim term",
    ylabel="Frequency",
    bins=0:500:maximum(ffs_per_go_term)+2000,
    legend = false,
    c = :black,
    size = (450, 300),
    lw = 2,
)

savefig(joinpath(@__DIR__, "funfams_per_go_slim_term.pdf"))

# Thesis

Int(median(ffs_per_go_term))
mean(ffs_per_go_term)
maximum(ffs_per_go_term)

@show GO_SLIM_TERMS[ffs_per_go_term .== maximum(ffs_per_go_term)] # GO:0016192

instance = length(descendants(go, go.terms["GO:0016192"]))
null = map(t->length(descendants(go, t)), collect(values(go.terms)))
PValue(instance, null)

#=
Multiple testing

Definitions:
- Type I error: rejection of a true H0, a false positive
- Power: probability that H0 is rejected when H1 is true

FDR methods:
- Controls expected proportion of rejected H0s that are true
- Less strict control of Type I errors
- More power
- e.g. Benjamini-Hochberg
    1. sort p-values P(1), P(2),...,P(J)
    2. reject first j hypotheses for which P(j) ≤ jq / J, where q is the FDR level and J is the number of hypotheses being tested (step-up procedure).

FWER methods:
- Controls probability of at least one Type I error
- More strict control of Type I errors
- Less power
- e.g. Bonferroni
- FWER = 1 − (1 − α)^J, where α is the significance threshold and J is the number of hypotheses being tested. With α = 0.05 and J = 100, FWER = 0.994.

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3859974/
=#

using HypothesisTests, MultipleTesting, Statistics, StatsPlots, PombeAgeingGenes, DataFrames, DataFramesMeta, CSVFiles

plotlyjs()

df = load(GrowthPhenotypesNoOutliers)


# Significance testing

null_dbns = Dict{String,Dict{String,Vector{Float64}}}()
null_dbns["YES_32"] = Dict{String,Vector{Float64}}()
null_dbns["EMM_32"] = Dict{String,Vector{Float64}}()

for g in groupby(@in(df, :condition, ["YES_32", "EMM_32"]), [:id, :condition])
    id = g.id[1]
    condition = g.condition[1]
    sizes = collect(g.size)
    null_dbns[condition][id] = sizes
end

df = by(df, [:id, :condition]) do g
    id = g.id[1]
    condition = g.condition[1]
    sizes = g.size
    null = null_dbns[startswith(condition, "YES") ? "YES_32" : "EMM_32"][id]
    p = try
        pvalue(UnequalVarianceTTest(g.size, null))
    catch ArgumentError
        1.
    end
    (mean_size = mean(sizes), p_value = p)
end

save(PombeAgeingGenes._fname * "_p_values.csv", df)


# Multiple testing correction

df[!,:p_corr] .= 0.

for g = groupby(df, :condition)
    g[:,:p_corr] = adjust(collect(g[:,:p_value]), BenjaminiHochberg())
end

plot()
stephist!(
    df.p_value,
    label="Uncorrected",
    bins=0:.01:1,
    c=:grey70,
)
stephist!(
    df.p_corr,
    label="Corrected",
    bins=0:.01:1,
    c=:black,
)
plot!(
    xlabel="P",
    ylabel="Frequency",
    size=(400,300),
)

savefig(joinpath(@__DIR__, "significance_testing_p_values_histogram.pdf"))

# Hits per strain and condition

hits_condition = by(
    df,
    :condition,
    hits = :p_corr => x->count(x .< .01),
)

sort!(hits_condition, :hits, rev=true)

fig1 = stephist(
    hits_condition.hits,
    xlabel="Hits per condition",
    ylabel="Frequency",
    legend=false,
    bins=0:10:maximum(hits_condition.hits),
    c=:black,
    size=(300,300),
)

savefig(joinpath(@__DIR__, "significance_testing_hits_per_condition_histogram.pdf"))


hits_id = by(
    df,
    :id,
    hits = :p_corr => x->count(x .< .01),
)

sort!(hits_id, :hits, rev=true)

fig2 = stephist(
    hits_id.hits,
    xlabel="Hits per strain",
    ylabel="Frequency",
    legend=false,
    bins=0:maximum(hits_id.hits),
    c=:black,
    size=(300,300),
)

savefig(joinpath(@__DIR__, "significance_testing_hits_per_id_histogram.pdf"))


# plot(
#     fig1,
#     fig2,
#     size=(700,300),
# )
#
# savefig(joinpath(@__DIR__, "significance_testing_hits_per_condition_id_histogram.pdf"))


# Thesis

df_sig = @where(df, :p_corr .< .01)

unique(df_sig.condition)
unique(df.condition)

unique(df_sig.id)
unique(df.id)

mean(hits_condition.hits)
median(hits_condition.hits)

mean(hits_id.hits)
median(hits_id.hits)

for row in eachrow(@where(hits_id, :hits .> 20))
    println(row.id, " & ", row.hits, " \\\\")
end



# const NULL_YES = @where(df, :condition .== "YES_32")
# const NULL_EMM = @where(df, :condition .== "EMM_32")

# df = by(@where(df, :condition .!= "YES_32", :condition .!= "EMM_32"), :id) do g
#     id = g.id[1]
#     condition = g.condition[1]
#
#     null = @where(
#         startswith(condition, "YES") ? NULL_YES : NULL_EMM,
#         :id .== id).size
#
#     by(g, :condition) do h
#         p = try
#             pvalue(UnequalVarianceTTest(h.size, null))
#         catch ArgumentError
#             1.
#         end
#         (meansize = mean(h.size), pvalue = p)
#     end
# end

# function correct!(df, method)
#     df[!,:corrected] = adjust(df[:pvalue], method)
#
#     df[!,:corrected_condition] .= 0.
#     for g = groupby(df, :condition)
#         g[:,:corrected_condition] = adjust(collect(g[:,:pvalue]), method)
#     end
#
#     df[!,:corrected_id] .= 0.
#     for g = groupby(df, :id)
#         g[:,:corrected_id] = adjust(collect(g[:,:pvalue]), method)
#     end
#
#     df = by(df, [:id, :condition]) do g
#         (pvalue = minimum(g.pvalue),
#          corrected = minimum(g.corrected),
#          corrected_condition = minimum(g.corrected_condition),
#          corrected_id = minimum(g.corrected_id),
#         )
#     end
#
#     return df
# end
#
# correct(df, method) = correct!(deepcopy(df), method)

# function plot_pvalues(df)
#     violin(["raw"], df[:pvalue])
#     violin!(["corrected globally"], df[:corrected])
#     violin!(["corrected per condition"], df[:corrected_condition])
#     violin!(["corrected per id"], df[:corrected_id])
#     plot!(legend=false, ylabel="p-value")
# end
# function plot_pvalues(df, method)
#     fig = plot_pvalues(correct(df, method))
#     method = string(method)[1:end-2]
#     plot!(title=string(method))
#     savefig("multiple_testing_correction_$(method).pdf")
#     fig
# end
#
# plot_pvalues(df, Bonferroni())
# plot_pvalues(df, BenjaminiHochberg())
#
# method = Bonferroni
# method = BenjaminiHochberg
#
# corrected = correct(df, method())

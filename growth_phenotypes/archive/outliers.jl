# include("common.jl")
using PombeAgeingGenes, DataFrames, DataFramesMeta

# Compare data set before and after removing outliers
df = load(GrowthPhenotypes)
df = @where(df, :condition .== "YES_32")
dfno = removeoutliers(df)
c = PlotThemes.wong_palette;

    # Kurtosis
k = by(df, [:condition, :id], df->kurtosis(df[:size]))
k = k[.!isnan.(k[:x1]),:]
kno = by(dfno, [:condition, :id], df->kurtosis(df[:size]))
kno = kno[.!isnan.(kno[:x1]),:]

opts = (α=.5, bins=minimum([k[:x1]; kno[:x1]]):.25:maximum([k[:x1]; kno[:x1]]))
fig = histogram(k[:x1]; label="With outliers", color=c[1], xlabel="Kurtosis (YES_32)", ylabel="Frequency", opts...)
histogram!(kno[:x1]; label="Outliers removed", color=c[2], opts...)
savefig("kurtosis_hist.pdf")

    # Skewness
s = by(df, [:condition, :id], df->skewness(df[:size]))
s = s[.!isnan.(s[:x1]),:]

sno = by(dfno, [:condition, :id], df->skewness(df[:size]))
sno = sno[.!isnan.(sno[:x1]),:]

opts = (α=.5, bins=minimum([s[:x1]; sno[:x1]]):.25:maximum([k[:x1]; kno[:x1]])) #norm=true)
fig = histogram(s[:x1]; label="With outliers", color=c[1], xlabel="Skewness (YES_32)", ylabel="Frequency", opts...)
histogram!(sno[:x1]; label="Outliers removed", color=c[2], opts...)
savefig("skewness_hist.pdf")

    # Entropy
df_entropy = ent(df)
dfno_entropy = ent(removeoutliers(df))

plotlyjs()
c = PlotThemes.wong_palette;
opts = (α=.5,)
violin(df_entropy[:x1], label="With outliers", c=c[1]; xticks=false, ylabel="Entropy", opts...)
violin!(dfno_entropy[:x1], label="Outliers removed", c=c[2]; opts...)
savefig("entropy.pdf")

pvalue(MannWhitneyUTest(df_entropy[:x1], dfno_entropy[:x1]))

# CDF of sizes per id for YES_32 with outliers removed
df = @where(loaddf(), :condition .== "YES_32")
removeoutliers!(df)

opts = (α=.1,)
let
    fig = plot()
    by(df, :id) do x
        v = sort(x[:size])
        plot!(v,cdf.(Normal(), zscore(v)); ylabel="Cumulative frequency", xlabel="Relative size (YES_32)", legend=false, color=c[2], opts...)
    end
    fig
    savefig("size_cdf.pdf")
end

df1 = by(df, :id, x->mean(x[:size]))
histogram(df1[:x1], xlabel="Mean relative size (YES_32)", ylabel="Frequency", legend=false, color=c[2])
savefig("size_hist.pdf")

mode(df1[:x1])
median(df1[:x1])
mean(df1[:x1])

# Violin plots before and after outlier removal for YES_32
plotviolin(index, data; kwargs...) = violin!([index], data, trim=false, linecolor=false; kwargs...)
plotlhs(index, data; kwargs...) = plotviolin(index, data, side=:left, color=c[1]; kwargs...)
plotrhs(index, data; kwargs...) = plotviolin(index, data, side=:right, color=c[2]; kwargs...)

df = loaddf()
scale = 3
c = PlotThemes.wong_palette;
ctr = 0
i = 0
outliers = String[]
fig = plot()
while ctr < 20
    global i += 1
    x = Float64.(@where(df, :condition .== "YES_32", :id .== ids[i])[:size])
    length(x) == 0 && continue
    F = Fence(x, scale)
    y = [x for x = x if !isoutlier(F, x)]
    if length(y) < length(x)
        global ctr += 1
        plotlhs(ctr, x, label="", xticks=(i, [ids[i]]))
        plotrhs(ctr, y, label="")
        push!(outliers, ids[i])
    end
end
xticks!(1:20, outliers, xrotation=45) #, bottom_margin=1mm)
ylabel!("Relative size")
xlabel!("Strains")
savefig(fig, "size_outliers_violin.pdf")

# # KS tests of normality
# # p < 0.05 == normally distributed
#
# using HypothesisTests, Distributions
#
# df = @where(loaddf(), :condition .== "YES_32")
# dfno = removeoutliers(df)
#
# ks = by(df, :id) do x
#     if size(x, 1) > 1
#         pvalue(ExactOneSampleKSTest(zscore(x[:size]), Normal()))
#     else
#         missing
#     end
# end
# dropmissing!(ks)
#
# ksno = by(dfno, :id) do x
#     if size(x, 1) > 1
#         pvalue(ExactOneSampleKSTest(zscore(x[:size]), Normal()))
#     else
#         missing
#     end
# end
# dropmissing!(ksno)
#
# c = PlotThemes.wong_palette;
# opts = (α=.5, bins=0:.025:1)
# histogram(ks[:x1]; label="With outliers", color=c[1], xlabel="KS p-value (YES_32)", ylabel="Frequency", legend=:topleft, opts...)
# histogram!(ksno[:x1]; label="Outliers removed", color=c[2], opts...)
# savefig("ks_hist.pdf")

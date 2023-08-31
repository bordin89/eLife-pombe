#=
Plot heatmaps that show (normalised) histograms for the frequency/probability of colonies in histogram bins.
=#

cd(@__DIR__)

using PombeAgeingGenes, DataFrames, Plots, Clustering, Distances, Statistics

@userplot HeatmapHistogram

recipetype(::Val{:heatmaphistogram}, args...; kwargs...) = HeatmapHistogram(args, kwargs)

@recipe function f(hh::HeatmapHistogram)#; kwargs...)
    hc, A, xticks, yticks = hh.args
    seriestype := :heatmap
    xticks --> (1:length(xticks), xticks)
    yticks --> (1:length(yticks), yticks[hc.order])
    c --> :viridis
    tickfontsize --> 5
    size --> (600,800)
    grid --> false
    A[hc.order,:]
end

plotlyjs()

df = load(GrowthPhenotypesWideform)

df = stack(df)

rename!(df, :variable=>:condition, :value=>:size)

# Histogramify
df[!,:bin] = round.(df[:,:size], digits=1)

g = by(df, [:condition, :bin], count = :bin => length)

wf = unstack(g, :condition, :bin, :count)

xticks = string.(names(wf)[2:end])
yticks = string.(wf[:condition])

# Heatmap of colony frequencies per bin (histogram)

A = convert(Matrix, wf[:,2:end])
A_0 = coalesce.(A, 0)
hcA = hclust(pairwise(Euclidean(), A_0'), linkage=:average)

heatmaphistogram(hcA, log10.(coalesce.(A, NaN)), xticks, yticks, colorbar_title = "Freq log10")
savefig(joinpath(@__DIR__, "heatmap_histogram_per_condition.pdf"))

# Heatmap of colony probabilities per bin (normalised histogram)

P = A ./ sum(coalesce.(A, 0), dims=2)
P_0 = coalesce.(P, 0.)
hcP = hclust(pairwise(Euclidean(), P_0'), linkage=:average)

heatmaphistogram(hcP, coalesce.(P, NaN), xticks, yticks, colorbar_title="Probability")
savefig(joinpath(@__DIR__, "heatmap_normalised_histogram_per_condition.pdf"))

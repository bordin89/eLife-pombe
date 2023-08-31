cd(@__DIR__)

using PombeAgeingGenes, DataFrames, DataFramesMeta

# include("growth_phenotypes_src.jl")

df = load(GrowthPhenotypes)

ids = collect(Set(String.(df[:id])))
conditions = collect(Set(String.(df[:condition])))

using StatsBase, Statistics, Combinatorics, DataStructures, HypothesisTests, Distributions
using Plots, StatsPlots

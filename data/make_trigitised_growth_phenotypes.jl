using PombeAgeingGenes, CSV

df = load(GrowthPhenotypesWideform)

dft = trigitise(df)

fname = "Jan2019_BBSRC_results_no_outliers_wideform_trigitised.csv"

CSV.write(joinpath(ENV["POMBEAGEINGGENES"], "data", fname), dft)


# Clustering

using Plots

A = Matrix(dft[2:end])

clusteredheatmap(A, c=:RdBu_r, clims=(-1.5,1.5))

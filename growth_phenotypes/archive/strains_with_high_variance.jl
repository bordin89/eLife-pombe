using PombeAgeingGenes, Statistics, Plots

df = load(GrowthPhenotypesWideform)

df[:var] = [var(df[i, 2:end]) for i = 1:size(df,1)]

@show sort(df, :var, rev=true)[ 1:10, [:id, :var]]

histogram(df[:var])

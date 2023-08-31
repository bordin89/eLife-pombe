using PombeAgeingGenes, StatsPlots

df = load(GrowthPhenotypesWideform)

@df df andrewsplot(:id, cols(2:size(df,2)))

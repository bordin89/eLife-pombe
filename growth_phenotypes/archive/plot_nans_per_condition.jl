using PombeAgeingGenes, DataFrames, Plots

df = load(GrowthPhenotypesNoOutliers)

# Mean sizes
# NB sizes for strain-condition pairs with ≤ `nrepeats` repeats are set to NaN
df = meansizes(df; nrepeats=2)

# Remove columns with lots of NaNs
df_nnan = by(df, :condition, nnan = :size => x->count(isnan.(x)))

histogram(df_nnan[:nnan],
    xlabel="NaNs per condition",
    ylabel="Frequency (# conditions)",
    bins=100,
    c=:grey50,
    legend=false,
    )
vline!([3000], c=:black, ls=:dot, lw=2)

savefig("sizes_nan_histogram.pdf")

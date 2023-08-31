using PombeAgeingGenes, DataFrames, Plots

df = load(GrowthPhenotypesNoOutliers)

# Mean sizes. NB sizes for strain-condition pairs with ≤ `nrepeats` repeats are set to NaN
df = meansizes(df; nrepeats=2)

# Remove columns with lots of NaNs
df_nnan = by(df, :condition, nnan = :size => x->count(isnan.(x)))
conditions_to_delete = @where(df_nnan, :nnan .> 3000)[:condition]
deleterows!(df, map(x->x ∈ conditions_to_delete, df[:condition]))

# Impute NaNs with mean size per condition
impute!(df)

# Wideform
wf = unstack(df, :id, :condition, :size)

# Coalesce missings
dropmissing(wf) # number of complete cases

nmissing = colwise(x->count(ismissing.(x)), wf)

histogram(nmissing,
    xlabel="missings per condition",
    ylabel="Frequency (# conditions)",
    bins=100,
    c=:grey50,
    legend=false,
    )

savefig("sizes_missing_histogram.pdf")

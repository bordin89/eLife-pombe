using PombeAgeingGenes, DataFrames, DataFramesMeta, Plots, CSV

df = load(GrowthPhenotypes)

df = unique(df[[:condition, :repeat_number, :phlox, :library_version]])

# conditions
unique(df[:condition])

# repeats per condition
repeats_per_condition = by(df, :condition, nrepeats = :repeat_number => length)

sort!(repeats_per_condition, :nrepeats)

CSV.write(joinpath(@__DIR__, "conditions_and_repeats.csv"), repeats_per_condition)

@where(repeats_per_condition, :nrepeats .== 2)

histogram(repeats_per_condition[:nrepeats], xlabel="Repeats per condition", ylabel="# conditions")
savefig(joinpath(@__DIR__, "repeats_per_condition.pdf"))

# Conditions that are not correlated with most of the other conditions
# See `clustermap_conditions.pdf` for reference

conditions = ["YES_Diamide_3mM",
    "YES_KCl_0.5M_SDS_0.04percent",
    "YES_Diamide_2mM",
    "YES_tunicamycin_1",
    "YES_NaCl_100mM_SDS0.04percent",
    "YES_tunicamycin_2",
    "EMM_DIP",
    "YES_LiCl_4mM_SDS_0.04percent",
    "YES_LiCl_7.5mM_SDS_0.04percent",
    "YES_EGTA_10mM",
    "YES_EGTA_10mM_5days"]

tmp = @in(repeats_per_condition, :condition, conditions)
CSV.write(joinpath(@__DIR__, "uncorrelated_conditions.csv"), tmp)

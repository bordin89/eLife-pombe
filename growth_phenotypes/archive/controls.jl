include("common.jl")

# Size of control strains w/ and w/o outliers
df = load(GrowthPhenotypes)
df = @where(df, :condition .== "YES_32")
dfno = removeoutliers(df)
df_controls = @where(df, map(x->x in controls, :id))
df_outliers = outliers(df_controls, keepoutliers=true)
dfno_controls = removeoutliers(df_controls)

plotlyjs()
c = PlotThemes.wong_palette;
opts = (trim=false, linecolor=false)
fig = @df df_outliers scatter(:id, :size, label="Outliers", markerstrokewidth=0, c=:grey50, α=0.25)

# TEMP while StatsPlots.jl issue #198 is open
@df df_controls violin(:id, :size, side=:left, label="With outliers", c=c[1]; opts...)
@df dfno_controls violin(:id, :size, side=:right, label="Outliers removed", c=c[2]; opts...)

@df df_controls violin!(:id, :size, side=:left, label="With outliers", c=c[1]; opts...)
@df dfno_controls violin!(:id, :size, side=:right, label="Outliers removed", c=c[2]; opts...)
plot!(ylabel="Relative size", xrotation=45, xlabel="Control strains", legend=true)
savefig("control_size_violin.pdf")

# CV
by(df_controls, :id, x->cv(x[:size]))
by(dfno_controls, :id, x->cv(x[:size]))

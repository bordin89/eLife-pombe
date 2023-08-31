using PombeAgeingGenes, Plots, Statistics, DataFrames, DataFramesMeta

plotlyjs()

df = load(GrowthPhenotypesWideform)

df = DataFrame(
    condition = names(df)[2:end],
    mean = map(mean, eachcol(df[:,2:end])),
    var = map(var, eachcol(df[:,2:end]))
)

sort!(df, :var)

scatter(
    df.mean,
    df.var,
    xlabel="Mean colony size per condition",
    ylabel="Variance",
    legend=false,
    c=:black,
    α=.25,
    size=(400,300),
    xlims=(-0.4, 0.4),
)

for r = eachrow(@where(df, :var .> 1.))
    plot!(annotation = (r.mean, r.var, text(string(r.condition), 5, :right, :top)))
end

for r = eachrow(@where(df, :mean .< -.2))
    pos = r.condition == :YES_Diamide_2mM ? :right : :hcenter
    plot!(annotation = [(r.mean, r.var, text(string(r.condition), 5, pos, :top))])
end

plot!()

# using GLM
#
# m = lm( @formula(var ~ mean), df)
#
# plot!(
#     df.mean,
#     predict(m, df),
# )

savefig(joinpath(@__DIR__, "colony_sizes_mean_var.pdf"))


filter(x->occursin("EGTA", x), string.(df.condition))
filter(x->occursin("SDS", x), string.(df.condition))
filter(x->occursin("Diamide", x), string.(df.condition))


df[map(x->occursin("SDS", x), string.(df.condition)), :]
df[map(x->occursin("Diamide", x), string.(df.condition)), :]

#=
# Histograms

df = load(GrowthPhenotypes)

opts = (ylabel="Frequency", c=:grey50, legend=false)

for c = ["YES_32", "EMM_32"]
    df_ = @where(df, :condition .== c)

    tmp = by(df_, :id, mean = :size => mean)
    histogram(tmp.mean; xlabel="Size (mean)", opts...)
    savefig(joinpath(@__DIR__, "size_mean_$(c).pdf"))

    tmp = by(df_, :id, var = :size => var)
    histogram(tmp[:,:var]; xlabel="Size (variance)", opts...)
    savefig(joinpath(@__DIR__, "size_var_$(c).pdf"))

    tmp = by(df_, :id, diff = :size => x->maximum(x)-minimum(x))
    histogram(tmp[:,:diff]; xlabel="Size (max - min)", opts...)
    savefig(joinpath(@__DIR__, "size_diff_$(c).pdf"))
end
=#

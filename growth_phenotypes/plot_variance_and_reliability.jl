using PombeAgeingGenes, Plots, DataFrames, DataFramesMeta, Combinatorics, StatsBase, Statistics, Random, GLM

plotlyjs()

# Plot histogram of colony size variance across all strains and conditions

df = load(GrowthPhenotypesNoOutliers)

df_var = by(df, [:id, :condition], var = :size => var)
filter!(r->!isnan(r.var), df_var) # strain-condition pairs with 1 repeat have var = NaN

mean(df_var.var)
quantile(df_var.var, .95)

histogram(
    filter(!isnan, df_var.var),
    xlabel="Colony size variance",
    ylabel="Strain-condition pairs",
    yscale=:log10,
    legend=false,
    size=(400,300),
)

savefig(joinpath(@__DIR__, "colony_sizes_variance.pdf"))


# Plot scatter of reproducibility between repeats in one condition

# function reliability(df)
#     xs = Float64[]
#     ys = Float64[]
#     for g in groupby(df, :id)
#         for (x, y) in combinations(g.size, 2)
#             push!(xs, x)
#             push!(ys, y)
#         end
#     end
#     return xs, ys
# end

function reliability(df)
    xs = Float64[]
    ys = Float64[]
    for g in groupby(df, :id)
        for xy in combinations(g.size, 2)
            shuffle!(xy)
            push!(xs, xy[1])
            push!(ys, xy[2])
        end
    end
    return xs, ys
end

function plot_reliability(xs, ys, n=1_000)
    idxs = sample(1:length(xs), n, replace=false)
    x = xs[idxs]
    y = ys[idxs]
    fig = scatter(
        x,
        y,
        xlabel="Colony size in repeat x",
        ylabel="Colony size in repeat y",
        c=:black,
        α=.1,
        lw=1,
        markersize=2,
        label="",
        size=(400,400),
        legend = false,
    )
    # y = x line
    plot!(
        collect(extrema(x)),
        identity,
        c=:black,
        label="y = x",
    )
    # linreg
    m = fit(LinearModel, x[:,:], y)
    @show "r² = $(round(r²(m), digits=3))"
    # plot!(
    #     x,
    #     # predict(m, df),
    #     predict(m, x[:,:]),
    #     c=:red,
    #     label="r² = $(round(r²(m), digits=3))",
    # )
    return fig
end

df = load(GrowthPhenotypesNoOutliers)

colonies_per_condition = by(df, :condition) do g
    size(g, 1)
end

sort!(colonies_per_condition, :x1)

# condition = "YES_glycerol"
condition = "YES_SDS_0.04percent"

Random.seed!(0)
xs, ys = reliability(@where(df, :condition .== condition))
plot_reliability(xs, ys, 20_000)
savefig(joinpath(@__DIR__, "colony_sizes_reliability_$(condition).pdf"))


n = 10000
idxs = sample(1:length(xs), n, replace=false)
x = xs[idxs]
y = ys[idxs]

cor(x, y)
cor(xs, ys)

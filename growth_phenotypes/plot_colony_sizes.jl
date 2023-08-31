#=
Plots cumulative frequency/probability distributions of colony sizes.
=#

using PombeAgeingGenes, Plots, DataFrames, Distributions

plotlyjs()

df = load(GrowthPhenotypesWideform)

A = Matrix(df[:,Not(:id)])
v = sort(vec(A))

N = fit(Normal, v)
std(N)^2


fig = stephist(
    v,
    xlabel="Colony size",
    ylabel="Frequency",
    bins=1000,
    legend=false,
    c=:black,
    size=(300, 300),
    right_margin=5Plots.PlotMeasures.mm,
)

##
# Add a second subplot that is zoomed in on the x-axis

fig2 = deepcopy(fig)

plot!(
    fig2,
    xlims=(-.5, .5),
)

plot(
    fig,
    fig2,
    size=(600, 300),
)

##

savefig(joinpath(@__DIR__, "colony_sizes_histogram.pdf"))

# Cum hist of colonies
plot(v, 1:length(v), xlabel="log2 relative size", ylabel="Cumulative frequency (colonies)", legend=false)

savefig(joinpath(@__DIR__, "colony_sizes_cumulative_histogram.pdf"))

#=
# Cum hist of colonies per id or condition

function _cum_sizes(df, by, f)
    plot()
    for g = groupby(df, by)
        x = sort(g[:size])
        if length(x) > 1
            plot!(x, f(1:length(x)), α=.5, c=:black)
        end
    end
    plot!(xlabel="Relative size per $(string(by))", ylabel="Cumulative probability", legend=false)
end
cumfreq_sizes(df, by) = _cum_sizes(df, by, identity)
cumprob_sizes(df, by) = _cum_sizes(df, by, x->(x ./ length(x)))

# Per condition
cumfreq_sizes(df, :condition)
savefig("colony_sizes_cumulative_frequency_per_condition.pdf")

cumprob_sizes(df, :condition)
savefig("colony_sizes_cumulative_probability_per_condition.pdf")

names(df)

# Per id

cumfreq_sizes(df, :id)
savefig("colony_sizes_cumulative_frequency_per_id.pdf")

cumprob_sizes(df, :id)
savefig("colony_sizes_cumulative_probability_per_id.pdf")
=#

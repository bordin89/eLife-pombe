#=
Bar plot of number of colonies per processing step.
=#

using JSON, Plots, Format

plotlyjs()

d = JSON.parsefile(joinpath(ENV["POMBEAGEINGGENES"], "data", "ncolonies.json"))

x = [
    "All",
    "Remove 'grid' and 'empty'",
    "Remove poor-quality colonies",
]

y = map(x->d[x], x)

labels = [(i-.5, y[i]+150_000, text(format(y[i], commas=true), 7)) for i = 1:length(x)]

bar(
    [
        "All data",
        "All colonies", # "Remove 'grid' and 'empty'",
        "High-quality colonies",# "Remove poor-quality colonies",
    ],
    y,
    # xticks = (0:3, ),
    xlabel="Growth phenotype data",
    ylabel="Colonies",
    c=:grey80,
    lw=0,
    # xrot=25,
    legend=false,
    # annotations=labels,
    grid=:y,
    # bottom_margin=-10Plots.PlotMeasures.mm,
    ylim=(0, 3_100_000),
    size=(400, 300),
)

savefig(joinpath(@__DIR__, "processing_steps.pdf"))

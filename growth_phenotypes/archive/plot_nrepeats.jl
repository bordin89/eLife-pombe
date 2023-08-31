include("growth_phenotypes_0_common.jl")

using PombeAgeingGenes, Plots

df = deserialize("$(fname)_no_outliers.jld")

nrepeats = by(df, [:condition, :id], nrepeats = :size => length)

# Use extrema to work out xlims for log10 xscale
extrema(nrepeats[:nrepeats])

histogram(nrepeats[:nrepeats],
    xlabel="Number of repeats per strain-condition pair",
    ylabel="Frequency",
    xscale=:log10,
    xlims=(.75, 2500),
    yscale=:log10,
    legend=false,
    c=:grey50,
    )

savefig("nrepeats_histogram.pdf")

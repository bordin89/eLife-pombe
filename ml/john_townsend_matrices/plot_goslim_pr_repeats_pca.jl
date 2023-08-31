#=
Plot PR curves of multiple repeats of GO Slim predictions.
=#

using PombeAgeingGenes, Plots

gr()
# plotlyjs() # not working

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts/ml/john_townsend_matrices/results/pca")

xs = filter(x->startswith(x, "x"), readdir(dir))

function loaddata(matrix_name)
    d = Dict()
    for repeat = readdir(joinpath(dir, matrix_name))
        !isdir(joinpath(dir, matrix_name, repeat)) && continue
        d[repeat] = Dict()
        # Load results for each GO term
        for file = readdir(joinpath(dir, matrix_name, repeat))
            !endswith(file, ".json") && continue
            goterm = replace(file, ".json"=>"")
            cv = loadcvresults(joinpath(dir,matrix_name,repeat,goterm*".json"))
            d[repeat][goterm] = (ŷs = vcat(cv["ŷs"]...),
                                 ys = vcat(cv["ys"]...))
        end
    end
    return d
end

function plot_PR_repeat(d)
    fig = plot()
    ŷss = []
    yss = []
    opts = (baseline=false, c=:grey70)

    # micro-averaged PR curves per repeat
    for r = keys(d)
        ŷs = vcat(getfield.(values(d[r]), :ŷs)...)
        ys = vcat(getfield.(values(d[r]), :ys)...)
        plot!(PR(ŷs, ys); label="", α=.25, lw=.5, opts...)
        push!(ŷss, ŷs)
        push!(yss, ys)
    end

    # micro-averaged PR curve of all repeats
    plot!(PR(vcat(ŷss...), vcat(yss...)); label="", opts...)

    plot!(size=(450,300), left_margin=10Plots.PlotMeasures.mm)
    fig
end

figs = Plots.Plot[]
for x in xs
    d = loaddata(x)
    fig = plot_PR_repeat(d)
    title!(x)
    push!(figs, fig)
end

plot(figs...,
    # layout=(3,4),
    # size=(1000,600),
    right_margin=5Plots.PlotMeasures.mm,
    )

savefig(joinpath(@__DIR__, "pca_pr_curves.png"))

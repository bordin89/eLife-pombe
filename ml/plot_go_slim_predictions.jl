#=
Plot PR curves of multiple repeats of GO Slim predictions.
=#

using PombeAgeingGenes
using Plots

c = PlotThemes.wong_palette;

function loaddata(dir)
    d = Dict()
    # dir = joinpath(dir), "predictions")
    # Load results for each GO term
    for file = readdir(dir)
        !endswith(file, ".json") && continue
        goterm = replace(file, ".json"=>"")
        cv = loadcvresults(joinpath(dir, goterm*".json"))
        d[goterm] = (ŷs = vcat(cv["ŷs"]...), ys = vcat(cv["ys"]...))
    end
    d
end

function plot_PR(d, dir)
    dir = joinpath(dir, "pr")
    ispath(dir) || mkpath(dir)

    fig = plot()
    opts = (baseline=false,)

    ŷs = getfield.(values(d), :ŷs)
    ys = getfield.(values(d), :ys)

    # micro-averaged PR curve of all repeats
    plot!(PR(vcat(ŷs...), vcat(ys...)); opts...) # label=label

    plot!(size=(450,300), left_margin=10Plots.PlotMeasures.mm)

    savefig(joinpath(dir, "micro-averaged_PR.pdf"))
    savefig(joinpath(dir, "micro-averaged_PR.png"))
    fig
end

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "go_slim", "RandomForestClassifier", "prediction", "ne_ff_cor2")

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "go_slim", "RandomForestClassifier", "prediction", "gp_ne_cor2")

data = loaddata(dir)

plot_PR(data, dir)





#=
series = Dict(
    "ne_ff" => "NE_FF",
    # "ne" => "NE",
    # "ff" => "FF",
    # "gp" => "GP",
    # "gp_ne" => "GP_NE",
    )

function loaddata()
    d = Dict()
    # Load results for each set of features
    for features = keys(series)
        d[features] = Dict()
        # Load results for each repeat
        for repeat = readdir(joinpath(dir, features))
            !isdir(joinpath(dir, features, repeat)) && continue
            d[features][repeat] = Dict()
            # Load results for each GO term
            for file = readdir(joinpath(dir, features, repeat))
                !endswith(file, ".json") && continue
                goterm = replace(file, ".json"=>"")
                # @show goterm
                cv = loadcvresults(joinpath(dir,features,repeat,goterm*".json"))
                d[features][repeat][goterm] = (ŷs = vcat(cv["ŷs"]...),
                                               ys = vcat(cv["ys"]...))
            end
        end
    end
    d
end

function plot_PR_repeat(data, dir)
    fig = plot()

    for (i, (features, label)) = enumerate(series)
        ŷss = []
        yss = []
        d = data[features]
        opts = (baseline=false, c=c[i])

        # micro-averaged PR curves per repeat
        for r = keys(d)
            ŷs = vcat(getfield.(values(d[r]), :ŷs)...)
            ys = vcat(getfield.(values(d[r]), :ys)...)
            plot!(PR(ŷs, ys); label="", α=.25, lw=.5, opts...)
            push!(ŷss, ŷs)
            push!(yss, ys)
        end

        # micro-averaged PR curve of all repeats
        plot!(PR(vcat(ŷss...), vcat(yss...)); label=label, opts...)
    end

    plot!(size=(450,300), left_margin=10Plots.PlotMeasures.mm)
    savefig(joinpath(dir, "pr", "micro-averaged_PR.pdf"))
    savefig(joinpath(dir, "pr", "micro-averaged_PR.png"))
    fig
end

=#

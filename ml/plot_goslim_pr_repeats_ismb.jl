#=
Plot PR curves of multiple repeats of GO Slim predictions.
=#

using PombeAgeingGenes
using Plots
plotlyjs()
gr()

c = PlotThemes.wong_palette;

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
        opts = (baseline=false, c=c[i], resolution=.001)

        # micro-averaged PR curves per repeat
        for r = keys(d)
            ŷs = vcat(getfield.(values(d[r]), :ŷs)...)
            ys = vcat(getfield.(values(d[r]), :ys)...)
            # plot!(PR(ŷs, ys); label="", α=.25, lw=.5, opts...)
            push!(ŷss, ŷs)
            push!(yss, ys)
        end

        # micro-averaged PR curve of all repeats
        plot!(PR(vcat(ŷss...), vcat(yss...)); label=label, aspect_ratio=:equal, opts...)
    end

    plot!(
        # size=(500,400),
        # left_margin=10Plots.PlotMeasures.mm,
    )
    savefig(joinpath(dir, "pr", "micro-averaged_PR_ismb.pdf"))
    savefig(joinpath(dir, "pr", "micro-averaged_PR_ismb.png"))
    fig
end

using DataStructures

series = OrderedDict(
    "ne" => "Networks",
    "ff" => "FunFams",
    "ne_ff" => "Networks and FunFams",
    # "gp" => "GP",
    # "gp_ne" => "GP_NE",
    )

model = "RandomForestClassifier"
commit = "final"
dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts","ml","go_slim",model,commit)
let p = joinpath(dir, "pr")
    ispath(p) || mkpath(p)
end

data = loaddata()

plotlyjs()
plot_PR_repeat(data, dir)

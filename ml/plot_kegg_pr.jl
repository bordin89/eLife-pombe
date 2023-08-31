#=
Plot PR curves of GO Slim predictions.
=#

using PombeAgeingGenes
using Plots
plotlyjs()

c = PlotThemes.wong_palette
c = c[[1, 3, 5]]

iter = enumerate([["gp", "Growth phenotypes"],
                  ["ne", "Network embeddings"],
                  ["gp_ne", "Both"]])

function loaddata()
    d = Dict()
    for features = ["gp", "ne", "gp_ne"]
        d[features] = Dict()
        for f = readdir(joinpath(dir, features))
            !endswith(f, ".json") && continue
            goterm = replace(f, ".json"=>"")
            d2 = loadcvresults("$dir/$features/$goterm.json")

            d[features][goterm] = (ŷs = vcat(d2["ŷs"]...), ys = vcat(d2["ys"]...))
        end
    end
    d
end

model = "RandomForestClassifier"
commit = "8e62edf"
# goterm = "GO0002181"
dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts","ml","kegg",model,commit)

data = loaddata()

# violin plots

fig = plot()
for (i, (features, label)) = iter
    aucs = []
    for (goterm, y) = data[features]
        push!(aucs, auc(PR(y.ŷs, y.ys)))
    end
    violin!(aucs, c=nothing)
end
plot!(legend=false, ylabel="AUPR",
    xticks=(1:3, ["Growth phenotypes", "Network embeddings", "Both"]))
savefig(joinpath(dir, "pr", "aupr_violin.pdf"))


# individual PR curves

for pathway = keys(data["gp"])
    fig = plot()
    for (i, (features, label)) = iter
        haskey(data[features], pathway) || continue
        y = data[features][pathway]
        pr = PR(y.ŷs, y.ys)
        plot!(pr, baseline=i == 1, label=label, c=c[i])
    end
    plot!(legend=true, size=(450,300), left_margin=10Plots.PlotMeasures.mm)
    savefig(joinpath(dir, "pr", pathway*".pdf"))
end

# micro-averaged PR curves

fig = plot()
for (i, (features, label)) = iter
    ŷs = vcat(getfield.(values(data[features]), :ŷs)...)
    ys = vcat(getfield.(values(data[features]), :ys)...)
    plot!(PR(ŷs, ys), baseline=false, c=c[i], label=label)
end
plot!(size=(450,300), left_margin=10Plots.PlotMeasures.mm)
savefig(joinpath(dir, "pr", "micro-averaged_PR.pdf"))

# micro-averaged PR curves and individual curves

fig = plot()
for (i, (features, label)) = iter
    for (goterm, y) = data[features]
        pr = PR(y.ŷs, y.ys)
        plot!(pr, baseline=false, label="", c=c[i], α=.1)
    end

    ŷs = vcat(getfield.(values(data[features]), :ŷs)...)
    ys = vcat(getfield.(values(data[features]), :ys)...)
    plot!(PR(ŷs, ys), baseline=false, c=c[i], label=label)
end
plot!(size=(450,300), left_margin=10Plots.PlotMeasures.mm)
savefig(joinpath(dir, "pr", "micro-averaged_individual_PR.pdf"))

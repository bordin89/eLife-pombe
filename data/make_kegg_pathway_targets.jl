#=
Make targets from genes in KEGG pathways for function prediction tasks.
=#

using PombeAgeingGenes, DataFrames, CSVFiles

df = DataFrame(id=String[], pathway=String[])

for pathway = KEGGPATHWAYS
    ids = keys(load(KEGGPathway(pathway)))

    for id = ids
        push!(df, (id, pathway))
    end
end

df[:value] = 1

targets = unstack(df, :id, :pathway, :value)
foreach(col->targets[col] = coalesce.(targets[col], 0), names(targets))

colwise(sum, targets[2:end])

save(joinpath(ENV["POMBEAGEINGGENES"], "data", "keggpathway_targets.csv"), targets)

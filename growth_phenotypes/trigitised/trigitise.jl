#=
Export trigitised wideform growth phenotypes for Maria's paper.
=#

using PombeAgeingGenes
using Plots
using DataFrames

df = load(GrowthPhenotypesWideform)

function f(df, r)
    df = trigitise!(deepcopy(df), ratio=r)
    A = Matrix(df[:,Not(:id)])
    clusteredheatmap(A, c=:RdBu, clims=(-2,2))
    plot!(xlabel="Conditions", ylabel="Genes")
    savefig(joinpath(@__DIR__, "growth_phenotypes_trigitised_$(Int(r*100))_pc.png"))
end

for r = 0.1:0.1:0.3
    f(df, r)
end

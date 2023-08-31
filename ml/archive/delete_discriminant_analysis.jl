#=
https://www.researchgate.net/publication/308015273_Linear_vs_quadratic_discriminant_analysis_classifier_a_tutorial

http://users.isr.ist.utl.pt/~wurmd/Livros/school/Bishop%20-%20Pattern%20Recognition%20And%20Machine%20Learning%20-%20Springer%20%202006.pdf
=#
cd(@__DIR__)

using PombeAgeingGenes, PombeAgeingGenes.GeneOntology, MultivariateStats, DataFrames, Statistics, DataFramesMeta, HypothesisTests, RDatasets, Random, MLBase, SparseArrays

using Plots, StatPlots

## Yeast

d = load(GrowthPhenotypesNoOutliers)
d = PombeAgeingGenes.wideform(d)
dropmissing!(d)

# Standardize
for col = names(d)[2:end]
    d[col] = PombeAgeingGenes.standardize(d[col])
end

# ageinggoslimterms = [:GO0006915, :GO0006914, :GO0005975, :GO0006325, :GO0006281, :GO0032200]
# goterm = ageinggoslimterms[6]

targets = load(GOSlimTargets)

commonids = intersect(d[:id], targets[:id])
X = @in(d, :id, commonids)
Y = @in(targets, :id, commonids)
sort!(X, :id)
sort!(Y, :id)

X = collect(convert(Array, X[2:end])')
Y = collect(convert(Array{Int}, Y[2:end])')
goterms = names(targets)[2:end]

plotlyjs()

# Run FisherLinearDiscriminant CV on all GO Slim terms
for (i, goterm) = enumerate(goterms)
    y = Y[i,:]

    ps, prs = fit(CrossValidation(5, 10), FisherLinearDiscriminant, X, y)
    mean(ps)

    fig = @df DataFrame(ps) violin(:metrics, :values, ylim=(0,1), legend=false, title=goterm)
    savefig("FisherLinearDiscriminant/$(goterm)_metrics.pdf")

    fig = @plotpr prs
    plot!(title=goterm)
    savefig("FisherLinearDiscriminant/$(goterm)_PR.pdf")
    # break
end

# Run KernelFisherDiscriminant CV on all GO Slim terms
# KernelFisherDiscriminant works really well on Iris data set, but not here

K = PombeAgeingGenes.Kernel(X)

for (i, goterm) = enumerate(goterms)
    i != 5 && continue

    println(i, " ", goterm)
    println("Fitting")

    y = Y[i,:]

    ps, prs = fit(CrossValidation(5, 10), KernelFisherDiscriminant, X, y, K=K)

    println("Plotting")

    # fig = @df DataFrame(ps) violin(:metrics, :values, ylim=(0,1), legend=false, title=goterm)
    # savefig("KernelFisherDiscriminant/$(goterm)_metrics.pdf")

    fig = @plotpr prs
    plot!(title=goterm)
    savefig("KernelFisherDiscriminant/$(goterm)_PR.pdf")
    # break
end

## Iris
using RDatasets, SparseArrays, MLBase, Random, PombeAgeingGenes, Statistics

iris = dataset("datasets", "iris")
indices = collect(1:size(iris, 1))
shuffle!(indices)

Y = string.(iris[5])[indices]
Y = labelencode(labelmap(Y), Y)
Y = Array{Int}(sparse(1:length(Y), Y, true))
yᵢ = 3

y = Y[:,yᵢ]

X = convert(Array, iris[1:4])[indices,:]'
X = PombeAgeingGenes.standardize(X; dims=2)

ps, prs = fit(CrossValidation(5,10), FisherLinearDiscriminant, X, y)
mean(ps)
fig = @plotpr prs
savefig("iris_fda.pdf")

K = Kernel(X)

ps, prs = fit(CrossValidation(5,10), KernelFisherDiscriminant, X, y, K=K)
mean(ps)
fig = @plotpr prs
savefig("iris_kfda.pdf")

M = fit(KernelFisherDiscriminant, X, y)
evaluate(M, X)
probability(M, X)
y
using StatPlots
d = DataFrame(ps)
@df d violin(:metrics, :values, ylim=(0,1))

cd(@__DIR__)
using PombeAgeingGenes
using Lasso, MLBase
using DataFrames, RDatasets, Random, MLBase, SparseArrays, Statistics

# standardize_(A::AbstractArray; dims=:) = (A .- mean(A, dims=dims)) ./ std(A, dims=dims)

function classweights_(y)
    yp = length(y) / (2 * count(y .== 1))
    yn = length(y) / (2 * count(y .== 0))
    map(x->x == 1 ? yp : yn, y)
end

function Lasso.cross_validate_path(path::RegularizationPath; gen=Kfold(length(y),10),
                                   select=:CVmin, fitargs...)
    m = path.m
    y = m.rr.y
    offset = m.rr.offset
    Xstandardized = m.pp.X
    cross_validate_path(path, Xstandardized, y; gen=gen, select=select, offset=offset,
                        standardize=false, fitargs...)
end

# function kthfold(X::AbstractArray, Y::AbstractArray, trainidxs::AbstractArray)
#     testidxs = setdiff(1:size(Y, 1), trainidxs)
#     x, tx = X[trainidxs,:], X[testidxs,:]
#     y, ty = Y[trainidxs], Y[testidxs]
#
#     M = fit(LassoPath, x, y, Binomial(), irls_maxiter=1000)
#     M̂ = cross_validate_path(M; gen=StratifiedKfold(y, 10), irls_maxiter=1000)
#     tŷ = predict(M, tx)[:,M̂]
#     tŷ, ty
# end
#
# function cv(X::AbstractArray, Y::AbstractVector, k::Integer=10)
#     ŷs = Vector{Float64}[]; ys = Vector{Float64}[]
#     gen = collect(StratifiedKfold(Y, k))
#
#     Threads.@threads for i = 1:k
#     # for i = 1:k
#         trainidxs = gen[i]
#         tŷ, ty = kthfold(X, Y, trainidxs)
#         push!(ŷs, tŷ); push!(ys, ty)
#     end
#
#     ŷs, ys
# end

function kthfold(X::AbstractArray, tX::AbstractArray, y::AbstractVector, ty::AbstractVector)
    M = fit(LassoPath, X, y, Binomial(), irls_maxiter=1000)
    M̂ = cross_validate_path(M; gen=StratifiedKfold(y, 10), irls_maxiter=1000)
    tŷ = predict(M, tX)[:,M̂]
    tŷ, ty
end

function cv(X::AbstractArray, y::AbstractVector, k::Integer=10)
    ŷs = Vector{Float64}[]
    ys = Vector{Float64}[]
    gen = collect(StratifiedKfold(Y, k))

    for i = 1:k
        trainidxs = gen[i]
        tŷ, ty = kthfold(traintestsplit(X, trainidxs)..., traintestsplit(y, trainidxs)...)
        push!(ŷs, tŷ)
        push!(ys, ty)
    end

    ŷs, ys
end

## Iris
iris = dataset("datasets", "iris")
indices = collect(1:size(iris, 1))
shuffle!(indices)

Y = string.(iris[5])[indices]
Y = labelencode(labelmap(Y), Y)
Y = Array{Float64}(sparse(1:length(Y), Y, true))
yᵢ = 3
Y = Y[:,yᵢ]

X = convert(Array, iris[1:4])[indices, :]
X = standardize(X; dims=1)

k = 10

ŷs, ys = cv(X, Y)

performances = [Performance(ŷ,y) for (ŷ,y) = zip(ŷs,ys)]
mean(accuracy, performances)
mean(precision, performances)
mean(PombeAgeingGenes.recall, performances)
mean(f1, performances)
mean(mcc, performances)
# mean(PombeAgeingGenes.fisher, performances)

pvalues = [NullDistribution((ŷ,y)->PR(ŷ,y).auc, ŷ, y)(PR(ŷ,y).auc) for (ŷ,y) = zip(ŷs,ys)]

using Plots, StatPlots

histogram(pvalues)

prs = [PR(ŷ, y, 1.) for (ŷ, y) = zip(ŷs, ys)]
mean(auc, prs)
PombeAgeingGenes.@plotprs prs


## Yeast
d = load(GrowthPhenotypesNoOutliers)
d = wideform(d)
dropmissing!(d)

Y = load(GeneOntology.GOSlimTargets)
for c = names(Y)
    Y[c] = coalesce.(Y[c], 0)
end

ageinggoslimterms = [:GO0006915, :GO0006914, :GO0005975, :GO0006325, :GO0006281, :GO0032200]
Y = Y[[:id; ageinggoslimterms]]

map(sum, eachcol(Y[2:end]))

commonids = intersect(d[:id], Y[:id])
X = @in(d, :id, commonids)
Y = @in(Y, :id, commonids)
sort!(X, :id)
sort!(Y, :id)

X = convert(Array{Float64}, X[2:end])
Y = convert(Array{Float64}, Y[2:end])

X = standardize_(X; dims=1)

indices = collect(1:size(X, 1))
shuffle!(indices)

X = X[indices,:]
Y = Y[indices,:]

Y_ = deepcopy(Y)
Y = Y_[:,6]

using Plots, LassoPlot
M = fit(LassoPath, X, Y, Binomial(), irls_maxiter=100)
plot(M, select=:AICc, showselectors=[:AICc], nCVfolds=100)
savefig("example_regularization_paths.pdf")

ŷs, ys = cv(X, Y)

performances = [Performance(ŷ,y) for (ŷ,y) = zip(ŷs,ys)]
mean(performances)

pvalues = [NullDistribution((ŷ,y)->PR(ŷ,y).auc, ŷ, y)(PR(ŷ,y).auc) for (ŷ,y) = zip(ŷs,ys)]

using Plots, StatPlots

histogram(pvalues)

prs = [PR(ŷ, y, 1.) for (ŷ, y) = zip(ŷs, ys)]
mean(auc, prs)
PombeAgeingGenes.@plotprs prs

using StatPlots


function dfforplotting(ps)
    n = length(ps)
    df = DataFrame(
        metric = [fill("Accuracy", n); fill("Precision", n); fill("Recall", n);
            fill("f1", n); fill("MCC", n)],
        value = [accuracy.(ps); precision.(ps); PombeAgeingGenes.recall.(ps); f1.(ps); mcc.(ps)])
end

tmp = dfforplotting(performances)

@df tmp violin(:metric, :value, legend=false, c=:grey75, tickfontsize=10)
savefig("example_performance.pdf")

using Plots
plotlyjs()
prs = [PR(ŷ, y, 1.) for (ŷ, y) = zip(ŷs, ys)]
mean(auc, prs)

# @plotpr PR(vcat(ŷs...), vcat(ys...), 1.)
PombeAgeingGenes.@plotprs prs
savefig("example_PR.pdf")

#=
Logisitc regression with structured multi-class prediction on Iris data set

https://stackoverflow.com/questions/49135107/logistic-regression-using-flux-jl
=#

using PombeAgeingGenes
using RDatasets, Statistics, Random, SparseArrays, LinearAlgebra, DataFrames, Plots
using Flux
using Flux: onehotbatch, onecold, throttle, @epochs, binarycrossentropy, crossentropy
using Base.Iterators: partition

###########################################################################################

function onehot(xs)
    labels = unique(xs)
    d = Dict(zip(labels, 1:length(labels)))
    ys = map(x->d[x], xs)
    sparse(1:length(ys), ys, true)
end

vecofvecs(xs; dims=2) = [xs[:,i] for i = 1:size(xs, dims)]

"""
    logloss_(ŷ, y)

Results are the same as `sklearn.metrics.log_loss`.
"""
function logloss_(ŷ::AbstractArray, y::AbstractArray; weight=1)
    mean(binarycrossentropy.(ŷ, y) .* weight)
end

l1(x) = sum(y->norm(y, 1), x)
l2(x) = sum(y->norm(y, 2)^2, x)

function loss(x, y, λ=.01)
    ŷ = m(x)
    penalty = λ * l1(params(m))
    logloss_(ŷ, y; weight=classweights_(y)) + penalty
    # logloss_(ŷ, y) + penalty
end

# function loss(x, y, λ=.005)
#     penalty = λ * l1(params(m))
#     crossentropy(m(x), y, weight=classweights_(y)) + penalty
# end

mutable struct History
    train_loss::Vector{Float64}
    train_accuracy::Vector{Float64}
    train_precision::Vector{Float64}
    val_loss::Vector{Float64}
    val_accuracy::Vector{Float64}
    val_precision::Vector{Float64}
end

History() = History([Float64[] for _ = fieldnames(History)]...)

function trainingperformance!(history::History)
    trainingperformance!(history, "train", x, y)
    trainingperformance!(history, "val", Tx, Ty)
end

function trainingperformance!(history::History, mode, x, y)
    ŷ = m(x)
    push!(getfield(history, Symbol("$(mode)_loss")), loss(x, y).data)
    push!(getfield(history, Symbol("$(mode)_accuracy")), accuracy(ŷ, y))
    push!(getfield(history, Symbol("$(mode)_precision")), precision(ŷ, y))
end

macro train(n, history, ex)
    esc(quote
        @progress for i = 1:$(n)
            @info "Epoch $i"
            $ex
            trainingperformance!($history)
        end
    end)
end

function plotperformance(history)
    fig = plot()

    for x = ["train", "val"], y = ["loss", "accuracy"]
        l = "$(x)_$(y)"
        plot!(getfield(history, Symbol(l)), label=l)
    end

    plot!(xlabel="Epochs")
end

function classweights_(y)
    n_samples = length(y)
    p = count(y .== 1)
    n = n_samples - p
    n_classes = 2
    p = n_samples / (n_classes * p)
    n = n_samples / (n_classes * n)
    [i == 1 ? p : n for i = y]
end

###########################################################################################

# accuracy(x, y) = mean(onecold(m(x)) .== onecold(y))

# # Iris
# iris = dataset("datasets", "iris")
# indices = collect(1:size(iris, 1))
# shuffle!(indices)
#
# x = collect(convert(Array, iris[1:4]))[indices, :]'
# y = onehot(string.(iris[5])[indices])'
#
# x, Tx = x[:,1:125], x[:,126:end]
# y, Ty = y[:,1:125], y[:,126:end]
#
# data = zip(vecofvecs(x), vecofvecs(y))
# m = Chain(Dense(4, 3, sigmoid))

# Yeast
df = load(GrowthPhenotypesNoOutliers)
df = wideform(df)
df = dropmissing(df)

Y = load(GeneOntology.GOSlimTargets)
for c = names(Y)
    Y[c] = coalesce.(Y[c], 0)
end

ageinggoslimterms = [:GO0006915, :GO0006914, :GO0005975, :GO0006325, :GO0006281, :GO0032200]
Y = Y[[:id; ageinggoslimterms]]

commonids = intersect(df[:id], Y[:id])
X = @in(df, :id, commonids)
Y = @in(Y, :id, commonids)
sort!(X, :id)
sort!(Y, :id)

X = collect(convert(Array{Float64}, X[2:end])')
Y = collect(convert(Array{Float64}, Y[2:end])')

X = standardize(X; dims=2)

indices = collect(1:size(X, 2))
shuffle!(indices)

X = X[:, indices]
Y = Y[:, indices]

tvsplit = 2000
nth = 3

x, Tx = X[:, 1:tvsplit], X[:, tvsplit+1:end]
y, Ty = Y[:, 1:tvsplit], Y[:, tvsplit+1:end]
y, Ty = y[nth:nth, :], Ty[nth:nth, :]

data = zip(vecofvecs(x), vecofvecs(y))

# m = Chain(Dense(size(x, 1), 32, relu),
#           Dense(32, size(y, 1), sigmoid))

m = Chain(Dense(size(x, 1), size(y, 1), sigmoid))

# Train
m(x)

accuracy(m(x), y)

evalcb = () -> @show(loss(x, y))
opt = ADAM(params(m))

history = History()

# a = @epochs 2 Flux.train!(loss, data, opt, cb = throttle(evalcb, 10))
@train 10 history Flux.train!(loss, data, opt, cb = throttle(evalcb, 10))

accuracy(m(x), y)

f1(m(x), y)

plotperformance(history)

ŷ = m(x)

precision(ŷ, y)

sum(y .== 1)

Flux.crossentropy(ŷ, y)
Performance(ŷ, Array{Bool}(y))

pr = PR(ŷ, Array{Bool}(y), true)
@plotpr(pr)

pr = PR(m(Tx), Bool.(Ty), true)
@plotpr(pr)


using Lasso, LassoPlot, Statistics

standardize(A; dims=:) = (A .- mean(A, dims=dims)) ./ std(A, dims=dims)





using MultivariateStats

Xp, Xn = X[:, Y[nth,:] .== 1], X[:, Y[nth,:] .== 0]
M = fit(LinearDiscriminant, Xp, Xn)

Yn = map(i->evaluate(M, Xn[:,i]), 1:size(Xn,2))
Yp = map(i->evaluate(M, Xp[:,i]), 1:size(Xp,2))

Ŷ = map(i->evaluate(M, X[:,i]), 1:size(X,2))

using StatPlots
histogram(Yn, label="N", xlabel="Score", ylabel="Freq")
histogram!(Yp, label="P")

Performance(Ŷ, Bool.(Y[nth,:]))

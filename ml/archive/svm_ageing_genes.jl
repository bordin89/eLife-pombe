using
    Distributed,
    PombeAgeingGenes,
    StatsBase,
    MLDataUtils,
    LIBSVM,
    DelimitedFiles

function hyperparameteroptim(X, y, param_grid::Dict; kwargs...)
    best_score = 0.
    best_params = nothing

    for gamma = param_grid[:gamma], cost = param_grid[:cost]
        params = (gamma=gamma, cost=cost)
        M = svmtrain(X, y; probability=true, params..., kwargs...)
        p, d = svmpredict(M, X)
        ŷ = d[y[1] == 1 ? 1 : 2,:]
        s = auc(PR(ŷ, y))

        if s > best_score
            best_params = params
            best_score = s
        end
    end

    best_params
end

function f(X, y; kwargs...)
    ŷs = Vector{Float64}[]
    ys = Vector{Float64}[]
    param_grid = Dict(:gamma => [0.1, 1., 10.], :cost => [1., 10., 100.])

    for ((Xtrain, ytrain), (Xtest, ytest)) = stratifiedkfolds((X, y), 5)
        params = hyperparameteroptim(Xtrain, ytrain, param_grid; kwargs...)
        @show params
        M = svmtrain(Xtrain, ytrain; cachesize=1000., params..., kwargs...)
        p, d = svmpredict(M, Xtest)
        ŷ = d[ytest[1] == 1 ? 1 : 2,:]
        push!(ŷs, ŷ)
        push!(ys, ytest)
    end

    PR(vcat(ŷs...), vcat(ys...)), Performance.(ŷs, ys)
end

using
    PombeAgeingGenes,
    StatsBase,
    MLDataUtils,
    DelimitedFiles

X, _ = load(ML, deepNF_embeddings=true, center=true)

ageing_genes = unique(vcat(values(GeneOntology.ageingterms)...))

y = float.([x in ageing_genes for x = X[:id]])

X = permutedims(convert(Matrix, X[2:end]))

@time results = f(X, y)

results = cat.(iterrange, results, dims=1)
writedlm("results.csv", results)

using Flux
using Flux: onehotbatch, onecold, crossentropy, throttle, @epochs
using Base.Iterators: repeated, partition

X = rand(50, 2500)
y = rand(Bool, 2500)
using MLDataUtils, Statistics
using RDatasets: dataset

#Classification C-SVM
iris = dataset("datasets", "iris")
Y = convert(Vector, iris[:, :Species])
Y = y
Y = onehotbatch(Y, unique(Y))
# X = permutedims(convert(Matrix, iris[:, 1:4]))
X
y
X, Y = shuffleobs((X, Y))

m = Chain(
    Dense(size(X,1), 256),
    Dense(256, 128),
    Dense(128, size(Y,1)),
    softmax)

batch_size = 128
mb_idxs = partition(1:size(X,2), batch_size)
data = [(X[:, i], Y[:, i]) for i in mb_idxs]

loss(x, y) = crossentropy(m(x), y)
accuracy(x, y) = mean(onecold(m(x)) .== onecold(y))
opt = ADAM()
evalcb = () -> @show(loss(X, Y))
@epochs 100 Flux.train!(loss, params(m), data, opt, cb = throttle(evalcb, 10))

m(X)

accuracy(X,Y)

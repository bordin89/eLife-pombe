#=
Plot bootstrap curves demonstrating how much additional AUC PR is provided by
additional growth phenotypes.
=#

using
    Distributed,
    PombeAgeingGenes,
    StatsBase,
    MLDataUtils,
    LIBSVM

addprocs(10, exeflags="--project=.")

@everywhere using
    PombeAgeingGenes,
    StatsBase,
    MLDataUtils,
    LIBSVM

using
    DelimitedFiles

@everywhere function hyperparameteroptim(X, y, param_grid::Dict; kwargs...)
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

function bootstrap_condition(X, y, nfeatures::Int, nboot::Int; kwargs...)
    @show nfeatures
    nconditions, nids = size(X)

    results = @distributed vcat for b = 1:nboot
        feature_idxs = StatsBase.sample(1:nconditions, nfeatures, replace=false)
        x = X[feature_idxs, :]
        ŷs = Vector{Float64}[]
        ys = Vector{Float64}[]
        param_grid = Dict(:gamma => [0.01, 0.1], :cost => [1., 10.])

        for ((Xtrain, ytrain), (Xtest, ytest)) = kfolds((x, y), 5)
            params = hyperparameteroptim(Xtrain, ytrain, param_grid; kwargs...)
            M = svmtrain(Xtrain, ytrain; cachesize=1000., params..., kwargs...)
            p, d = svmpredict(M, Xtest)
            ŷ = d[ytest[1] == 1 ? 1 : 2,:]
            push!(ŷs, ŷ)
            push!(ys, ytest)
        end

        ŷs, ys = vcat(ŷs...), vcat(ys...)
        auc(PR(ŷs, ys))
    end
end

function f(X, y, iterrange; nboot=10)
    map(nfeatures->bootstrap_condition(X, y, nfeatures, nboot), iterrange)
end


X, _ = load(ML)

ageing_genes = unique(vcat(values(GeneOntology.ageingterms)...))

y = Int.([x in ageing_genes for x = X[:id]])

X = permutedims(convert(Matrix, X[2:end]))

nconditions = size(X, 1)
iterrange = 1:5:nconditions-1

@time results = f(X, y, iterrange, nboot=10)

results = cat.(iterrange, results, dims=1)
writedlm("results.csv", results)


X = rand(50, 2500)
y = rand(Bool, 2500)

nconditions = size(X, 1)
iterrange = 25:5:nconditions-1

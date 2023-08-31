#=
Use growth phenotypes to predict growth phenotypes.
=#

using PombeAgeingGenes, Lasso, Statistics, Plots, Random, Distributions, StatsPlots

Random.seed!(0)

struct Error{T<:AbstractFloat}
    train::T
    test::T
    shuffle::T
    random::T
end

function Statistics.mean(errors::AbstractVector{Error})
    Error(map(x->mean(getfield.(errors, x)), [:train, :test, :shuffle, :random])...)
end

function cv(X, y, k=5)
    errors = Error[]

    for ((Xtrain, ytrain), (Xtest, ytest)) = kfolds((X, y), k)
        Xtrain, ytrain, Xtest, ytest = map(collect, (Xtrain, ytrain, Xtest, ytest))
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        # 𝒩 = fit(Normal, Xtrain)
        # Xrandom = clamp.(rand(𝒩, size(Xtest)...), 0., 2.)
        𝒩 = vec(mapslices(x->fit(Normal, x), Xtrain, dims=1))
        # Xrandom = clamp.(permutedims(hcat(map(_->rand.(𝒩), 1:size(Xtrain,1))...)), 0., 2.)

        M = fit(LassoPath, Xtrain, ytrain, Normal())
        Mi = cross_validate_path(M)

        ŷtrain = clamp.(predict(M, Xtrain)[:,Mi], 0., 2.)
        ŷtest = clamp.(predict(M, Xtest)[:,Mi], 0., 2.)
        # ŷrand = clamp.(predict(M, Xrandom)[:,Mi], 0., 2.)

        n = 10
        shuffleerrs = Vector{Float64}(undef, n)
        randomerrs = Vector{Float64}(undef, n)

        for i = 1:n
            Xrandom = clamp.(permutedims(hcat(map(_->rand.(𝒩), 1:size(Xtrain,1))...)), 0., 2.)
            # Xrandom = clamp.(rand(𝒩, size(Xtrain)...), 0., 2.)
            ŷrand = clamp.(predict(M, Xrandom)[:,Mi], 0., 2.)
            randomerrs[i] = mse(ŷrand, ytrain)
            shuffleerrs[i] = mse(ŷtrain, shuffle(ytrain))
        end

        error = Error(
            mse(ŷtrain, ytrain),
            mse(ŷtest, ytest),
            minimum(shuffleerrs),
            minimum(randomerrs))

        push!(errors, error)
    end

    mean(errors)
end

function mse(ŷs, ys)
    mean(abs2.(ŷs .- ys))
end

function loo(A)
    n = size(A,1)
    errors = Error[]

    for i = 1:n
        X = A[setdiff(1:n, i),:]
        y = A[i,:]
        error = cv(X, y)
        push!(errors, error)
    end

    errors
end

df = load(GrowthPhenotypesWideform)
A = permutedims(Matrix(df[2:end]))

n = size(A,1)
i = 30
X = A[setdiff(1:n, i),:]
y = A[i,:]
err = cv(X, y)

errs = loo(A)

mean(errs)

trainerr = getfield.(errs, :train)
testerr = getfield.(errs, :test)
shuffleerr = getfield.(errs, :shuffle)
randomerr = getfield.(errs, :random)

violin([trainerr, testerr, shuffleerr, randomerr],
    ylabel="MSE",
    xticks=(1:4, ["Train", "Test", "Shuffle", "Random"]),
    legend=false,
    )

savefig(joinpath(@__DIR__, "predict_growth_phenotypes.pdf"))

using HypothesisTests, Random, Statistics

function bootstrapttest(x::AbstractVector{T}, y::AbstractVector{T},
                        b::Integer=10^4) where T<:Real
    n = length(x)
    m = length(y)
    x̄ = mean(x)
    ȳ = mean(y)
    z̄ = mean([x; y])
    x′ = x .- x̄ .+ z̄
    y′ = y .- ȳ .+ z̄
    tt = EqualVarianceTTest(x, y).t
    tts = Vector{Float64}(undef, b)

    for i = 1:b
        tts[i] = EqualVarianceTTest(sample(x′, n), sample(y′, m)).t
    end

    p = count(tts .≥ tt) / length(tts)
end

function bootstrapconfint(x::AbstractVector{T}, y::AbstractVector{T},
                          b::Integer=10^4, q::AbstractFloat=0.95) where T<:Real
    n = length(x)
    m = length(y)
    ts = mean(x) - mean(y)
    tss = Vector{Float64}(undef, b)

    for i = 1:b
        tss[i] = mean(sample(x, n)) - mean(sample(y, m))
    end

    ts, quantile(tss, q/2), quantile(tss, 1-(q/2))
end

# Training and testing
bootstrapconfint(testerr, trainerr)
bootstrapttest(testerr, trainerr)

# Shuffled and random
bootstrapconfint(shuffleerr, randomerr)
bootstrapttest(shuffleerr, randomerr)

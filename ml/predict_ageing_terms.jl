#=
Use network embeddings to predict ageing genes.
=#

using PombeAgeingGenes, Lasso

# , Statistics, Plots, Random, Distributions, StatsPlots

# Random.seed!(0)

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
        𝒩 = vec(mapslices(x->fit(Normal, x), Xtrain, dims=1))

        # TODO why does this not work with Binomial?
        M = fit(LassoPath, Xtrain, ytrain,
            Normal(),
            # Binomial(),
            # wts=calculateweights(ytrain)),
            # irls_maxiter=10000,
            # irls_tol=1e-5,
            # cd_maxiter=10000
            )
        Mi = cross_validate_path(M)

        ŷtrain = clamp.(predict(M, Xtrain)[:,Mi], 0., 2.)
        ŷtest = clamp.(predict(M, Xtest)[:,Mi], 0., 2.)

        n = 10
        shuffleerrs = Vector{Float64}(undef, n)
        randomerrs = Vector{Float64}(undef, n)

        for i = 1:n
            Xrandom = clamp.(permutedims(hcat(map(_->rand.(𝒩), 1:size(Xtrain,1))...)), 0., 2.)
            ŷrandom = clamp.(predict(M, Xrandom)[:,Mi], 0., 2.)
            randomerrs[i] = auc(PR(ŷrandom, ytrain))
            shuffleerrs[i] = auc(PR(ŷtrain, shuffle(ytrain)))
        end

        error = Error(
            auc(PR(ŷtrain, ytrain)),
            auc(PR(ŷtest, ytest)),
            mean(shuffleerrs),
            mean(randomerrs))

        push!(errors, error)
    end

    mean(errors)
end

function mse(ŷs, ys)
    mean(abs2.(ŷs .- ys))
end

function calculateweights(y)
    n = length(y)
    np = count(y .== 1)
    nn = count(y .== 0)
    wp = np / n
    wn = nn / n
    [i == 1 ? wp : wn for i = y]
end

X, Y, _ = load(ML, Matrix, growthphenotypes=false, networkembeddings=true, Y=:ageing)
Y = vec(Y)

err = cv(X, Y)

errs = [cv(X, Y) for _ = 1:10]

trainerr = getfield.(errs, :train)
testerr = getfield.(errs, :test)
shuffleerr = getfield.(errs, :shuffle)
randomerr = getfield.(errs, :random)

violin([trainerr, testerr, shuffleerr, randomerr],
    ylabel="PR AUC",
    xticks=(1:4, ["Train", "Test", "Shuffle", "Random"]),
    legend=false,
    )

savefig(joinpath(@__DIR__, "predict_ageing_terms.pdf"))

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
bootstrapconfint(trainerr, testerr)
bootstrapttest(trainerr, testerr)

# Shuffled and random
bootstrapconfint(randomerr, shuffleerr)
bootstrapttest(randomerr, shuffleerr)

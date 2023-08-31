#=
Try to get KFDA working with hyperparameter tuning.
=#
using PombeAgeingGenes, LinearAlgebra, MLBase

using PombeAgeingGenes: testidxs

X, Y, goterms = load(ML)

σs = [2, 4, 8]
Ks = [Kernel((x,y)->rbf2(x,y,σ), X) for σ = σs]

y = Y[5,:]

ŷs = Vector{Float64}[]
ys = Vector{Float64}[]

# outer loop: performance
for (i, outertrainidxs) = enumerate(StratifiedKfold(y, 5))
    X_, tX = traintestsplit(X, outertrainidxs)
    y_, ty = traintestsplit(y, outertrainidxs)
    outertestidxs = testidxs(y, outertrainidxs)

    best_kernel = []

    # inner loop: hyperparameter selection
    for (j, innertrainidxs) = enumerate(StratifiedKfold(y_, 5))
        _trainidxs = outertrainidxs[innertrainidxs]
        _testidxs = testidxs(y_, innertrainidxs)

        X__, tX_ = traintestsplit(X_, innertrainidxs)
        y__, ty_ = traintestsplit(y_, innertrainidxs)
        # K__ = Kernel((x,y)->rbf2(x,y,8), X__)

        prs = []

        for K in Ks
            K__ = K[_trainidxs, _trainidxs]

            M = fit(KernelFisherDiscriminant, X__, y__, K__)

            # K___ = Kernel((x,y)->rbf2(x,y,8), X__, tX_)
            K___ = K[_trainidxs, _testidxs]

            ŷ = probability(M, K___)
            push!(prs, PR(ŷ, y__).auc)
        end
        _, best = findmax(prs)
        push!(best_kernel, best)
        println(best)
    end

    K = Ks[mode(best_kernel)]
    M = fit(KernelFisherDiscriminant, X_, y_, K[outertrainidxs, outertrainidxs])
    ŷ = probability(M, K[outertrainidxs, outertestidxs])
    push!(ŷs, ŷ)
    push!(ys, ty)
    # break
end

Performance.(ŷs, ys)
PR(vcat(ŷs...), vcat(ys...))

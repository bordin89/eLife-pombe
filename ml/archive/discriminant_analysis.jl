#=
https://www.researchgate.net/publication/308015273_Linear_vs_quadratic_discriminant_analysis_classifier_a_tutorial

http://users.isr.ist.utl.pt/~wurmd/Livros/school/Bishop%20-%20Pattern%20Recognition%20And%20Machine%20Learning%20-%20Springer%20%202006.pdf
=#

using Distributed
addprocs(3)
@everywhere begin
    cd(@__DIR__)
    push!(LOAD_PATH, ENV["GIT"], ENV["POMBEAGEINGGENES"])
    using PombeAgeingGenes, PombeAgeingGenes.GeneOntology, DataFrames, Serialization, Plots, StatPlots
    plotlyjs()
end

cd(@__DIR__)
push!(LOAD_PATH, ENV["GIT"], ENV["POMBEAGEINGGENES"])
using PombeAgeingGenes, PombeAgeingGenes.GeneOntology, DataFrames, Serialization, Plots, StatPlots
plotlyjs()

## Yeast

@everywhere X, Y, goterms = load(ML)

X, Y, goterms = load(ML)

# Run FisherLinearDiscriminant CV on all GO Slim terms
# @sync @distributed for i = 1:length(goterms)
for i = 1:length(goterms)
    p = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/FisherLinearDiscriminant"
    isdir(p) || throw(Base.IOError("Directory does not exist: $p", 1))

    goterm = goterms[i]
    y = Y[i,:]

    ps, prs = fit(CrossValidation(5), FisherLinearDiscriminant, X, y)

    fig = @df DataFrame(ps) violin(:metrics, :values, ylim=(0,1), legend=false, title=goterm)

    savefig("$p/$(goterm)_metrics.pdf")

    fig = @plotpr prs
    plot!(title=goterm)
    savefig("$p/$(goterm)_PR.pdf")
    # break
end

# Run KernelFisherDiscriminant CV on all GO Slim terms
# KernelFisherDiscriminant works really well on Iris data set, but not here

using PombeAgeingGenes: rbf

rbf2(x::AbstractVector, y::AbstractVector, σ=1) = exp(-norm(x-y, 2)^2 / (2σ^2))

K = Kernel(X)

σs = [2,4,8,16]
Ks = [Kernel((x,y)->rbf2(x,y,σ), X) for σ = σs]

K = Ks[3]

K = Kernel((x,y)->rbf2(x,y,10), X)
k = K.k

if !isposdef(k)
    eig = eigen(k)
    k += I * abs(eig.values[1])
end

K = Kernel(k, rbf2)


i = 5
y = Y[i,:]

using MLBase

k = 5

for (i, trainidxs) = enumerate(StratifiedKfold(y, k))
    X_, tX = traintestsplit(X, trainidxs, dims=2)
    y_, ty = traintestsplit(y, trainidxs)

    # Hyperparameter estimation
    # TODO kernel regularization parameter
    pss = []; prss = []

    for j = 1:length(σs)
        ps, prs = fit(CrossValidation(5,1), KernelFisherDiscriminant, X_, y_, K=Ks[j])
        push!(pss, ps); push!(prss, prs)
    end

    _, best_K = findmax(getfield.(prss, :auc))

    ps, prs = fit(CrossValidation(5,1), KernelFisherDiscriminant, X, y, K=Ks[best_K])
end












for (i, goterm) = enumerate(goterms)
    i != 5 && continue

    println(i, " ", goterm)
    println("Fitting")

    y = Y[i,:]

    pss = []
    prss = []

    # TODO Do this on the training data
    for j = 1:length(σs)
        ps, prs = fit(CrossValidation(5), KernelFisherDiscriminant, X, y, K=Ks[j])
        push!(pss, ps)
        push!(prss, prs)
    end

    _, idx = findmax(getfield.(prss, :auc))

    best_sigma = σs[idx]




    println("Plotting")

    # fig = @df DataFrame(ps) violin(:metrics, :values, ylim=(0,1), legend=false, title=goterm)
    # savefig("KernelFisherDiscriminant/$(goterm)_metrics.pdf")

    fig = @plotpr prs
    plot!(title=goterm)
    savefig("KernelFisherDiscriminant/$(goterm)_PR.pdf")
    # break
end

## Iris

# @everywhere begin
#     using RDatasets, SparseArrays, MLBase, Random, PombeAgeingGenes, Statistics
#     plotlyjs()
# end
#
# iris = dataset("datasets", "iris")
# indices = collect(1:size(iris, 1))
# shuffle!(indices)
#
# Y = string.(iris[5])[indices]
# Y = labelencode(labelmap(Y), Y)
# Y = Array{Int}(sparse(1:length(Y), Y, true))'
#
# X = convert(Matrix, iris[1:4])[indices,:]'
# X = PombeAgeingGenes.standardize(X; dims=2)
#
# @everywhere X, Y = $X, $Y
#
# @distributed for i = 1:3
#     y = Y[i,:]
#     ps, prs = fit(CrossValidation(5,10), FisherLinearDiscriminant, X, y)
#     push!(z, ps)
#
#     fig = @plotpr prs
#     # plot!(title=goterm)
#     savefig("kill_$i.pdf")
# end
#
# ps, prs = fit(CrossValidation(5,10), FisherLinearDiscriminant, X, y)
# mean(ps)
# fig = @plotpr prs
# savefig("iris_fda.pdf")
#
# K = Kernel(X)
#
# ps, prs = fit(CrossValidation(5,10), KernelFisherDiscriminant, X, y, K=K)
# mean(ps)
# fig = @plotpr prs
# savefig("iris_kfda.pdf")
#
# M = fit(KernelFisherDiscriminant, X, y)
# evaluate(M, X)
# probability(M, X)
#
# using StatPlots
# d = DataFrame(ps)
# @df d violin(:metrics, :values, ylim=(0,1))

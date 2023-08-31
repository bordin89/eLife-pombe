#=
https://arxiv.org/pdf/1804.03515.pdf
=#
using PombeAgeingGenes, DecisionTree, Distributed, JSON, Plots

addprocs(10, exeflags="--project=.")
nworkers()
@everywhere using DecisionTree

dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/RandomForestClassifier/go_slim/gp_ne"
isdir(dir) || mkpath(dir)

X, Y, goterms = load(ML, Matrix, networkembeddings=true)
# y = [i == 1 for i = Y[5,:]]

function rfc(Xtrain, ytrain, Xtest; kwargs...)
    model = DecisionTree.fit!(RandomForestClassifier(; n_trees=1000, kwargs...), Xtrain, ytrain)
    ŷ = DecisionTree.predict_proba(model, Xtest)[:,2]
end

grid = Dict(
    :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :partial_sampling => [.5, .7, 1.],
    )

function f(X, y, grid)
    ŷs, ys, p, pr = @crossvalidate X y begin
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        params = @gridsearch rfc grid
        ŷ = rfc(Xtrain, ytrain, Xtest; params...)
    end
end

for i = 1:size(Y,1)
    goterm = goterms[i]
    @show goterm
    fname = "$dir/$(goterm)"
    y = [i == 1 for i = Y[i,:]]
    ŷs, ys, p, pr = f(X, y, grid)
    savefig(plot(pr, legend=:topright), fname * ".pdf")
    writecvresults(fname * ".json", ŷs, ys)
    # break
end

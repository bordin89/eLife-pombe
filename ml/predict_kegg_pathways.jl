#=
https://arxiv.org/pdf/1804.03515.pdf
=#
using PombeAgeingGenes, DecisionTree, Distributed, JSON, Plots

addprocs(10, exeflags="--project=.")
nworkers()
@everywhere using DecisionTree

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "kegg_pathways",
               "RandomForestClassifier", "c0f6b9b", "ne")
isdir(dir) || mkpath(dir)

X, Y, pathway = load(ML, Matrix, Y=:kegg, networkembeddings=true, growthphenotypes=false)

function rfc(Xtrain, ytrain, Xtest; kwargs...)
    model = DecisionTree.fit!(RandomForestClassifier(; n_trees=100, kwargs...), Xtrain, ytrain)
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
    pathway = KEGGPATHWAYS[i]
    @show pathway
    fname = joinpath(dir, pathway)
    y = [i == 1 for i = Y[i,:]]
    count(y) < 5 && continue
    ŷs, ys, p, pr = f(X, y, grid)
    savefig(plot(pr, legend=:topright), fname * ".pdf")
    writecvresults(fname * ".json", ŷs, ys)
    break
end

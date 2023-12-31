using PombeAgeingGenes, DecisionTree, JSON, Plots, Random, DataFrames, CSV, CodecZlib

# using Distributed
# if haskey(ENV, "N_PROCS")
#     addprocs(parse(Int, ENV["N_PROCS"]), exeflags="--project=$(ENV["GIT"])/pombage")
#     println("nprocs() = ", nprocs())
# end
# @everywhere using DecisionTree

ENV["GKSwstype"] = "nul"

const N_TREES = haskey(ENV, "N_TREES") ? parse(Int, ENV["N_TREES"]) : 500

function setupdir(features)
    if haskey(ENV, "HOSTNAME") && occursin("myriad", ENV["HOSTNAME"])
        dir = features
    else
        dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "go_slim",
                       "RandomForestClassifier", features)
    end
    isdir(dir) || mkpath(dir)
    dir
end

function rfc(Xtrain, ytrain, Xtest; kwargs...)
    model = DecisionTree.fit!(RandomForestClassifier(; n_trees=N_TREES, kwargs...),
                              Xtrain, ytrain)
    ŷ = DecisionTree.predict_proba(model, Xtest)[:,2]
end

function cv(X, y, grid)
    ŷs, ys, p, pr = @crossvalidate X y begin
        Xtrain = permutedims(Xtrain)
        Xtest = permutedims(Xtest)
        params = @gridsearch rfc grid
        ŷ = rfc(Xtrain, ytrain, Xtest; params...)
    end
end

function cvgoterms(X, Y, goterms; dir, runnumber=0)
    for i = 1:size(Y,1)
        try
            goterm = string(goterms[i])
            println("goterm = ", goterm)
            y = [j == 1 for j = Y[i,:]]
            Random.seed!(runnumber+i)
            cvgoterm(X, y, goterm; dir=dir)
        catch ArgumentError
            continue
        end
    end
end

function cvgoterm(X, y, goterm; dir)
    fname = joinpath(dir, goterm)
    ŷs, ys, p, pr = cv(X, y, grid)
    savefig(plot(pr, legend=:topright), fname*".pdf")
    writecvresults(fname*".json", ŷs, ys)
end

function repeats(X, Y, goterms, nrepeats::Integer; dir::AbstractString)
    repeats(X, Y, goterms, 1, nrepeats; dir=dir)
end

function repeats(X, Y, goterms, start::Integer, stop::Integer; dir::AbstractString)
    for i = start:stop
        dir_i = joinpath(dir, string(i))
        isdir(dir_i) || mkpath(dir_i)
        cvgoterms(X, Y, goterms, dir=dir_i; runnumber=i)
    end
end

function load_john_townsend_matrix(fp, T=DataFrame)
    X = open(fp) do io
    	CSV.read(GzipDecompressorStream(io)) |> DataFrame
    end
    rename!(X, :Column1=>:id)
    dropmissing!(X)
    for col in names(X[:,Not(:id)])
    	if eltype(X[:,col]) != Float64
    		X[!,col] = parse.(Float64, X[:,col])
    	end
    end

    Y = load(GeneOntology.GOSlimTargets)
    commonids = sort(X[:, :id] ∩ Y[:, :id])
    X = sort!(X[map(id->id ∈ commonids, X[:, :id]), :], :id)
    Y = sort!(Y[map(id->id ∈ commonids, Y[:, :id]), :], :id)

    if T <: AbstractMatrix
        return permutedims.(convert.(T{Float64}, (X[:,Not(:id)], Y[:,Not(:id)])))..., names(Y)[2:end]
    end

    return X, Y
end

#=
Predict GO Slim terms using:
    - John Townsend's matrices
=#

include("src.jl")

# data_dir = "$(ENV["POMBEAGEINGGENES"])/data/john_townsend_matrices"
base_dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/john_townsend_matrices/results"

filename = ARGS[1]

dir = joinpath(base_dir, split(basename(filename), '.')[1])
isdir(dir) || mkpath(dir)

X, Y, goterms = load_john_townsend_matrix(filename, Matrix)

const grid = Dict(
    # :n_subfeatures => [floor(sqrt(size(X,1))), 10, 25, 50],
    :n_subfeatures => [10, 25, 50],
    # :partial_sampling => [.5, .75, 1.],
    )

# cvgoterms(X, Y, goterms; dir=dir)
repeats(X, Y, goterms, 5; dir=dir)

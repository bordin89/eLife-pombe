#=
Predict GO Slim terms using:
  - John Townsend's matrices

    export JULIA_NUM_THREADS=8
    export N_TREES=100
    bsub16 -n $JULIA_NUM_THREADS -J "cor_PCA_ML" "julia --project=$GIT/pombage $POMBEAGEINGGENES/Scripts/ml/john_townsend_matrices/predict_go_slim_terms_pca.jl"
=#

include("src.jl")

# data_dir = "$(ENV["POMBEAGEINGGENES"])/data/john_townsend_matrices"
base_dir = "$(ENV["POMBEAGEINGGENES"])/Scripts/ml/john_townsend_matrices/results/pca"

# filename = ARGS[1]
filename = joinpath(ENV["POMBEAGEINGGENES"], "data/john_townsend_matrices/xCor.csv.gz")

dir = joinpath(base_dir, split(basename(filename), '.')[1])
isdir(dir) || mkpath(dir)

X, Y, goterms = load_john_townsend_matrix(filename, Matrix)

using MultivariateStats

m = fit(PCA, X)
X = transform(m, X)

# clusteredheatmap(X, c=:RdBu)

const grid = Dict(
    :n_subfeatures => [Int(floor(sqrt(size(X,1)))), 5, 10],
    :partial_sampling => [.5, .75, 1.],
    )

# N_TREES = 10

# cvgoterms(X, Y, goterms; dir=dir)
repeats(X, Y, goterms, 5; dir=dir)

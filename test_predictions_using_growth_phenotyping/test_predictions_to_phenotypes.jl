#=
Take a predicted GO term
Get other close GO terms with high Resnik semantic similarity
Get all genes associated with that set of terms
Pull out phenotypes of those genes
Calculate mean pairwise Pearson correlation coefficient
Compare to bootstrap null distribution
=#

using PombeAgeingGenes
using CSV
using DataFrames
using DataFramesMeta
using Statistics
using Combinatorics
using PombeAgeingGenes
using LinearAlgebra
using Random
using UnicodePlots
using MultipleTesting

# Predictions

predictions = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "predictions", "predictions_combined.csv")
predictions = CSV.File(predictions) |> DataFrame
predictions = @where(predictions, :probability .≥ 0.7)


# Growth phenotypes

phenotypes = load(GrowthPhenotypesWideform)

PombeAgeingGenes.trigitise!(phenotypes)

M = permutedims(Matrix(phenotypes[:, Not(:id)]))


# Test

function calculate_mean_of_cor(M)
    return calculate_mean_upper_triangle(cor(M))
end

function calculate_mean_upper_triangle(M)
    t = zero(eltype(M))
    n = 0
    for j in 1:size(M, 2), i in 1:j
        x = M[i,j]
        isnan(x) && continue
        t += x
        n += 1
    end
    return t / n
end

function bootstrap(M, n_with_go; n_boot = 1000)
    rand_inds = Vector{Int}(undef, n_with_go)
    n = size(M, 2)
    map(1:n_boot) do _
        rand!(rand_inds, 1:n)
        M_sub = @views M[:, rand_inds]
        calculate_mean_of_cor(M_sub)
    end
end

function calculate_p_value(x, bs)
    return (count(bs .> x) + 1) / (length(bs) + 1)
end

function remove_high_variance_conditions(M; pc = 0.2)
    conditions_vars = vec(var(M, dims = 2))
    idxs = sortperm(conditions_vars, rev = true)
    to_remove = length(idxs) * pc
    idxs_to_keep = idxs .> to_remove
    return M[idxs_to_keep, :]
end

function calculate_cor_for_go_term(predictions, M, gene_ids, go_id)
    gene_ids_with_go_term = predictions.gene_id[predictions.go_id .== go_id]
    idxs = [gene_id in gene_ids_with_go_term for gene_id in gene_ids]
    count(idxs) == 0 && return (missing, missing)
    MM = @views M[:, idxs]
    MM = @views remove_high_variance_conditions(MM)
    x = calculate_mean_of_cor(MM)
    bs = bootstrap(M, size(MM, 2))
    p = calculate_p_value(x, bs)
    return x, p
end

gene_ids = phenotypes.id
go_ids = unique(predictions.go_id)

# go_id = go_ids[1]
# calculate_cor_for_go_term(predictions, M, gene_ids, go_id)

cache = Dict{String,Tuple{Float64,Float64}}()

n = 100
# n = length(go_ids)

for go_id in go_ids[1:n]
    if !haskey(cache, go_id)
        res = calculate_cor_for_go_term(predictions, M, gene_ids, go_id)
        ismissing(res[1]) && ismissing(res[2]) && continue
        cache[go_id] = res
    end
end

res = collect(values(cache))
xs = getfield.(res, 1)
ps = getfield.(res, 2)
qs = adjust(ps, BenjaminiHochberg())

count(isnan, xs)
count(ismissing, xs)

boxplot(xs)
boxplot(ps)
boxplot(qs)

count(qs .< 0.05)
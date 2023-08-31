cd(@__DIR__)

using PombeAgeingGenes, PombeAgeingGenes.GeneOntology, MultivariateStats, DataFrames,
    Statistics, DataFramesMeta, HypothesisTests, MLDataUtils
using Plots

function idswithgoterm(goterm::Symbol)  # in the form `:GO0006915`
    df = load(GOSlimTargets, goterm)
    Set(string.(df[df[goterm] .== 1, :][:id]))
end

function plotdecomposition(A, pc1::Int)
    scatter(A[pc1,:], A[pc1+1,:], α=.25, c=:grey, markerstrokewidth=0)
    plot!(legend=false)
end
function plotdecomposition(A, pc1::Int, geneindexes::AbstractVector)
    An = A[:,.!geneindexes]
    Ap = A[:,geneindexes]
    scatter(An[pc1,:], An[pc1+1,:], α=.25, c=:grey, markerstrokewidth=0, label="N")
    scatter!(Ap[pc1,:], Ap[pc1+1,:], c=:orange, markerstrokewidth=0, label="P")
end

function plotpca(A, pc1=1, geneindexes=nothing; kwargs...)
    geneindexes === nothing ?
        plotdecomposition(A, pc1) :
        plotdecomposition(A, pc1, geneindexes)
    plot!(xlabel="PC $pc1", ylabel="PC $(pc1+1)"; kwargs...)
end

function plottsne(A, geneindexes=nothing; kwargs...)
    geneindexes === nothing ?
        plotdecomposition(A') :
        plotdecomposition(A', geneindexes)
    plot!(axis=:off, grid=false; kwargs...)
end

ageinggoslimterms = [:GO0006915, :GO0006914, :GO0005975, :GO0006325, :GO0006281, :GO0032200]

"""
    select(A[; xlo[, xhi[, ylo[, yhi[, returninds]]]]])

Select ``xy`` points in `A` within the box formed by `xlo`, `xhi`, `ylo` and `yhi`.
"""
function select(A; indices=false,
                xlo=minimum(A[1,:]), xhi=maximum(A[1,:]),
                ylo=minimum(A[2,:]), yhi=maximum(A[2,:]))
    inds = (xlo .≤ A[1,:] .≤ xhi) .& (ylo .≤ A[2,:] .≤ yhi)
    indices ? inds : A[:, inds]
end

function idsincluster(selection, fname, df=df)
    ids = df[:id][selection]
    writedlm(fname, ids)
end

#=
http://www.scholarpedia.org/article/Mutual_information
=#
using StatsBase, Discretizers

function probability(xs)
    r = round.(xs, digits=2)
    f = values(counter(r))
    p = f ./ length(xs)
end

"""
    nbins(xs)

Return the number of bins with which to discretize `xs`.
"""
function nbins(xs::AbstractVector)::Int
    n = length(xs)
    w = 2 * iqr(xs) / cbrt(n)
    lo, hi = extrema(xs)
    fd = (hi - lo) / w
    sturges = log(2, n) + 1
    bins = isinf(fd) ? sturges : max(fd, sturges)
    ceil(Int, bins)
end

"""
    discretizer(xs...)

Return a discretizer for `xs`.
"""
function discretizer(xs::AbstractVector{T}...) where T<:Real
    n = maximum(nbins.(xs))
    d = DiscretizeUniformWidth(n)
    LinearDiscretizer(binedges(d, vcat(xs...)))
end

"""
    DiscretizedVector(xs)
    DiscretizedVector(nbins, xs)
    DiscretizedVector(binedges, xs)
    DiscretizedVector(discretizer, xs)

Discretize `xs` into a DiscretizedVector. Binning can be controlled by supplying `nbins`,
`binedges` or `discretizer`.
"""
struct DiscretizedVector{T<:Real,D<:Integer} <: AbstractVector{D}
    d::LinearDiscretizer
    xs::Vector{T}
end
function DiscretizedVector(d::LinearDiscretizer, xs::AbstractVector{<:Real})
    DiscretizedVector{eltype(xs),Int}(d, xs)
end
function DiscretizedVector(binedges::AbstractVector, xs::AbstractVector{<:Real})
    DiscretizedVector(LinearDiscretizer(binedges), xs)
end
function DiscretizedVector(nbins::Integer, xs::AbstractVector{<:Real})
    DiscretizedVector(binedges(DiscretizeUniformWidth(nbins), xs), xs)
end
function DiscretizedVector(xs::AbstractVector{<:Real})
    DiscretizedVector(nbins(xs), xs)
end

Base.size(v::DiscretizedVector) = size(v.xs)
Base.getindex(v::DiscretizedVector, i::Integer) = encode(v.d, v.xs[i])
Base.IndexStyle(v::DiscretizedVector) = IndexLinear()

"""
    discretize(xs...)
    discretize(d, xs...)

Discretize Vectors `xs`, optionally supplying the discretizer `d`.
"""
function discretize(d::LinearDiscretizer, xs::AbstractVector{<:Real}...)
    map(x->DiscretizedVector(d, x), xs)
end
function discretize(xs::AbstractVector{<:Real}...)
    discretize(discretizer(xs...), xs...)
end

"""
    frequencies(v)

Bin frequencies of a DiscretizedVector.
"""
function frequencies(v::DiscretizedVector{T,D}) where {T,D}
    f = zeros(D, v.d.nbins)

    for i = v
        f[i] += one(D)
    end

    f
end

"""
    mi(xs, ys)

Mutual information of frequency vectors `xs` and `ys` in bits.

Ref: https://github.com/rmaestre/Mutual-Information/blob/master/it_tool.py
"""
function mi(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    n = length(xs)
    m = length(ys)
    n == m || throw(DimensionMismatch("Inconsistent array length."))

    mi = 0.

    for i = Set(xs), j = Set(ys)
        xi = xs .== i
        yj = ys .== j

        Px  = sum(xi) / n
        Py  = sum(yj) / n
        Pxy = sum(xi .& yj) / n

        P = Pxy * log2(Pxy / (Px * Py))

        if isfinite(P)
            mi += P
        end
    end

    mi
end
function mi(xs::DiscretizedVector, ys::DiscretizedVector)
    mi(frequencies(xs), frequencies(ys))
end

xs = [1.1,1.3,0.9,1.,2.,2.3,3.1]
ys = [1,1,1,2,2,2,3.,1]
dx, dy = discretize(xs, ys)
fx, fy = frequencies.((dx, dy))
mi(fx, fy)


# function frequencies(xs::AbstractVector{T}, bins::Int, vmin::T, vmax::T) where T<:Real
#     h = fit(Histogram, xs, range(vmin, vmax+vmin, length=bins+1), closed=:left)
#     fx = h.weights
# end
# function frequencies(xs::AbstractVector{T}, ys::AbstractVector{T}...) where T<:Real
#     vmin, vmax = extrema(vcat(xs, ys...))
#     b = maximum(bins.((xs, ys...)))
#     frequencies.((xs, ys...), b, vmin, vmax)
# end

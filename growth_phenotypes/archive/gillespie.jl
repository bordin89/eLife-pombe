#=
https://link.springer.com/article/10.1007%2Fs11538-017-0356-4
https://www.khanacademy.org/science/biology/ecology/population-growth-and-regulation/a/exponential-logistic-growth
=#
cd(@__DIR__)
using Distributions, DataFrames
using Plots, StatPlots

Initial = Normal(10_000, 1000)

"SSA for cell cycle with one stage."
function ssa(T::Float64)
    θ = 180
    λ = 1/θ
    K = 100_0000
    N = round(Int, rand(Initial))
    t = 0.

    data = [(t=t, N=N)]

    while t < T
        A = λ * ((K-N)*N)/K
        τ = log(1 / (1-rand())) / A
        t += τ
        N += 2
        push!(data, (t=t, N=N))
    end

    data
end

data = ssa(1000.)
plot(getindex.(data, :t), getindex.(data, :N), xlabel="t (min)", ylabel="Cells")

d = [ssa(1000.) for i = 1:100]
plot([getindex.(i, :t) for i = d], [getindex.(i, :N) for i = d], α=.5, c=:grey, legend=false, xlabel="t (min)", ylabel="Cells")
savefig("gillespie_carrying_capacity.png")

mids = [i[round(Int, length(i)/2)].t for i = d]

function getsize(x)
    for i = x
        i.t > 500 && return i.N
    end
end

mids = getsize.(d)

using Statistics
coefvar(xs) = std(xs) / mean(xs)

coefvar(mids)

"SSA for cell cycle with four stages."
function f(T::Float64)
    θ = 180
    λ = 1/θ
    λ /= 4
    N = round(Int, rand(Initial))
    N /= 4
    G2 = N
    M  = N
    G1 = N
    S  = N
    X = [G2, M, G1, S]
    k = length(X)
    t = 0.
    data = [(t=t, G2=G2, M=M, G1=G1, S=S)]

    while t < T
        a = X .* λ
        A = sum(a)

        τ = log(1 / (1-rand())) / A
        t += τ

        j = sum(cumsum(a) .< (A * rand())) + 1

        X[j] -= 1

        if j != k
            X[j+1] += 1
        end
        if j == k
            X[1] += 2
        end

        push!(data, (t=t, G2=X[1], M=X[2], G1=X[3], S=X[4]))
    end
    data
end

data = f(10000.)

df = DataFrame(t = Float64[], G2 = Int[], M = Int[], G1 = Int[], S = Int[]) #, G0 = Int[])
for i = data
    push!(df, i)
end
df[:total] = map(x->sum(values(x[2:end])), eachrow(df))

@df df plot(:t, :G2, label="G2", legend=:topleft, ylabel="Cells", xlabel="t (min)")
@df df plot!(:t, :M, label="M")
@df df plot!(:t, :G1, label="G1")
@df df plot!(:t, :S, label="S")
@df df plot!(:t, :total, label="Total")
savefig("gillespie_four_poisson_processes.pdf")

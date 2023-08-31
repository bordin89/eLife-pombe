#=
Pombe doubling times

https://www.sciencedirect.com/science/article/pii/S0960982215015092?via%3Dihub
2.34 ± 0.18 hr = 140.4 ± 10.8 min
λ = 1 / 140.4 = 0.00712 probability of doubling per min


https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5526333/
3 hr = 180 min

https://link.springer.com/article/10.1007%2Fs11538-017-0356-4
https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005843
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436811/

Methodology
sample x cells from a normal (?) dbn

sample times at which these cells will double from a uniform (?) dbn between [0,140]

loop over time from t=t₁ to t=T
    for all cells `y` whose doubling time == t
        create `y` new cells
        sample their doubling time from a normal distribution N(140, 10)
        update doubling time for all `y` cells

Count how many cells are present
=#
cd(@__DIR__)

using Distributions, Plots #, StatPlots

θ = 180
Division = Exponential(θ)
Initial = Normal(10_000, 1000)
Synchrony = Uniform(0., θ)

function pin()
    n = round(Int, rand(Initial))
    colony = round.(Int, rand(Synchrony, n))
end

function celldivision!(colony::AbstractVector, t::Int)
    # Indexes of cells dividing at time t
    dividing_cells = colony .== t
    # Set next division time for mother cells
    colony[dividing_cells] .= round.(Int, rand.(Division)) .+ t
    # Add new daughter cells and set next division time
    append!(colony, (round.(Int, rand(Division, sum(dividing_cells))) .+ t))
end

function colony(T=1000)
    # Initial colony plated by the pin
    c = pin()
    sizehint!(c, round(Int, 2^((T/θ)) * length(c)))

    # Divide cells
    size = Vector{Int}(undef, T+1)

    for t = 1:T
        size[t] = length(c)
        celldivision!(c, t)
    end

    size[T+1] = length(c)

    # Get total colony size at time T
    size
end

c = colony(1000)

cs = [colony(1000) for i = 1:100]

finals = [c[end] for c = cs]

using Statistics
coefvar(xs) = std(xs) / mean(xs)

coefvar(finals)

plot(cs, α=.5, c=:grey, legend=false, xlabel="t (min)", ylabel="Cells")
savefig("poisson_process.pdf")
plot(cs, α=.5, c=:grey, legend=false, xlabel="t (min)", ylabel="Cells", yscale=:log10)
savefig("poisson_process_logy.pdf")

#=
Singular value decomposition of input data and plot decay of singular values.
=#

using PombeAgeingGenes, LinearAlgebra, Plots

df = load(GrowthPhenotypesWideform)
A = Matrix(df[2:end])

x = load(ML, growthphenotypes=false, networkembeddings=true)
B = Matrix(x[1][2:end])

s = svd(A)
t = svd(B)

plot(s.S, label="Growth phenotypes", xlabel="ith singular value", ylabel="Singular value") #, y_scale=:log10)
plot!(t.S, label="Network embeddings")
savefig(joinpath(@__DIR__, "svd.pdf"))

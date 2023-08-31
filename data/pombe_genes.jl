#=
Parse STRING v10.5 database to obtain a sorted list of the 5100 S. pombe genes present.
=#

using DelimitedFiles

x = readdlm("Data/deepNF_embeddings/4896.protein.links.detailed.v10.5.txt")
y = Array{String}(x[2:end,1:2][:])
z = map(x->x[6:end-2], unique(y))
sort!(z)
writedlm("Data/deepNF_embeddings/pombe_genes.csv", z)

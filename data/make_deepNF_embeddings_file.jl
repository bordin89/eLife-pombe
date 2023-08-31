using PombeAgeingGenes, MAT, DelimitedFiles, CSV

dir = joinpath(ENV["POMBEAGEINGGENES"], "data", "network_embeddings")
embeddings = read(matopen(joinpath(dir, "yeast_MDA_arch_256-256_embeddings.mat")), "embeddings")
genes = string.(vec(readdlm(joinpath(dir, "pombe_genes.csv")))) # row labels

df = DataFrame(embeddings)
df.id = genes

# Select genes in Bioneer library
X, Y = load(ML) # Genes in Bioneer library
idx = [i ∈ X.id for i = genes]
df = DataFrame(embeddings[idx, :])
df.id = genes[idx]


df = df[:,sort(names(df))]
CSV.write(joinpath(dir, "network_embeddings.csv"), df)

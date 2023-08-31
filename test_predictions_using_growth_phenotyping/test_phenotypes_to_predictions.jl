#=
Cluster phenotypes by hierarchical clustering
For various clusterings by cutting the tree
Look at GO terms predicted for genes in each cluster
Is the distribution of semantic similarities closer than expected by chance
Compare to bootstrap null distribution
=#

using PombeAgeingGenes
using PombeAgeingGenes.GeneOntology
using Clustering
using Distances
using DataStructures
using CSV
using DataFrames
using DataFramesMeta
using OBOParse
using Statistics
using Combinatorics

# Predictions

predictions = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "predictions", "predictions_combined.csv")
predictions = CSV.File(predictions) |> DataFrame
predictions = @where(predictions, :probability .≥ 0.7)


# Growth phenotypes

phenotypes = load(GrowthPhenotypesWideform)
# PombeAgeingGenes.trigitise!(phenotypes)
M = Matrix(phenotypes[:, Not(:id)])'
gene_ids = phenotypes.id

# Information content

ontology = load(GO)

alt_ids = begin
    d = Dict{String,String}()
    for (go_id, go_term) in ontology.terms
        if haskey(go_term.tagvalues, "alt_id")
            for alt_id in go_term.tagvalues["alt_id"]
                d[alt_id] = go_id
            end
        end
    end
    d
end

annotations = load(GOAnnotations, ecs = EVIDENCE_CODES_ALL)

N = length(annotations)

term_counts = counter(getfield.(annotations, :go))

for (go_id, freq) in copy(term_counts)
    if haskey(alt_ids, go_id)
        alt_id = alt_ids[go_id]
        term_counts[alt_id] = freq
    end
end

function information_content(
    go_id::AbstractString,
    term_counts::AbstractDict,
    N::Integer,
    ontology::Ontology,
)
    freq = term_counts[go_id]

    if !haskey(ontology.terms, go_id)
        go_id = alt_ids[go_id]
    end

    for go_term in descendants(ontology, ontology.terms[go_id])
        go_id = string(go_term)
        freq += term_counts[go_id]
    end

    p = freq / N
    ic = -log(p)
    return ic
end

d_ic = Dict{String,Float64}()

for go_id in unique(getfield.(annotations, :go))
    if !haskey(d_ic, go_id)
        ic = information_content(go_id, term_counts, N, ontology)
        d_ic[go_id] = ic
    end

    if !haskey(ontology.terms, go_id)
        go_id = alt_ids[go_id]
        d_ic[go_id] = ic
    end

    for a_term in ancestors(ontology, ontology.terms[go_id])
        a_id = string(a_term)

        if !haskey(d_ic, a_id)
            d_ic[a_id] = information_content(a_id, term_counts, N, ontology)
        end
    end
end

max_ic = maximum(collect(values(d_ic)))

d_ic = DefaultDict(max_ic, d_ic)

function resnik_semantic_similarity(go1, go2, ontology, d_information_content)
    a1 = ancestors(ontology, go1)
    a2 = ancestors(ontology, go2)

    intersection = intersect!(a1, a2)
    length(intersection) == 0 && return missing

    max_ic = 0.0

    for go_term in intersection
        go_id = string(go_term)
        ic = d_information_content[go_id]
        if ic > max_ic
            max_ic = ic
        end
    end

    return max_ic
end

function calculate_average_resnik_semantic_similarity(go_ids, ontology, d_information_content)
    go_terms = [ontology.terms[go_id] for go_id in go_ids]

    t = 0.0
    n = 0

    for namespace in [
        "biological_process",
        "molecular_function",
    ]
        for (go1, go2) in combinations(go_terms, 2)
            go1.namespace != namespace && continue
            go2.namespace != namespace && continue

            rss = resnik_semantic_similarity(go1, go2, ontology, d_information_content)

            if !ismissing(rss)
                t += rss
                n += 1
            end
        end
    end

    return t / n
end

function calculate_average_resnik_semantic_similarity_for_clusters(clusters, gene_ids, predictions, ontology, alt_ids, d_ic)
    map(unique(clusters)) do cluster
        genes_in_cluster = Set(gene_ids[clusters .== cluster])
        go_ids_in_cluster = unique!(predictions.go_id[[gene_id in genes_in_cluster for gene_id in predictions.gene_id]])
        for i in 1:length(go_ids_in_cluster)
            go_id = go_ids_in_cluster[i]
            if haskey(alt_ids, go_id)
                go_ids_in_cluster[i] = alt_ids[go_id]
            end
        end

        rss = calculate_average_resnik_semantic_similarity(go_ids_in_cluster, ontology, d_ic)
        # TODO bootstrap
    end
end

# Cluster

C = cor(M)
D = pairwise(CorrDist(), M, dims = 2)
any(isnan.(D))

h = hclust(D, linkage = :ward)

# UnicodePlots.heatmap(
#     C[h.order, h.order],
# )

# Plots.heatmap(
#     M'[h.order, :],
#     c = cgrad(:RdBu, rev=true),
#     clims = (-1., 1.),
# )

clusters = cutree(h, k = 13)

counter(clusters)

rss = calculate_average_resnik_semantic_similarity_for_clusters(clusters, gene_ids, predictions, ontology, alt_ids, d_ic)

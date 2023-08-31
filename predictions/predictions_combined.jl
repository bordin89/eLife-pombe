#=
Combine ML and CAFT predictions
=#

using CSV
using DataFrames
using DataFramesMeta
using PombeAgeingGenes
using PombeAgeingGenes.GeneOntology
using DataStructures
using OBOParse

# Load data
caft = CSV.File(joinpath(@__DIR__, "predictions_caft.csv")) |> DataFrame
DataFrames.rename!(caft, "probability"=>"probability_caft")

ml = CSV.File(joinpath(@__DIR__, "predictions_ml.csv")) |> DataFrame
DataFrames.rename!(ml, "probability"=>"probability_ml")

# Functions predicted by both methods

df_intersection = join(
    caft,
    ml,
    on = [
        :gene_id,
        :gene_name,
        :go_id,
        :go_name,
        :go_namespace,
        :pombase_evidence_code,
        :is_priority_unstudied_gene,
        :is_conserved_in_vertebrates,
        :is_unknown_function_gene,
    ],
    makeunique = true,
    kind = :inner,
)

# Set the probability to p(ML) + p(CAFT), or 1
# df_intersection[!, :probability] = ones(size(df_intersection, 1))
df_intersection[!, :probability] = map(zip(df_intersection.probability_ml, df_intersection.probability_caft)) do (m, c)
    min(m + c, 1)
end

histogram(df_intersection.probability)

# Functions predicted by only one method

df_union = join(
    caft,
    ml,
    on = [
        :gene_id,
        :gene_name,
        :go_id,
        :go_name,
        :go_namespace,
        :pombase_evidence_code,
        :is_priority_unstudied_gene,
        :is_conserved_in_vertebrates,
        :is_unknown_function_gene,
    ],
    makeunique = true,
    kind = :outer,
)

# Set the probability to the probability of the method that predicted it

df_union[!, :probability] = map(zip(df_union.probability_ml, df_union.probability_caft)) do (m, c)
    if ismissing(m) && !ismissing(c)
        c
    elseif !ismissing(m) && ismissing(c)
        m
    else
        missing
    end
end

# Keep the predicted functions with the highest probability

df = vcat(df_intersection, df_union)
unique!(df, [:gene_id, :go_id])
sort!(df, :probability, rev = true)

for col in [:probability_ml, :probability_caft]
    df[!, col] = coalesce.(df[!, col], 0.0)
end
df[!, :pombase_evidence_code] = coalesce.(df.pombase_evidence_code, "")
disallowmissing!(df)

histogram(df.probability)

# Sanity check that the probabilities from `df_intersection` were kept
#
# join(
#     df_intersection,
#     df,
#     on = [:gene_id, :go_id, :probability],
#     kind = :inner,
#     makeunique = true,
# )

# Add GO term information content
#
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2518162/

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

using UnicodePlots

histogram(collect(values(d_ic)))

max_ic = maximum(collect(values(d_ic)))

df[!, :information_content] = map(df.go_id) do go_id
    get(d_ic, go_id, max_ic)
end

histogram(df.information_content)

using StatsPlots

violin(
    [
        [d_ic[x] for x in getfield.(annotations, :go)],
        df.information_content,
    ],
    ylabel = "Information content",
    xticks = (
        1:2,
        ["Annotations", "Predictions"],
    ),
    legend = false,
    c = :grey80,
)

savefig(joinpath(@__DIR__, "information_content_annotations_vs_predictions.pdf"))

mean([d_ic[x] for x in getfield.(annotations, :go)])
mean(df.information_content)

median([d_ic[x] for x in getfield.(annotations, :go)])
median(df.information_content)

using HypothesisTests

pvalue(
    MannWhitneyUTest(
        [d_ic[x] for x in getfield.(annotations, :go)],
        df.information_content,
    ),
    tail = :right,
)

# Write CSV

df = df[:, [
    :gene_id,
    :gene_name,
    :go_id,
    :go_name,
    :go_namespace,
    :probability,
    :probability_ml,
    :probability_caft,
    :information_content,
    :pombase_evidence_code,
    :is_priority_unstudied_gene,
    :is_conserved_in_vertebrates,
    :is_unknown_function_gene,
]]

CSV.write(
    joinpath(@__DIR__, "predictions_combined.csv"),
    df,
)

# df = CSV.File(
#     joinpath(@__DIR__, "predictions_combined.csv"),
# ) |> DataFrame
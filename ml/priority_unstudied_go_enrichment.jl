using PombeAgeingGenes
using JSON
using .GeneOntology
using OBOParse
using DataFrames
using DataFramesMeta
using HypothesisTests
using MultipleTesing


function loaddata(dir::AbstractString, goterms::AbstractVector{<:AbstractString}=String[])
    df = DataFrame(id = String[], go = String[], y = Bool[], score = Float64[])

    if isempty(goterms)
        goterms = unique(map(x->split(x, ".")[1], readdir(dir)))
    end

    for goterm = goterms
        ŷs, ys, yids = try
            loaddata(dir, goterm)
        catch
            continue
        end

        goterm = replace(goterm, "GO"=>"GO:")

        for (ŷ, y, id) in zip(ŷs, ys, yids)
            push!(df, (id = id, go = goterm, y = y, score = ŷ))
        end
    end

    return df
end

function loaddata(dir::AbstractString, goterm::AbstractString)
    d = JSON.parsefile(joinpath(dir, goterm*".json"))
    ŷs = Array{Array{Float64}}(d["ŷs"])
    ys = Array{Array{Bool}}(d["ys"])
    yids = Array{Array{String}}(d["yids"])
    vcat(ŷs...), vcat(ys...), vcat(yids...)
end

# Known annotations
function get_annotations_and_ancestors(ontology, ecs=GeneOntology.EVIDENCE_CODES[:experimental])
    annotations = load(GOAnnotations, ecs=ecs)
    s = Set{Tuple{String,String}}()
    for x in annotations
        push!(s, (x.id, x.go))
        if haskey(ontology.terms, x.go)
            for a in ancestors(ontology, ontology[x.go])
                push!(s, (x.id, string(a)))
            end
        end
    end
    return s
end

go = load(GO)

dir = joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "go_slim", "RandomForestClassifier", "prediction", "b3dfb89", "ne_ff")

df = loaddata(dir)

# Remove known annotations with experimental or curated evidence codes

annotations_trusted = get_annotations_and_ancestors(
    go,
    EVIDENCE_CODES_TRUSTED,
)

filter!(r->(r.id, r.go) ∉ annotations_trusted, df)

# Remove low probability predictions
filter!(r->r.score > 0.1, df)

# Priority unstudied

fp = joinpath(ENV["POMBEAGEINGGENES"], "data", "priority_unstudied_genes.tsv")

priority_unstudied_ids = Set(readlines(fp)[2:end])

df = @in(df, :id, priority_unstudied_ids)

# Functional enrichment

function proteins_annotated_with_gos(annotations)
    gos = unique(getfield.(annotations, :go))
    d = Dict{String, Int}()
    for go in gos
        has_go = count(x->x.go == go, annotations)
        d[go] = has_go
    end
    return d
end

function go_enrichment(has_prediction, n_genes_predicted, has_annotation, n_genes_annotated)
    pvalue(
        PowerDivergenceTest(
            reshape(
                [
                    has_prediction, n_genes_predicted - has_prediction,
                    has_annotation, n_genes_annotated - has_annotation,
                ],
                (2, 2)
            )
        ), tail=:right
    )
end

annotations = load(GOAnnotations, ecs=EVIDENCE_CODES[:experimental])

go_2_n_proteins = proteins_annotated_with_gos(annotations)

n_genes_predicted = length(unique(tmp.id))
n_genes_annotated = length(unique(getfield.(annotations, :id)))

df_enrichment = DataFrame(
    go = String[],
    p = Float64[],
    has_prediction = Int[],
    has_annotation = Int[],
)

for g in groupby(tmp, :go)
    go = g.go[1]

    haskey(go_2_n_proteins, go) || continue
    has_annotation = go_2_n_proteins[go]

    has_prediction = size(g, 1)

    p = go_enrichment(
        has_prediction,
        n_genes_predicted,
        has_annotation,
        n_genes_annotated,
    )

    push!(
        df_enrichment,
        (
            go,
            p,
            has_prediction,
            has_annotation,
        )
    )
end

sort!(df_enrichment, :p)

df_enrichment.q = adjust(df_enrichment.p, Bonferroni())

@where(df_enrichment, :q .< 0.01)

df_enrichment

histogram(df_enrichment.q)

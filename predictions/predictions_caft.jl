#=
Convert FunFam-based CAFA 4 predictions to PomBase IDs
=#

using PombeAgeingGenes
using PombeAgeingGenes.GeneOntology
using CSV
using DataFrames
using OBOParse
using DataFramesMeta
using CodecZlib
using DataStructures

ontology = load(GO)

# Alternate GO IDs

alt_ids = begin
    d = Dict{String,String}()
    for (s,t) in ontology.terms
        if haskey(t.tagvalues, "alt_id")
            for alt in t.tagvalues["alt_id"]
                d[alt] = s
            end
        end
    end
    d
end

# Taxon constraints

taxa = [
    "Eukaryota",
    "Fungi or Dictyostelium",
    "Fungi",
    "Ascomycota",
    "Schizosaccharomycetes",
    "Schizosaccharomycetales",
    "Schizosaccharomycetaceae",
    "Schizosaccharomyces",
    "Schizosaccharomyces pombe",
]

function load_taxon_constraints(path, taxa; not_in=false)
    df = CSV.File(path, delim=',') |> DataFrame
    if not_in
        taxa = setdiff(Set(df[:,:taxon_label]), taxa)
    end
    df = @in(df, :taxon_label, taxa)
    s = Set(string.(rstrip.(df[:,:defined_class])))
    for gt in collect(s)
        if haskey(alt_ids, gt)
            gt = alt_ids[gt]
            push!(s, gt)
        end
        for a in descendants(ontology, ontology[gt], [:is_a, :regulates, :part_of])
            push!(s, string(a))
        end
    end
    return s
end

never_in_taxon = load_taxon_constraints(
    joinpath(ENV["POMBEAGEINGGENES"], "data", "never_in_taxon.tsv"),
    taxa,
)

# CAFT predictions

# Model 1: FunFam
# caft_predictions_dir = joinpath(ENV["POMBEAGEINGGENES"], "data", "caft", "TEAM1_MODEL1-0.1.1")
# caft_predictions_files = ["1.BP.pred", "1.MF.pred"]

# Model 2: FunFam + Pfam FunFam
caft_predictions_dir = joinpath(ENV["POMBEAGEINGGENES"], "data", "caft", "TEAM1_MODEL2-0.1.1")
caft_predictions_files = ["2.BP.pred", "2.MF.pred"]

struct PredictedFunction
    gene_id::String
    go_id::String
    probability::Float64
end

function read_caft_predictions(f)
    xs = PredictedFunction[]
    for line in eachline(f)
        startswith(line, "END") && break
        (startswith(line, "AUTHOR") || startswith(line, "MODEL")) && continue
        l = split(line)
        fp = PredictedFunction(l[1], l[2], parse(Float64, l[3]))
        push!(xs, fp)
    end
    return xs
end

caft_predictions = [
    read_caft_predictions(joinpath(caft_predictions_dir, caft_predictions_files[1]));
    read_caft_predictions(joinpath(caft_predictions_dir, caft_predictions_files[2]))
]

# Remove 'never in taxon'

filter!(x->x.go_id ∉ never_in_taxon, caft_predictions)

# Map CAFA IDs to PomBase IDs

function cafa_to_uniprot_map()
    d = Dict{String,String}()
    for line in eachline(joinpath(ENV["POMBEAGEINGGENES"], "data/cafa4/mapping.284812.map"))
        l = split(line)
        d[l[1]] = l[2]
    end
    return d
end

function uniprot_to_pombase_map()
    d = Dict{String,String}()
    for line in eachline(joinpath(ENV["POMBEAGEINGGENES"], "data/cafa4/uniprot_to_pombase.map"))
        l = split(line)
        d[l[1]] = l[2]
    end
    return d
end

function map_identifiers(predictions)
    d1 = cafa_to_uniprot_map()
    d2 = uniprot_to_pombase_map()
    mapped = typeof(predictions)()
    for p in predictions
        haskey(d1, p.gene_id) || continue
        uniprot_id = d1[p.gene_id]
        haskey(d2, uniprot_id) || continue
        pombase_id = d2[uniprot_id]
        push!(mapped, PredictedFunction(pombase_id, p.go_id, p.probability))
    end
    return mapped
end

caft_predictions = map_identifiers(caft_predictions)

# Filter out known pombe annotations

# annotations = begin
#     xs = load(
#         GOAnnotations,
#         ecs = EVIDENCE_CODES[:experimental], # TODO decide which evidence codes to trust
#     )
    
#     s = Set(map(x->(x.id, x.go), xs))

#     for (id, gt) in collect(s)
#         if haskey(alt_ids, gt)
#             gt = alt_ids[gt]
#             push!(s, (id, gt))
#         end

#         for a in ancestors(ontology, ontology[gt], [:is_a, :regulates, :part_of])
#             push!(s, (id, string(a)))
#         end
#     end
#     s
# end

# filter!(x->(x.gene_id, x.go_id) ∉ annotations, caft_predictions)

# Remove root terms of ontologies
filter!(x->x.go_id != "GO:0008150", caft_predictions) # remove `biological_process` term
filter!(x->x.go_id != "GO:0003674", caft_predictions) # remove `molecular_function` term
filter!(x->x.go_id != "GO:0005575", caft_predictions) # remove `cellular_component` term

# DataFrame

df = DataFrame(caft_predictions)
sort!(df, :probability, rev=true)

# Add evidence code

ec_exp = Set(EVIDENCE_CODES[:experimental])
ec_curated = Set(EVIDENCE_CODES[:curated])
ec_automatic = Set(EVIDENCE_CODES[:automatic])

# xs = load(
#         GOAnnotations,
#         ecs = EVIDENCE_CODES_ALL,
# )

# filter!(x->x.id == "SPCC320.03", xs)

annotations = begin
    xs = load(
        GOAnnotations,
        ecs = EVIDENCE_CODES_ALL,
    )

    d = Dict{Tuple{String,String}, Int}()

    for x in xs
        id = x.id
        gt = x.go
        ec = x.ec

        best_ec = if ec in ec_exp
            1
        elseif ec in ec_curated
            2
        elseif ec in ec_automatic
            3
        else
            4
        end

        if haskey(alt_ids, gt)
            gt = alt_ids[gt]
        end

        k = (id, gt)

        if haskey(d, k)
            curr_ec = d[k]
            
            if best_ec < curr_ec
                d[k] = best_ec
            end
        else
            d[k] = best_ec
        end

        for a in ancestors(ontology, ontology[gt], [:is_a, :regulates, :part_of])
            k = (id, string(a))

            if haskey(d, k)
                curr_ec = d[k]
                
                if best_ec < curr_ec
                    d[k] = best_ec
                end
            else
                d[k] = best_ec
            end
        end
    end
    d
end

df[!, :pombase_evidence_code] = map(zip(df.gene_id, df.go_id)) do id_gt
    if haskey(annotations, id_gt)
        v = annotations[id_gt]
        
        if v == 1
            "experimental"
        elseif v == 2
            "curated"
        elseif v == 3
            "automatic"
        else
            missing
        end
    else
        missing
    end
end




# annotations = begin
#     xs = load(
#         GOAnnotations,
#         ecs = EVIDENCE_CODES_ALL,
#     )

#     d = Dict{Tuple{String,String}, String}()

#     for x in xs
#         id = x.id
#         gt = x.go
#         ec = x.ec

#         best_ec = if ec in ec_exp
#             "experimental"
#         elseif ec in ec_curated
#             "curated"
#         elseif ec in ec_automatic
#             "automatic"
#         else
#             ""
#         end

#         d[(id, gt)] = best_ec


#         if haskey(alt_ids, gt)
#             gt = alt_ids[gt]
#             d[(id, gt)] = best_ec
#         end

#         for a in ancestors(ontology, ontology[gt], [:is_a, :regulates, :part_of])
#             d[(id, string(a))] = best_ec
#         end
#     end
#     d
# end

# df[!, :pombase_evidence_code] = map(zip(df.gene_id, df.go_id)) do id_gt
#     if haskey(annotations, id_gt)
#         annotations[id_gt]
#     else
#         missing
#     end
# end

# Normalise probabilities
df.probability .-= minimum(df.probability)
df.probability ./= maximum(df.probability)

# Remove low probability predictions

probability_threshold = 0.1

df = @where(df, :probability .≥ probability_threshold)

# Add gene name

function id_to_name_mapper()
    d = Dict{String,String}()
    open(GzipDecompressorStream, joinpath(ENV["POMBEAGEINGGENES"], "data", "funfam",
                                          "v4_2_0", "pombase_pombe_genome", "peptide.fa.gz")
    ) do io
        for line in eachline(io)
            if startswith(line, ">")
                l = split(line, ":pep")
                id = l[1][2:end]
                name = if occursin("|", l[2])
                    split(l[2], "|")[2]
                else
                    l[2]
                end
                d[id] = name
            end
        end
    end
    return d
end

id_to_name = id_to_name_mapper()

df[!, :gene_name] = [id_to_name[id] for id in df.gene_id]

# Add GO description

df[!, :go_name] = map(df.go_id) do go
    if haskey(alt_ids, go)
        go = alt_ids[go]
    end
    ontology[go].name
end

df[!, :go_namespace] = map(df.go_id) do go
    if haskey(alt_ids, go)
        go = alt_ids[go]
    end

    n = ontology.terms[go].namespace

    return if n == "biological_process"
        "BP"
    elseif n == "molecular_function"
        "MF"
    elseif n == "cellular_component"
        "CC"
    end
end


# Priority unstudied genes

priority_unstudied_genes = Set(
    CSV.read(
        joinpath(ENV["POMBEAGEINGGENES"], "data", "priority_unstudied_genes.tsv")
    )[:,1]
)

df[!, :is_priority_unstudied_gene] = [id in priority_unstudied_genes for id in df.gene_id]

# Vertebrate genes

vertebrate_genes = Set(
    CSV.read(
        joinpath(ENV["POMBEAGEINGGENES"], "data", "pombe_genes_conserved_in_vertebrates.tsv")
    )[:,1]
)

df[!, :is_conserved_in_vertebrates] = [id in vertebrate_genes for id in df.gene_id]

# Unknown function genes

unknown_function_genes = Set(
    CSV.read(
        joinpath(ENV["POMBEAGEINGGENES"], "data", "pombe_genes_unknown_function.tsv")
    )[:,1]
)

df[!, :is_unknown_function_gene] = [id in unknown_function_genes for id in df.gene_id]

df = df[:, [
    :gene_id,
    :gene_name,
    :go_id,
    :go_name,
    :go_namespace,
    :probability,
    :pombase_evidence_code,
    :is_priority_unstudied_gene,
    :is_conserved_in_vertebrates,
    :is_unknown_function_gene,
]]


# CSV

CSV.write(
    joinpath(@__DIR__, "predictions_caft.csv"),
    df,
)

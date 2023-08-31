#=
Predict GO terms for all pombe genes.
=#
using Revise
using PombeAgeingGenes, CSV, DataStructures, DataFrames, OBOParse, DataFramesMeta

ontology = load(GeneOntology.GO)

alt_ids = (function alt_ids_f()
    d = Dict{String,String}()
    for (s,t) in ontology.terms
        if haskey(t.tagvalues, "alt_id")
            for alt in t.tagvalues["alt_id"]
                d[alt] = s
            end
        end
    end
    return d
end)()

# taxon constraints
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
    df = CSV.read(path)
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
    joinpath(ENV["POMBEAGEINGGENES"], "Data", "never_in_taxon.tsv"),
    taxa,
)

only_in_taxon = load_taxon_constraints(
    joinpath(ENV["POMBEAGEINGGENES"], "Data", "only_in_taxon.tsv"),
    taxa,
    not_in=true,
)

#=
s = Set(map(x->x.go, load(GeneOntology.GOAnnotations)))
s ∩ never_in_taxon
s ∩ only_in_taxon
=#

annotations = begin
    s = Set(map(x->(x.id, x.go), load(GeneOntology.GOAnnotations)))
    for (id, gt) in collect(s)
        gt = get(alt_ids, gt, gt)
        if gt in never_in_taxon || gt in only_in_taxon
            continue
        end
        if haskey(alt_ids, gt)
            gt = alt_ids[gt]
            push!(s, (id, gt))
        end
        for a in ancestors(ontology, ontology[gt], [:is_a, :regulates, :part_of])
            push!(s, (id, string(a)))
        end
    end
    s
end

# Priority unstudied genes
priority_unstudied_genes = Set(CSV.read(joinpath(ENV["POMBEAGEINGGENES"], "Data", "priority_unstudied_genes.tsv"))[:,1])

# Vertebrate genes
vertebrate_genes = Set(CSV.read(joinpath(ENV["POMBEAGEINGGENES"], "Data", "pombe_genes_conserved_in_vertebrates.tsv"))[:,1])

# funfam inclusion thresholds
tcs = load(FunFamInclusionThresholds)
tcs_dict = Dict(map(x->(x.id, x.tc), tcs))

# pombe funfam hits
hits = load(FunFamHits)

# only keep sequences that meet the funfam inclusion threshold
# filter!(x->x.score ≥ tcs_dict[x.ff], hits)

# GO terms associated with FFs
funfam_goterms = load(FunFamGOTerms)
# funfam_goterms = load(FunFamGOTermsIEA)

funfam_goterms_dict = (function funfam_goterms_dict_f()
    d = DefaultDict{String,Vector{String}}([])
    for (ff,gt) in funfam_goterms
        if gt in never_in_taxon || gt in only_in_taxon
            continue
        end
        push!(d[ff], gt)
    end
    return d
end)()

function funfam_goterm_predictions(hits,
                                   funfam_goterms_dict,
                                   ontology,
                                   annotations,
                                   tcs_dict,
                                   priority_unstudied_genes,
                                   vertebrate_genes;
                                   include_uniprotkb_kw_iea,
                                   )
    predictions = []
    for hit in hits
        goterms = funfam_goterms_dict[hit.ff]
        for goterm in goterms
            goterm = get(alt_ids, goterm, goterm) # TEMP
            funfam_inclusion_threshold = tcs_dict[hit.ff]
            is_in_funfam = hit.score ≥ funfam_inclusion_threshold
            push!(predictions, (
                identifier = hit.id,
                funfam = hit.ff,
                goterm = goterm,
                goterm_name = ontology[goterm].name,
                ontology = ontology[goterm].namespace,
                is_in_pombe_annotations = (hit.id, goterm) ∈ annotations,
                is_priority_unstudied_gene = hit.id ∈ priority_unstudied_genes,
                is_conserved_in_vertebrates = hit.id ∈ vertebrate_genes,
                bit_score = hit.score,
                funfam_inclusion_threshold = funfam_inclusion_threshold,
                is_in_funfam = is_in_funfam,
                funfam_inclusion_delta = is_in_funfam ? 0. : funfam_inclusion_threshold - hit.score,
                )
            )
        end
    end
    return predictions
end

predictions = funfam_goterm_predictions(
    hits,
    funfam_goterms_dict,
    ontology,
    annotations,
    tcs_dict,
    priority_unstudied_genes,
    vertebrate_genes;
    include_uniprotkb_kw_iea=true
    )

# Genes with predictions
predictions_set = Set(map(x->(x[1], x[2]), predictions))
Set(map(x->x[1], predictions))

df = DataFrame(predictions)
@show(df[1:5, :])

# TEMP only use BP ontology terms
df = df[df[:,:ontology] .== "biological_process", :]

# NB sorting the dataframe before removing duplicate rows ensures that we always keep the
# annotation that has the highest bit_score
sort!(df, :bit_score, rev=true)
unique!(df, [:identifier, :goterm])

# counters
counter(df[:,:ontology])
# counter(df[:,:is_goterm_in_slim_or_descendant_of_slim_term])
counter(df[:,:is_in_funfam])

# fraction of these predictions below inclusion threshold
let x = df[:,:is_in_funfam]
    count(x), count(x) / length(x)
end

# fraction of these predictions that are already known
let df = df[df[:,:is_in_funfam] .== true, :]
    x = df[:,:is_in_pombe_annotations]
    count(x), count(x) / length(x)
end

# predcitions for priority unstudied genes
let df = @where(df, :is_priority_unstudied_gene .== true)
    x = df[:,:is_in_funfam]
    count(x), count(x) / length(x)
end

# remove predictions already in pombe annotations
df = @where(df, :is_in_pombe_annotations .== false)

# write results
CSV.write(joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "funfam", "funfam_goterm_predictions.csv"), df)


# new_predictions = CSV.read(joinpath(ENV["POMBEAGEINGGENES"], "Scripts", "ml", "gene_goterm_predictions.csv"))
# new_predictions[!, :goterm] = map(x->replace(x, "GO"=>"GO:"), new_predictions[:, :goterm])
# new_predictions_set = Set(zip(new_predictions.identifier, new_predictions.goterm))
#
# # Common predictions
# new_predictions_set ∩ predictions_set
#
# setdiff(new_predictions_set, predictions_set)
#
# # Common IDs
# Set(new_predictions[:, :identifier]) ∩ Set(map(x->x[1], predictions))
#
# # Priority unstudied genes
# priority_unstudied_genes = Set(CSV.read(joinpath(ENV["POMBEAGEINGGENES"], "Data", "priority_unstudied_genes.tsv"))[:,1])
#
# priority_unstudied_genes ∩ Set(map(x->x[1], predictions)) # FF
# priority_unstudied_genes ∩ Set(new_predictions[:, :identifier]) # RF
#
# filter(x->x[1] in priority_unstudied_genes, predictions_set)

using PombeAgeingGenes, Plots, OBOParse, .GeneOntology, DataStructures, CodecZlib

plotlyjs()

go = load(GO)

ontologies = ["biological_process", "molecular_function", "cellular_component"]

ec_classes = [:experimental, :curated, :automatic, :bad]

ids = Set(protein_coding_genes())

annotations = load(GOAnnotations; ecs=EVIDENCE_CODES_ALL)
filter!(x->x.id in ids, annotations)

function get_annotations_per_ec_class(annotations, terms_in_ontology)
    xs = filter(x->x.go in terms_in_ontology, annotations)
    d = Dict{String,Set{String}}()
    seen = Set{String}() # fast and frugal tree
    for class in ec_classes
        ecs = Set(GeneOntology.EVIDENCE_CODES[class])
        ys = filter(x->x.ec in ecs, xs)
        ids = Set(getfield.(ys, :id))
        setdiff!(ids, seen)
        d[string(class)] = ids
        union!(seen, ids)
    end
    d["none"] = setdiff(ids, seen)
    return d
end

function get_annotations_per_ontology(annotations, ontology)
    terms_in_ontology = Set{String}()
    for (k,v) in go.terms
        if v.namespace == ontology
            push!(terms_in_ontology, k)
        end
    end
    return get_annotations_per_ec_class(annotations, terms_in_ontology)
end

function faft!(A, i, d)
    h = 0 # calculate the cumulative sum
    for (j, class) in enumerate(ec_classes)
        h += length(d[string(class)])
        A[i,j] = h
    end
    return A
end

# Fast and frugal tree matrix
A = fill(length(ids), 4, 5)

# All ontologies
d = get_annotations_per_ec_class(annotations, Set(keys(go.terms)))
faft!(A, 1, d)

# Individual ontologies
for (i, ontology) in enumerate(ontologies)
    d = get_annotations_per_ontology(annotations, ontology)
    faft!(A, i+1, d)
end

bar(
    A[:,end:-1:1],
    xlabel="Ontologies",
    ylabel="Genes",
    xticks=(1:size(A,1), ["All"; uppercasefirst.(replace.(ontologies, "_"=>" "))]),
    c=reshape(repeat([:grey80; PlotThemes.wong_palette[[end,1,3,5]]], inner=size(A,1)), size(A)...),
    # c=reshape(repeat([:grey90, :grey70, :grey50, :grey30, :black], inner=size(A,1)), size(A)...),
    labels=permutedims(uppercasefirst.(reverse([string.(ec_classes); "None"]))),
    lw=0,
    grid=:y,
    size=(600,300),
)

savefig(joinpath(@__DIR__, "evidence_codes_of_pombe_annotations.pdf"))

# Thesis

# All ontologies together
d1 = get_annotations_per_ontology(annotations, "biological_process")
d2 = get_annotations_per_ontology(annotations, "molecular_function")
d3 = get_annotations_per_ontology(annotations, "cellular_component")
# BP, MF and CC
length(d1["experimental"] ∪ d2["experimental"] ∪ d3["experimental"])
length(d1["experimental"] ∪ d2["experimental"] ∪ d3["experimental"]) / 5137
# BP and MF
length(d1["experimental"] ∪ d2["experimental"])
length(d1["experimental"] ∪ d2["experimental"]) / 5137

# Ontologies separately

# BP
d = get_annotations_per_ontology(annotations, "biological_process")
length(d["curated"])+length(d["experimental"])
(length(d["curated"])+length(d["experimental"])) / 5137
length(d["none"]) + length(d["bad"])
(length(d["none"]) + length(d["bad"])) / 5137

# MF
d = get_annotations_per_ontology(annotations, "molecular_function")
length(d["curated"])+length(d["experimental"])
(length(d["curated"])+length(d["experimental"])) / 5137
length(d["none"]) + length(d["bad"])
(length(d["none"]) + length(d["bad"])) / 5137

# CC
d = get_annotations_per_ontology(annotations, "cellular_component")
length(d["experimental"])
length(d["experimental"]) / 5137
length(d["none"]) + length(d["bad"])
(length(d["none"]) + length(d["bad"])) / 5137

# Experimental CC annotations from high-throughput studies
terms_in_ontology = Set([k for (k,v) in go.terms if v.namespace == "cellular_component"])
d = Dict{String,Set{String}}(x => Set() for x in ["ht", "lt"])
for a in load(GOAnnotations; ecs=EVIDENCE_CODES[:experimental])
    if a.go in terms_in_ontology
        if a.ec in Set(["HTP", "HDA", "HMP", "HGI", "HEP"])
            push!(d["ht"], a.id)
        else
            push!(d["lt"], a.id)
        end
    end
end
setdiff!(d["ht"], d["lt"])
d

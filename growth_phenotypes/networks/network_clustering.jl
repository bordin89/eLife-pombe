using PombeAgeingGenes

using Statistics, DataFrames, LightGraphs, Plots, LinearAlgebra, Clustering,
	HypothesisTests, MultipleTesting

# Graph idea

function threshold_edges(A; edge_threshold)
    A = copy(A)
    A[A .< edge_threshold] .= zero(eltype(A))
    return A
end

function jaccard(x::Set{T}, y::Set{T}) where T
    length(x) == 0 && length(y) == 0 && return 1.
    sx, sy, si = 0, length(y), 0
	for el in x
		sx += 1
		if el ∈ y
			si += 1
		end
	end
	su = sx + sy - si
	return si / su
end

# function jaccard(x::Set{T}, y::Set{T}) where T
#     length(x) == 0 && length(y) == 0 && return 1.
#     return length(x ∩ y) / length(x ∪ y)
# end

function jaccard(g::Graph)
    n = nv(g)
    A = zeros(n, n)
    Threads.@threads for j in 1:n
        neighbors_j = Set(neighbors(g, j))
        for i in 1:j-1
            neighbors_i = Set(neighbors(g, i))
            A[i,j] = jaccard(neighbors_i, neighbors_j)
        end
    end
    return A .+ A'
end

mutable struct OverrepresentationTest
	p::Float64
	goterm::Symbol
	cluster_id::Int64
end

function overrepresentation_test(g, clust, cluster_id, annotations, goterm)
	# GO term annotations for vertices not/in cluster
	vertices_in_cluster = vertices(g)[clust.assignments .== cluster_id]
	not_in_cluster = annotations[Not(vertices_in_cluster), goterm]
	in_cluster = annotations[vertices_in_cluster, goterm]

	# Fisher's one-sided p-value that a cluster is enriched for a GO term
	#
	# -  X1 X2   X1 = not in cluster, X2 = in cluster
	# –– –– ––
	# Y1 a  b    Y1 = doesn't have term
	# Y2 c  d    Y2 = has term
	a = count(not_in_cluster .== 0)
	c = count(not_in_cluster .== 1)
	b = count(in_cluster .== 0)
	d = count(in_cluster .== 1)
	p = pvalue(FisherExactTest(a, b, c, d))#, tail=:right)
	return OverrepresentationTest(p, goterm, cluster_id)
end

function overrepresentation_tests(g, clust, annotations)
	tests = OverrepresentationTest[]
	cluster_ids = sort(unique(clust.assignments))
	goterms = names(annotations[:, Not(:id)])

	for cluster_id in cluster_ids, goterm in goterms
		push!(tests, overrepresentation_test(g, clust, cluster_id, annotations, goterm))
	end

	# Multiple testing correction
	setfield!.(tests, :p, adjust(getfield.(tests, :p), BenjaminiHochberg()))

	return tests
end

df = load(GrowthPhenotypesWideform)
ids = df[:, :id]

@assert issorted(df[:, :id])

A = Matrix(df[:, 2:end])

A = cor(A')
# clusteredheatmap(A, c=:RdBu_r, clims=(-1,1))

# TODO might have to tune the threshold
g = Graph(threshold_edges(A-I, edge_threshold=.65))

# connected_components(g)
# dd = degree(g)
# histogram(dd)

@time j = jaccard(g)
@assert issymmetric(j)

clust = dbscan(1 .- j, .7, 2)
unique(clust.assignments)

annotations = load(GeneOntology.GOSlimTargets)

tests = overrepresentation_tests(g, clust, annotations)

corr = filter(x->x.p < 1, tests)


x = ids[clust.assignments .== 0]
for i = x
	println(i)
end

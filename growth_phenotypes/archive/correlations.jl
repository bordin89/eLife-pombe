#=
Seaborn clustermaps of srains and conditions
=#
include("common.jl")

using Statistics, PyCall

@pyimport seaborn as sns
@pyimport matplotlib.pyplot as plt

plotlyjs()

A = df[:size]

@where(df, :size .== 0)

@where(df, :id .== "SPAC20H4.08", :condition .== "YES_Xilose_2_percent_0.1_glucose")


tmp = by(df, [:id, :condition], x->size(x,1))

sort(tmp, :x1)

tmp2 = unstack(tmp, :id, :condition, :x1)

tmp3 = convert(Matrix, tmp2[2:end])
tmp4 = log10.(coalesce.(tmp3, 0) .+1)

sns.clustermap(tmp4, cmap="viridis")
plt.show()

heatmap(tmp4)

histogram(@where(df, :id .== "SPBC29B5.01", :condition .== "EMM_3AT")[:size])

tmp = sort(df, :size, rev=true)

@show tmp[1:10, :][[:id, :condition, :size]]


@show sort(A)[end-10:end]
count(A .== 0.)

extrema(A)
B = A
fig = histogram(B, xlabel="Colony size", ylabel="Frequency")
savefig(fig, "colony_size_histogram_all.pdf")

histogram(B[B .< 10], xlabel="Colony size", ylabel="Frequency")
savefig("colony_size_histogram_less10.pdf")

histogram(B[0 .< B .< 10], xlabel="Colony size", ylabel="Frequency")
savefig("colony_size_histogram_all_0_between_10.pdf")

tmp = load(GrowthPhenotypesNoOutliersWideform)
A = convert(Matrix{Float64}, tmp[2:end])

dropmissing!(tmp)


using Colors
colormap("RdBu", 1000, mid=.)
heatmap(log10.(A) .+ .4, size=(1000,1000), c=:RdBu)

@show sort(A[:])
B = log10.(A)
heatmap(log10.(A), c=colormap("RdBu", 1000, mid=.58))
savefig("heatmap_strains_conditions.pdf")

x = string.(tmp[1])
y = string.(names(tmp)[2:end])
B = cor(A)
extrema(B)
C = cor(A')

ax = sns.clustermap(B, xticklabels=y, cmap="viridis")
plt.show(ax)

ax = sns.clustermap(C) #, xticklabels=x)
plt.show(ax)

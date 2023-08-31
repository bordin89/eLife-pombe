include("common.jl")

# Pairplots
# df = loaddf()
df = @where(df, :condition .== "YES_32")
removeoutliers!(df)
df

df[:hashstring] = string.(df[:hash])

df = unstack(df, :id, :hash, :size)
df = map(x->coalesce.(x, 0.), eachcol(df))

@df df cornerplot(cols(2:8), markercolor=:Blues, compact=true, size=(1000,1000), xlims=(-.2,2), ylims=(-.2,2), xrotation=45, histpct=0.00001, linecolor=:black)

# Heatmap
gr()
df = loaddf()
df = @where(df, :condition .== "YES_32")
removeoutliers!(df)
df = unstack(df, :id, :hash, :size)
df = map(x->coalesce.(x, NaN), eachcol(df))
A = hcat([convert(Array, df[i]) for i in 2:ncol(df)]...)
heatmap(A, color=:plasma, xlabel="Repeat", ylabel="KO")

# Positional effects
df = loaddf()
df = @where(df, :condition .== "YES_32")
h1 = df[:hash][1]
h2 = df[:hash][end]
df1 = @select(@where(df, :hash .== h1), :id, :size, :assay_plate_col, :assay_plate_row, :assay_plate)
df2 = @select(@where(df, :hash .== h2), :id, :size, :assay_plate_col, :assay_plate_row, :assay_plate)
df = join(df1, df2, on=[:id, :assay_plate_col, :assay_plate_row, :assay_plate], kind=:inner)

cols = distinguishable_colors(length(controls)+2, colorant"white")[2:end]
opts = (α=.5, markerstrokewidth=false)
@df df scatter(:size, :size_1, c=:grey75, legend=false, label="Other strains", xlabel="Repeat x", ylabel="Repeat y"; opts...)
for (control, c) = zip(controls, cols)
    tmp = @where(df, :id .== control)
    @df tmp scatter!(:size, :size_1, label=control, c=c; opts...)
end
plot!([0,1.5], [0,1.5], c=:black, label="y=x")
plot!(legend=:bottomright)
savefig("repeat_repeat_scatter.pdf")

# Variation between repeats
df = loaddf()

df = @where(df, :condition .== "YES_32")
df[:hash] = string.(df[:hash])
countmap(df[:hash])

@df df violin(:hash, :size, xticks=false, xlabel="YES_32 repeats", ylabel="Relative size", legend=false)

dfno = removeoutliers(df)

@df dfno violin(:hash, :size, xticks=false, xlabel="YES_32 repeats", ylabel="Relative size", legend=false)

# Residuals
using Combinatorics, Plots.PlotMeasures

function residual(x::T, y::T; a::T=one(T), b::T=zero(T)) where T<:Number
    r = y - (a*x + b)
end

df = loaddf()
df = @where(df, :condition .== "YES_32")
removeoutliers!(df)
df = unstack(df, :id, :hash, :size)

function residuals_(df, id)
    df = @where(df, :id .== id)
    cols = names(df)[2:end]  # remove :id

    res = Float64[]

    for (i, j) = combinations(cols, 2)
        i == j && continue

        r1, r2 = Symbol.((i, j))
        vs = df[r1][1], df[r2][1]

        any(ismissing.(vs)) && continue

        push!(res, abs(residual(vs...)))
    end
    res
end

res_controls = [residuals_(df, i) for i = controls]
violin(res_controls, xticks=(1:length(controls), controls), xr=45, legend=false, ylabel="Absolute residual", c=false, bottom_margin=2mm)
savefig("residuals_violin.pdf")

sig = []
for (i,j) = combinations(1:length(res_controls), 2)
    t = MannWhitneyUTest(res_controls[i], res_controls[j])
    p = pvalue(t)
    p < 0.05 && push!(sig, (controls[i], controls[j], p))
end
sig

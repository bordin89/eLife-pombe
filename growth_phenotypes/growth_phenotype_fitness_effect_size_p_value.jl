#=
Export a table for Maria's paper that includes growth phenotype:
- fitness
- effect size
- P value
=#

using PombeAgeingGenes
using DataFrames
using CodecZlib
using Statistics
using DataFramesMeta
using HypothesisTests
using MultipleTesting
using CSV
using UnicodePlots
using DataStructures


df_clean = load(GrowthPhenotypes)[:, [:condition, :id, :size]]

df_fitness = by(df_clean, [:id, :condition]) do g
    xs = g.size
    return (
        fitness_mean = round(mean(xs), digits=3),
        fitness_median = round(median(xs), digits=3),
    )
end


df_normalised = load(GrowthPhenotypesNoOutliers)

df_effect_size = by(df_normalised, [:id, :condition]) do g
    xs = g.size
    return (
        effect_size = round(xs[argmax(abs.(xs))], digits=3),
        effect_size_mean = round(mean(xs), digits=3),
        effect_size_median = round(median(xs), digits=3),
    )
end

# P value for `fitness`

function null_distributions(df)
    d = DefaultDict{String, Dict{String, Vector{Float64}}}(Dict{String, Vector{Float64}})

    for g in groupby(df, [:id, :condition])
        id = g.id[1]
        condition = g.condition[1]

        d[startswith(condition, "YES") ? "YES_32" : "EMM_32"][id] = g.size
    end
    
    return d
end

function calculate_p_values(df)
    d = null_distributions(df)

    return by(df, [:id, :condition]) do g
        id = g.id[1]
        condition = g.condition[1]
        control_condition = startswith(condition, "YES") ? "YES_32" : "EMM_32"
        null = d[control_condition][id]
        xs = g.size

        p = if length(xs) ≤ 1 || length(null) ≤ 1 || iszero(var(xs))
            missing
        else
            pvalue(UnequalVarianceTTest(xs, null))
        end
  
        return (p_value = p,)
    end
end

function correct_p_values(df)
    df_sub = dropmissing(df, :p_value)

    df_sub[!, :p_corrected] = fill(NaN, size(df_sub, 1))

    for g = groupby(df_sub, :condition)
        g.p_corrected = adjust(collect(g.p_value), BenjaminiHochberg())
    end

    join_on = [:id, :condition, :p_value]
    for col in [:fitness_mean, :fitness_median]
        if col in names(df)
            push!(join_on, col)
        end
    end

    df = join(
        df,
        df_sub,
        on = join_on,
        kind = :outer,
    )

    for col in [:p_value, :p_corrected]
        df[!, col] = coalesce.(df[!, col], NaN)
    end
    
    dropmissing!(df)

    return df
end

df = load(GrowthPhenotypes)[:, [:id, :condition, :size]]

# df_p = calculate_p_values(df[1:1_000_000, :])
df_p = calculate_p_values(df)
df_p = correct_p_values(df_p)

count(isnan, df_p.p_value)
count(isnan, df_p.p_corrected)

count(<(.05), filter(!isnan, df_p.p_value))
count(<(.05), filter(!isnan, df_p.p_corrected))

histogram(filter(!isnan, df_p.p_value))
histogram(filter(!isnan, df_p.p_corrected))

df = join(df_fitness, df_effect_size, on=[:id, :condition], kind=:inner)
df = join(df, df_p, on=[:id, :condition], kind=:inner)

DataFrames.rename!(df, "id"=>"gene_id")
CSV.write(joinpath(@__DIR__, "growth_phenotype_fitness_effect_size_p_value.csv"), df)


#=
Control conditions

Table with the fitness for the mutants under the control conditions.
This should be a separated table because you don’t normalize against control,
but this way you can see what genes have an effect in growth in unstressed conditions.
=#

get_controls(df) = df[map(x->x in Set(["YES_32", "EMM_32"]), df.condition), :]

df = get_controls(load(GrowthPhenotypes)[:, [:condition, :id, :size]])

function calculate_p_values(df)
    d = null_distributions(df)

    return by(df, [:id, :condition]) do g
        id = g.id[1]
        condition = g.condition[1]
        control_condition = startswith(condition, "YES") ? "YES_32" : "EMM_32"
        null = d[control_condition]["wt"]
        xs = g.size

        p = if length(xs) ≤ 1 || length(null) ≤ 1 || iszero(var(xs))
            missing
        else
            pvalue(UnequalVarianceTTest(xs, null))
        end

        return (
            fitness_mean = round(mean(xs), digits=3),
            fitness_median = round(median(xs), digits=3),
            p_value = p,
        )
    end
end

df_p = calculate_p_values(df)
df_p = correct_p_values(df_p)

count(isnan, df_p.p_value)
count(isnan, df_p.p_corrected)

count(<(.05), filter(!isnan, df_p.p_value))
count(<(.05), filter(!isnan, df_p.p_corrected))

histogram(filter(!isnan, df_p.p_value))
histogram(filter(!isnan, df_p.p_corrected))

df = df_p
DataFrames.rename!(df, "id"=>"gene_id")
CSV.write(joinpath(@__DIR__, "growth_phenotype_control_conditions_fitness.csv"), df)





# # NULL_YES = @where(df, :condition .== "YES_32")
# # NULL_EMM = @where(df, :condition .== "EMM_32")

# # df_p = by(df, :id) do g
# #     id = g.id[1]
# #     condition = g.condition[1]

# #     # null = @where(
# #     #     startswith(condition, "YES") ? NULL_YES : NULL_EMM,
# #     #     :id .== id
# #     # ).size

# #     return by(g, :condition) do h
# #         null = @where(
# #             startswith(condition, "YES") ? NULL_YES : NULL_EMM,
# #             :id .== id
# #         ).size

# #         p = try
# #             pvalue(UnequalVarianceTTest(h.size, null))
# #         catch ArgumentError
# #             1.
# #         end

# #         return (
# #             p_value = p,
# #         )
# #     end
# # end


# # function calculate_p_values(df)
# #     # NULL_YES = @where(df, :condition .== "YES_32")
# #     # NULL_EMM = @where(df, :condition .== "EMM_32")

# #     return by(df, [:id, :condition]) do g
# #         id = g.id[1]
# #         condition = g.condition[1]

# #         # df_condition = startswith(condition, "YES") ? NULL_YES : NULL_EMM
# #         control_condition = startswith(condition, "YES") ? "YES_32" : "EMM_32"

# #         df_null = @views df[(df.condition .== control_condition) .& (df.id .== id), :]
# #         # null = @views df_condition[df_condition.id .== id, :size]

# #         # p = try
# #         #     pvalue(UnequalVarianceTTest(g.size, df_null.size))
# #         # catch ArgumentError
# #         #     1.
# #         # end
# #         pvalue(UnequalVarianceTTest(g.size, df_null.size))

# #         return (p_value = p,)
# #     end
# # end
# d = DefaultDict{String, Dict{String, Vector{Float64}}}(Dict{String, Vector{Float64}})

# d["a"]["b"] = [1.]
# d



#=

# P value for `effect_size`

df = load(GrowthPhenotypesNoOutliers)

NULL_YES = @where(df, :condition .== "YES_32")
NULL_EMM = @where(df, :condition .== "EMM_32")

df_p = by(@where(df, :condition .!= "YES_32", :condition .!= "EMM_32"), :id) do g
    id = g.id[1]
    condition = g.condition[1]

    null = @where(
        startswith(condition, "YES") ? NULL_YES : NULL_EMM,
        :id .== id).size

    return by(g, :condition) do h
        p = try
            pvalue(UnequalVarianceTTest(h.size, null))
        catch ArgumentError
            1.
        end

        return (
            p_value = p,
        )
    end
end
=#

#====
# NULL_YES = @where(df, :id .== "wt", :condition .== "YES_32").size
# NULL_EMM = @where(df, :id .== "wt", :condition .== "EMM_32").size

# df = by(df, [:id, :condition]) do g
#     xs = g.size

#     null = startswith(g.condition[1], "YES") ? NULL_YES : NULL_EMM

#     p = try
#         pvalue(UnequalVarianceTTest(xs, null))
#     catch ArgumentError
#         1.
#     end

#     return (
#         fitness_mean = round(mean(xs), digits=3),
#         fitness_median = round(median(xs), digits=3),
#         p_value = p,
#     )
# end

# df[!, :p_corrected] .= NaN
# df[!, :p_corrected] = fill(NaN, size(df, 1))

# for g = groupby(df, :condition)
#     g[:, :p_corrected] = adjust(collect(g.p_value), BenjaminiHochberg())
# end

=#
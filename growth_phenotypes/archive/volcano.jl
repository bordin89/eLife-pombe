include("common.jl")

# Volcano plot
df = loaddf()
df = @where(df, :condition .== "YES_32")

df_p = by(df, :condition) do x
    wt_dist = @where(x, :id .== "wt")[:size]
    by(x, :id) do y
        s = size(y, 1)
        if s > 1
            pvalue(UnequalVarianceTTest(y[:size], wt_dist))  # 2-sided
        else
            missing
        end
    end
end
rename!(df_p, :x1=>:p)
df_p[:bf] = bonferroni(df_p[:p])
histogram(df_p[:p])

df_fc = by(df, :condition) do x
    wt_dist = @where(x, :id .== "wt")[:size]
    wt = mean(wt_dist)
    by(x, :id) do y
        mean(y[:size]) / wt
    end
end
rename!(df_fc, :x1=>:fc)
histogram(df_fc[:fc])

df_j = join(df_p, df_fc, on=[:id, :condition])
df_j[:logp] = -log10.(df_j[:p])

plotlyjs()
fc = .15
opts = (α=.25, markerstrokewidth=false)
@df @where(df_j, :bf.≥.05) scatter(:fc, :logp, c=:grey50; opts...)
@df @where(df_j, :bf.<.05, :fc.≥1-fc, :fc.≤1+fc) scatter!(:fc, :logp, c=:grey50; opts...)
@df @where(df_j, :bf.<.05, :fc.<1-fc) scatter!(:fc, :logp, c=:red; opts...)
@df @where(df_j, :bf.<.05, :fc.>1+fc) scatter!(:fc, :logp, c=:green; opts...)
plot!(ylabel="-log10(p)", xlabel="FC")

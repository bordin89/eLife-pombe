#=
Plot bootstrap curves demonstrating how much additional information is provided by
additional growth phenotypes.
=#

using StatsBase, PombeAgeingGenes, DecisionTree, Plots

X = rand(20, 1000)

X, Y = load(ML, Matrix)

iterrange = 1:5:nconditions-1

function bootstrap(X, iterrange; nboot=10)
    nconditions, nids = size(X)
    results = Vector{Float64}[]

    for (i, n) = enumerate(iterrange)
        push!(results, Float64[])

        for b = 1:nboot
            feature_idx = sample(1:nconditions, n, replace=false)
            target_idx = sample(setdiff(1:nconditions, feature_idx))

            x = X[feature_idx, :]
            # y = X[target_idx, :]
            y = Y[1, :]

            ŷs = Vector{Float64}[]
            ys = Vector{Float64}[]

            k = 5

            for ((Xtrain, ytrain), (Xtest, ytest)) = kfolds((x, y), k)
                # model = DecisionTreeRegressor() #; kwargs...)
                model = RandomForestRegressor(n)
                DecisionTree.fit!(model, collect(Xtrain'), collect(ytrain))
                ŷ = DecisionTree.predict(model, collect(Xtest'))
                push!(ŷs, ŷ)
                push!(ys, ytest)
            end

            d = msd.(ŷs, ys)
            push!(results[i], mean(d))
        end
    end

    results
end

results = bootstrap(X, iterrange, nboot=2)

plot(iterrange,
    mean.(results),
    fillrange=(maximum.(results), minimum.(results)),
    xticks=1:2:nconditions-1,
    lc=:black,
    legend=false,
    xlabel="# conditions",
    ylabel="Mean squared error",
    fillalpha=.5,
    )


fig = plot()

for (i,n) = enumerate(iterrange)
    ys = results[i]
    xs = repeat([n], length(ys))
    scatter!(fig, xs, ys)
end

plot!(legend=false)

n = 10
feature_idx = sample(1:nconditions, n, replace=false)
x = X[feature_idx, :]
target_idx = sample(setdiff(1:nconditions, feature_idx))
y = X[target_idx, :]

ŷs = Vector{Float64}[]
ys = Vector{Float64}[]

k = 5
Model = DecisionTreeRegressor

for ((Xtrain, ytrain), (Xtest, ytest)) = kfolds((x, y), k)
    model = DecisionTreeRegressor() #; kwargs...)
    DecisionTree.fit!(model, collect(Xtrain'), collect(ytrain))
    ŷ = DecisionTree.predict(model, collect(Xtest'))
    push!(ŷs, ŷ)
    push!(ys, ytest)
end

d = msd.(ŷs, ys)

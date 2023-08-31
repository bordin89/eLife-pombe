using PombeAgeingGenes, DecisionTree, Plots

X, Y = load(ML, Matrix, deepNF_embeddings=true)
X = permutedims(X)
y = Y[5,:]
y = float.(y)

model = RandomForestRegressor(n_trees=1000)

DecisionTree.fit!(model, X, y)

ŷ = DecisionTree.predict(model, X)
pr = PR(ŷ, y)
@plotpr pr

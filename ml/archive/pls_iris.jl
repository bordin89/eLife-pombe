#=
Test PLS on Iris data set.
=#
using PombeAgeingGenes, RDatasets, Plots
iris = dataset("datasets", "iris")
indices = collect(1:size(iris, 1))
shuffle!(indices)

Y = string.(iris[5])[indices]
Y = labelencode(labelmap(Y), Y)
Y = Array{Float64}(sparse(1:length(Y), Y, true))
X = convert(Array, iris[1:4])[indices, :]
X, Xtest = X[1:120,:], X[121:end,:]
Y, Ytest = Y[1:120,:], Y[121:end,:]

X = collect(X')
Y = collect(Y')

Xtest = collect(Xtest')
Ytest = collect(Ytest')

M = PLS(X,Y,3)
p = xtransform(M, X)
y = evaluate(M, Xtest)

scatter(p[1,:], p[2,:], c=getindex.(argmax(Y, dims=1), 1), legend=true)

q = xtransform(M, Xtest)

scatter(q[1,:], q[2,:], c=getindex.(argmax(Ytest, dims=1), 1), legend=true)

testY = evaluate(M, Xtest)
Ytest
testerror = norm(Ytest - testY) / sqrt(size(Ytest, 1))

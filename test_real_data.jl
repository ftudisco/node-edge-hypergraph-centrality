using LinearAlgebra
using SparseArrays
using PyPlot

include("centrality_tools.jl")
include("real_data_tools.jl")

B,W,edges = read_tags_data("tags-math-sx");
labels = read_tags_data_node_labels("tags-math-sx")

n,m = size(B)

centralities = Dict()

αs = [-10 -1 ]

a = 2
b = 10
f = x -> x.^a
g = x -> x.^(1/a)
ϕ = x -> x.^b
ψ = x -> x.^(1/b)

# B,labels,names = read_walmart_data();
x,y = compute_centrality(B,W,f,g,ϕ,ψ,maxiter=500);

centralities[(a,b)] = Dict("x"=>x, "y"=>y)

xid = sortperm(vec(x),rev=true);
xx = x[xid]

yid = sortperm(vec(y),rev=true);
yy = y[yid]

w = vec(diag(W))
ww = w[yid]

δ = vec(sum(B,dims=1))
δid = sortperm(δ,rev=true)
δδ = δ[yid]

figure()
subplot(121)
semilogy(xx[xx.>1e-10],"o")
subplot(122)
semilogy(yy,"x")
semilogy(vec(1:m)[δδ.==1],yy[δδ.==1],"x")
# subplot(133)
# # PyPlot.semilogy(δδ[yy.>1e-10],"+")
# plot(δ./maximum(δ),y./maximum(y),"+")

yy = yy./maximum(yy)
figure()
subplot(231); hist(yy[δδ.==1],bins=50,range=(0,1),density=true)
subplot(232); hist(yy[δδ.==2],bins=50,range=(0,1),density=true)
subplot(233); hist(yy[δδ.==3],bins=50,range=(0,1),density=true)
subplot(234); hist(yy[δδ.==4],bins=50,range=(0,1),density=true)
subplot(235); hist(yy[δδ.==5],bins=50,range=(0,1),density=true)

figure()
loglog(w./maximum(w),y./maximum(y),"+")
xlabel("edge weight")
ylabel("edge centrality")
figure()
plot(δ./maximum(δ),y./maximum(y),"+")


# print top k edges
k = 35
top_k_edges = Dict()
for e in 1:k
    top_k_edges[e] = Set(findnz(B[:,yid[e]])[1])
end

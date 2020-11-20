using LinearAlgebra
using SparseArrays
# using PyPlot
using PyCall
using StatsBase
wordcloud = pyimport("wordcloud");
sb = pyimport("seaborn")

include("centrality_tools.jl")
include("real_data_tools.jl")

B,W,edges = read_tags_data("tags-math-sx");
labels = read_tags_data_node_labels("tags-math-sx")

# B,labels,names = read_walmart_data();

n,m = size(B)

centralities = Dict()

αs = [-10 -1 2 10]
βs = [-10 -1 2 10]

for a in αs
    for b in βs
        f = x -> x.^a
        g = x -> x.^(1/a)
        ϕ = x -> x.^b
        ψ = x -> x.^(1/b)
        println("a = $a, b = $b")
        x,y = compute_centrality(B,W,f,g,ϕ,ψ,maxiter=500);
        centralities[(a,b)] = Dict("x"=>x, "y"=>y)
    end
end


figure()
for (i,a) in enumerate(αs)
    for (j,b) in enumerate(βs)
        subplot(4,4,4*(i-1)+j)
        x = centralities[(a,b)]["x"]
        y = centralities[(a,b)]["y"]
        yid = sortperm(vec(y),rev=true);
        yy = y[yid]
        yy = yy./maximum(yy)
        semilogy(yy,"x")
        title("a=$a, b=$b")
        tight_layout()
    end
end

figure()
x = centralities[(10,2)]["x"]
y = centralities[(10,2)]["y"]
x = x./maximum(x)
y = y./maximum(y)
for (i,a) in enumerate([-10 -1 2])
    xx = centralities[(a,2)]["x"]
    yy = centralities[(a,2)]["y"]
    xx = xx./maximum(xx)
    yy = yy./maximum(yy)
    subplot(2,3,i)
    plot(x,xx,"o")
    subplot(2,3,3+i)
    plot(y,yy,"o")
end


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
subplot(231); sb.histplot(data=yy[δδ.==1],element="step",stat="probability")
subplot(232); sb.histplot(data=yy[δδ.==2],element="step",stat="probability") #hist(yy[δδ.==2],bins=50,range=(0,1),density=true)
subplot(233); sb.histplot(data=yy[δδ.==3],element="step",stat="probability") #hist(yy[δδ.==3],bins=50,range=(0,1),density=true)
subplot(234); sb.histplot(data=yy[δδ.==5],element="step",stat="probability") #hist(yy[δδ.==4],bins=50,range=(0,1),density=true)
subplot(235); sb.histplot(data=yy[δδ.==6],element="step",stat="probability") #hist(yy[δδ.==5],bins=50,range=(0,1),density=true)

figure()
loglog(w./maximum(w),y./maximum(y),"+")
xlabel("edge weight")
ylabel("edge centrality")
figure()
plot(δ./maximum(δ),y./maximum(y),"+")


# print top k edges
k = 100
top_k_edges = Dict()
top_k_edges_all = []
text = ""
for e in 1:k
    idx_nodes = findnz(B[:,yid[e]])[1]
    lab_nodes = Set()
    for (i,id) in enumerate(idx_nodes)
        push!(lab_nodes,labels[id])
        push!(top_k_edges_all, labels[id])
        global text *= labels[id]*", "
    end
    top_k_edges[e] = lab_nodes #Set(findnz(B[:,yid[e]])[1])
end
text = text[1:end-2];


freq = StatsBase.countmap(top_k_edges_all);
wc = wordcloud.WordCloud(background_color="white").generate_from_frequencies(freq);
# snp = pyimport("seaborn")

figure()
subplot(121)
PyPlot.imshow(wc, interpolation="spline16")
PyPlot.axis("off")
subplot(122)
a = sb.histplot(data=top_k_edges_all)
xticks(rotation=90)
tight_layout()

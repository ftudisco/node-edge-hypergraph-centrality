using LinearAlgebra
using SparseArrays
using PyPlot
using PyCall
using StatsBase
wordcloud = pyimport("wordcloud");
sb = pyimport("seaborn")

include("centrality_tools.jl")
include("real_data_tools.jl")


# given two vectors u and v returns the sequence of correlation values between top k ranked entries of u and v
function correlation_map(u,v)
    m = length(u);
    spearman = [];
    kendall = [];
    idu = sortperm(u,rev=true)
    for i in 1:floor(log2(m))
        k = Int64(2^i)
        append!(  spearman , corspearman(vec(u[idu[1:k]]), vec(v[idu[1:k]]))  )
        append!(  kendall , corkendall(vec(u[idu[1:k]]), vec(v[idu[1:k]]))  )
        println("$k...")
    end
    spearman[isnan.(spearman)].=0
    kendall[isnan.(kendall)].=0
    return spearman, kendall
end






B,W,labels,names,edges = read_walmart_data();

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


x = centralities[(10,2)]["x"]
y = centralities[(10,2)]["y"]
x = x./maximum(x)
y = y./maximum(y)
correlations = Dict()
for a in αs
    for b in βs
        xx = centralities[(a,b)]["x"]
        xx = xx./maximum(xx)
        sp_x,ke_x = correlation_map(x,xx)
        yy = centralities[(a,b)]["y"]
        yy = yy./maximum(yy)
        sp_y,ke_y = correlation_map(y,yy)
        println("a = $a, b = $b")
        correlations[(a,b)] = Dict("sp_x"=>sp_x, "sp_y"=>sp_y, "ke_x"=>ke_x, "ke_y"=>ke_y)
    end
end


ksy = []
for i in 1:floor(log2(length(y)))
    k = Int64(2^i)
    append!(ksy,k)
end
ksx = []
for i in 1:floor(log2(length(x)))
    k = Int64(2^i)
    append!(ksx,k)
end

figure()
for (i,a) in enumerate(αs)
    for (j,b) in enumerate(βs)
        
        if (a==10 && b==2) 
                        
        else
        
            # subplot(221)
            # semilogx(correlations[(a,b)]["sp_x"], "x",label = "a=$a, b=$b")
            # ylabel("Spearman corr")
            # title("node centrality")
            subplot(121)
            semilogx(ksx,correlations[(a,b)]["ke_x"], label = "α=$a, β=$b", linewidth=2)
            ylabel("Kendall-τ correlation")
            title("node centrality")
            xlabel("Number top ranked elements")

            # subplot(222)
            # semilogx(correlations[(a,b)]["sp_y"], "x",label = "a=$a, b=$b")
            # title("edge centrality")
            subplot(122)
            semilogx(ksy,correlations[(a,b)]["ke_y"], label = "α=$a, β=$b", linewidth=2)
            title("edge centrality")
            xlabel("Number top ranked elements")
            
        end
    end
end
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
tight_layout()





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
        ylabel("sorted edge centrality")
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
    title("x centrality 10,2 vs $a,2")
    subplot(2,3,3+i)
    plot(y,yy,"o")
    title("y centrality 10,2 vs $a,2")
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
title("plot of y vs edge size")
# subplot(133)
# # PyPlot.semilogy(δδ[yy.>1e-10],"+")
# plot(δ./maximum(δ),y./maximum(y),"+")

yy = yy./maximum(yy)
figure()
subplot(231); ax = sb.histplot(data=yy[δδ.==1],element="step",stat="density"); ax.set_xlim(0, 1);
subplot(232); ax = sb.histplot(data=yy[δδ.==2],element="step",stat="density"); ax.set_xlim(0, 1);#hist(yy[δδ.==2],bins=50,range=(0,1),density=true)
subplot(233); ax = sb.histplot(data=yy[δδ.==3],element="step",stat="density"); ax.set_xlim(0, 1); #hist(yy[δδ.==3],bins=50,range=(0,1),density=true)
subplot(234); ax = sb.histplot(data=yy[δδ.==4],element="step",stat="density"); ax.set_xlim(0, 1); #hist(yy[δδ.==4],bins=50,range=(0,1),density=true)
subplot(235); ax = sb.histplot(data=yy[δδ.==5],element="step",stat="density"); ax.set_xlim(0, 1); #hist(yy[δδ.==5],bins=50,range=(0,1),density=true)

figure()
subplot(121)
loglog(w./maximum(w),y./maximum(y),"+")
xlabel("edge weight")
ylabel("edge centrality")
subplot(122)
plot(δ./maximum(δ),y./maximum(y),"+")
xlabel("edge size")
ylabel("edge centrality")

# print top k edges
k = 100
top_k_edges = Dict()
top_k_edges_all = []
top_k_edges_all_id =[]
text = ""
for e in 1:k
    idx_nodes = findnz(B[:,yid[e]])[1]
    lab_nodes = Set()
    for (i,id) in enumerate(idx_nodes)
        push!(lab_nodes,names[labels[id]])
        push!(top_k_edges_all, names[labels[id]])
        push!(top_k_edges_all_id, id)
        global text *= names[labels[id]]*", "
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

freq_ids = StatsBase.countmap(top_k_edges_all_id);
p = sortperm(collect(values(freq_ids)));
[collect(values(freq_ids))[p] collect(keys(freq_ids))[p] labels[collect(keys(freq_ids))[p]]]


#node centrality ranking vs most popular nodes in highest ranked edges
figure()
top_nodes_by_edges = collect(keys(freq_ids))[p];
m = length(top_nodes_by_edges)
top_nodes_by_centrality = xid[1:m]
sb.scatterplot(x=top_nodes_by_centrality,y=top_nodes_by_edges,size=top_nodes_by_centrality)

#edge centrality ranking vs sum of node centrality within each edge
BX = Diagonal(x)*B
sum_xcentrality_on_edges = sum(BX,dims=1)
yx = sum_xcentrality_on_edges./maximum(sum_xcentrality_on_edges)
figure()
sb.scatterplot(x=y,y=vec(yx))

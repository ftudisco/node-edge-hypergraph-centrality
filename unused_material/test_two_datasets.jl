using LinearAlgebra
using SparseArrays
using PyPlot
using PyCall
using StatsBase

wordcloud = pyimport("wordcloud");
sb = pyimport("seaborn")

include("centrality_tools.jl")
include("real_data_tools.jl")


αs = [-10 -1 2 10]
βs = [-10 -1 2 10]

###### Walmart ############ ===================================================
B,W,labels,names,edges = read_walmart_data();
n_walmart,m_walmart = size(B)
centralities_walmart = centrality_table(αs,βs,B,W,maxiter=500)

x = centralities_walmart[(10,2)]["x"]
y = centralities_walmart[(10,2)]["y"]
correlations_walmart = correlation_table(x,y,centralities_walmart)
###### Stackoverflow ######
B,W,edges = read_tags_data("tags-math-sx");
labels = read_tags_data_node_labels("tags-math-sx")
n_math,m_math = size(B)
centralities_math = centrality_table(αs,βs,B,W,maxiter=500)

x = centralities_math[(10,2)]["x"]
y = centralities_math[(10,2)]["y"]
correlations_math = correlation_table(x,y,centralities_math)
# =============================================================================



###### FIGURE 1: correlations centralities edges and nodes =====================================
linestyle_str = [
          ("solid", "solid"),      # Same as (0, ()) or "-"
          ("dotted", "dotted"),    # Same as (0, (1, 1)) or "."
          ("dashed", "dashed"),    # Same as "--"
          ("dashdot", "dashdot")]  # Same as "-."
     
kx1 = [Int64(2^i) for i in 1:floor(log2(n_walmart)) ]; append!(kx1,n_walmart)
kx2 = [Int64(2^i) for i in 1:floor(log2(n_math))    ]; append!(kx2,n_math)
ky1 = [Int64(2^i) for i in 1:floor(log2(m_walmart)) ]; append!(ky1,m_walmart)
ky2 = [Int64(2^i) for i in 1:floor(log2(m_math))    ]; append!(ky2,m_math)

figure()
for (i,(a,b)) in enumerate(keys(correlations_math))
    if (a==10 && b==2) 
                        
    else
        subplot(221)
        println(i%5+1)
        semilogx(kx1,correlations_walmart[(a,b)]["sp_x"], label="α=$a, β=$b", linewidth=2, linestyle = linestyle_str[(i-1)%4+1][2])
        
        ylabel("Kendall-τ correlation")
        title("Walmart trips, node centrality")
        
        subplot(223)
        semilogx(kx2,correlations_math[(a,b)]["sp_x"], label="α=$a, β=$b", linewidth=2, linestyle = linestyle_str[(i-1)%4+1][2])
        xlabel("Number top ranked elements")
        ylabel("Kendall-τ correlation")
        title("Math stackoverflow, node centrality")

        subplot(222)
        semilogx(ky1,correlations_walmart[(a,b)]["sp_y"], label="α=$a, β=$b", linewidth=2, linestyle =linestyle_str[(i-1)%4+1][2])
        title("Walmart trips, edge centrality")
        
        if i==length(keys(correlations_math))
            l = legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        end
        
        subplot(224)
        semilogx(ky2,correlations_math[(a,b)]["sp_y"], label="α=$a, β=$b", linewidth=2, linestyle = linestyle_str[(i-1)%4+1][2])
        xlabel("Number top ranked elements")
        title("Math stackoverflow, edge centrality")
    end
end
# legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
tight_layout()
# ========================================================================================================



# Figure 2: Walmart edge centrality vs weights and size #============================================
B,W,labels,names,edges = read_walmart_data();
w = vec(diag(W))
δ = vec(sum(B,dims=1))
x = centralities_walmart[(10,2)]["x"]
y = centralities_walmart[(10,2)]["y"]


figure()
subplot(121)
loglog(w./maximum(w),y./maximum(y),"+")
xlabel("edge weight")
ylabel("edge centrality")
subplot(122)
plot(δ./maximum(δ),y./maximum(y),"+")
xlabel("edge size")
ylabel("edge centrality")
tight_layout()
# ================================================================================================

# Figure 2: Math stack edge centrality vs weights and size #============================================
B,W,edges = read_tags_data("tags-math-sx");
w = vec(diag(W))
δ = vec(sum(B,dims=1))
x = centralities_math[(10,2)]["x"]
y = centralities_math[(10,2)]["y"]


figure()
subplot(121)
loglog(w./maximum(w),y./maximum(y),"+")
xlabel("edge weight")
ylabel("edge centrality")
subplot(122)
plot(δ./maximum(δ),y./maximum(y),"+")
xlabel("edge size")
ylabel("edge centrality")
tight_layout()  
# ================================================================================================


# Figure 3: Wordcloud of most popular nodes, WARNING: define choose B first ========================================
B,W,edges = read_tags_data("tags-math-sx");
degree = vec(sum(B,dims=2))
did = sortperm(degree,rev=true)
labels = read_tags_data_node_labels("tags-math-sx")
x = centralities_math[(10,2)]["x"]; x = x./maximum(x)
xid = sortperm(vec(x),rev=true);
y = centralities_math[(10,2)]["y"]
yid = sortperm(vec(y),rev=true);
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
        push!(lab_nodes,labels[id])
        push!(top_k_edges_all, labels[id])
        push!(top_k_edges_all_id, id)
        global text *= labels[id]*", "
    end
    top_k_edges[e] = lab_nodes #Set(findnz(B[:,yid[e]])[1])
end
text = text[1:end-2];


freq = StatsBase.countmap(top_k_edges_all);
wc = wordcloud.WordCloud(background_color="white").generate_from_frequencies(freq);
# snp = pyimport("seaborn")

figure()
PyPlot.imshow(wc, interpolation="spline16")
PyPlot.axis("off")

figure()
a = sb.histplot(data=top_k_edges_all)
xticks(rotation=90)
ylabel("Number of occurences in the top 100 edges")
tight_layout()

freq_centrality = Dict()
for i = 1 : length(values(freq))
    freq_centrality[labels[xid[i]]] = x[xid[i]]
end
wc = wordcloud.WordCloud(background_color="white").generate_from_frequencies(freq_centrality);

figure()
PyPlot.imshow(wc, interpolation="spline16")
PyPlot.axis("off")

freq_degree = Dict()
for i = 1 : length(values(freq))
    freq_degree[labels[did[i]]] = degree[did[i]]
end
wc = wordcloud.WordCloud(background_color="white").generate_from_frequencies(freq_degree);

figure()
PyPlot.imshow(wc, interpolation="spline16")
PyPlot.axis("off")
# ================================================================================================================




# Figure 5: kendall compare edge centrality vs edge count map
f_id = StatsBase.countmap(top_k_edges_all_id)
keys1 = [k for k in keys(f_id)]
sortperm_vals = sortperm([v for v in values(f_id)], rev=true)
top_edges_sorted = keys1[sortperm_vals]


degree = vec(sum(B,dims=2))
did = sortperm(degree,rev=true)
c_centrality = []
c_degree = []
for i = 1:length(top_edges_sorted)
    append!(c_centrality, corkendall(top_edges_sorted[1:i],xid[1:i]) )
    append!(c_degree, corkendall(top_edges_sorted[1:i],did[1:i]) )
end
figure()
plot(c_centrality)
plot(c_degree)















freq_ids = StatsBase.countmap(top_k_edges_all_id);
p = sortperm(collect(values(freq_ids)));
[collect(values(freq_ids))[p] collect(keys(freq_ids))[p] labels[collect(keys(freq_ids))[p]]]


#node centrality ranking vs most popular nodes in highest ranked edges
figure()
top_nodes_by_edges = collect(keys(freq_ids))[p];
m = length(top_nodes_by_edges)
top_nodes_by_centrality = xid[1:m]
sb.scatterplot(x=top_nodes_by_centrality,y=top_nodes_by_edges,size=top_nodes_by_centrality)








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

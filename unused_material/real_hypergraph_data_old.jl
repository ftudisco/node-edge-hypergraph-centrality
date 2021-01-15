using LinearAlgebra
using SparseArrays
using StatsBase


include("centrality_tools.jl")
include("real_data_tools.jl")

# choose dataset :

dataset_name = "stackoverflow" 
name_title = "Math stackexchange co-tags"

# dataset_name = "walmart"
# name_title = "Walmart trips"  



if dataset_name == "walmart"
    B,W,labels,names,edges = read_walmart_data();
elseif dataset_name == "stackoverflow"
    B,W,edges = read_tags_data("tags-math-sx");
    labels = read_tags_data_node_labels("tags-math-sx")
end

n,m = size(B)


mappings = Dict("linear"    => (x -> x, x -> x, x -> x, x -> x), 
            "log-exp"       => (x -> x, x -> x.^(1/10), x -> log.(x), x -> exp.(x)),
            "max"   => (x -> x, x -> x.^(1/5), x -> x.^15, x -> x.^(1/15))
            )

centralities = Dict();
for maps in mappings
    @show maps[1]
    x,y = compute_centrality(B,W,maps[2],maxiter=500);
    centralities[maps[1]] = Dict("x"=>x, "y"=>y)
end

include("make-top-centrality-table.jl")
top_labels_by_centrality = make_centrality_table(centralities,labels, number_of_top_nodes = 10)



ftlabels = 15
ftticks = 12
fttitle = 18

# node scatter plot =============================================================
figure()
suptitle("$name_title node centrality", fontsize=fttitle)
subplot(131)
x = centralities["linear"]["x"];
xx = centralities["log-exp"]["x"]; 
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("linear",fontsize=ftlabels)
ylabel("log-exp",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)
subplot(132)
x = centralities["max"]["x"];
xx = centralities["log-exp"]["x"]; 
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("max",fontsize=ftlabels)
ylabel("log-exp",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)
subplot(133)
x = centralities["linear"]["x"];
xx = centralities["max"]["x"]; 
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("linear",fontsize=ftlabels)
ylabel("max",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)


# edge scatter plot =============================================================
figure()
suptitle("$name_title edge centrality", fontsize=fttitle)
subplot(131)
x = centralities["linear"]["y"]; x = x[1:5000]
xx = centralities["log-exp"]["y"]; xx = xx[1:5000]
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("linear",fontsize=ftlabels)
ylabel("log-exp",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)
subplot(132)
x = centralities["max"]["y"]; x = x[1:5000]
xx = centralities["log-exp"]["y"]; xx = xx[1:5000]
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("max",fontsize=ftlabels)
ylabel("log-exp",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)
subplot(133)
x = centralities["linear"]["y"]; x = x[1:5000]
xx = centralities["max"]["y"]; xx = xx[1:5000]
plot(x./maximum(x),xx./maximum(xx), "o")
xlabel("linear",fontsize=ftlabels)
ylabel("max",fontsize=ftlabels)
xticks(fontsize=ftticks)
yticks(fontsize=ftticks)



# Scatterplot edge centrality vs edge weights ==============================
w = vec(diag(W)) #edge weights
δ = vec(sum(B,dims=1)) #edge degree

f = figure()
suptitle("$name_title edge centrality vs edge weights", fontsize=fttitle)
for (i,maps) in enumerate(keys(mappings))
    y = centralities[maps]["y"]
    subplot(1,length(keys(mappings)),i)
    plot(y./maximum(y), w./maximum(w),"o")
    if i == 1 
        ylabel("edge weight",fontsize=ftlabels)
    end
    xlabel("$maps edge centrality",fontsize=ftlabels)
    xticks(fontsize=ftticks)
    yticks(fontsize=ftticks)
end



















# scatterplot centralities ====================================
figure()
x = centralities["linear"]["x"]
y = centralities["linear"]["y"]
a = collect(keys(mappings)); list = a[a.!="linear"];
for (i,maps) in enumerate(list)
    subplot(2,length(list),i)
    xx = centralities[maps]["x"]; 
    plot(x./maximum(x),xx./maximum(xx), "o")
    xlabel("linear node centrality")
    ylabel("$maps node centrality")
    subplot(2,length(list),length(list)+i)
    yy = centralities[maps]["y"]; 
    plot(y./maximum(y),yy./maximum(yy), "o")
    xlabel("linear edge centrality")
    ylabel("$maps edge centrality")
end
tight_layout()



# scatterplot centralities ALL ====================================
figure()
for (i,maps) in enumerate(keys(mappings))
    x = centralities[maps]["x"]
    a = collect(keys(mappings)); list = a[a.!=maps];
    for (j,maps2) in enumerate(list)
        subplot(length(mappings),length(list),(i-1)*length(list)+j)
        xx = centralities[maps2]["x"]; 
        plot(x./maximum(x),xx./maximum(xx), "o")
        xlabel("$maps",fontsize=15)
        ylabel("$maps2",fontsize=15)
        xticks(fontsize=13)
        yticks(fontsize=13)
    end
end
figure()
for (i,maps) in enumerate(keys(mappings))
    y = centralities[maps]["y"]
    a = collect(keys(mappings)); list = a[a.!=maps];
    for (j,maps2) in enumerate(list)
        subplot(length(mappings),length(list),(i-1)*length(list)+j)
        yy = centralities[maps2]["y"]; 
        plot(y./maximum(y),yy./maximum(yy), "o")
        xlabel("$maps edge centrality",fontsize=15)
        ylabel("$maps2 edge centrality",fontsize=15)
        xticks(fontsize=13)
        yticks(fontsize=13)
    end
end
tight_layout()








# Kendall tau plot ==========================================
ksy = [Int64(2^i) for i in 1:floor(log2(m)) ]; append!(ksy,m)

figure()
for maps in mappings
        
    if (maps[1]=="linear") 
                    
    else
        map_name = maps[1]
        subplot(121)
        semilogx(ksx,correlations[map_name]["ke_x"], label = "$map_name", linewidth=2)
        ylabel("Kendall-τ correlation")
        title("node centrality")
        xlabel("Number top ranked elements (k)")
        subplot(122)
        semilogx(ksy,correlations[map_name]["ke_y"], label = "$map_name", linewidth=2)
        title("edge centrality")
        xlabel("Number top ranked elements (k)")
        
    end

end
legend() #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
tight_layout()
# ========================================================




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

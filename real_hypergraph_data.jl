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
    B,w,labels,names,edges = read_walmart_data();
elseif dataset_name == "stackoverflow"
    B,w,edges = read_tags_data("tags-math-sx");
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
    x,y = compute_centrality(B,maps[2],edge_weights = w, maxiter=100);
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

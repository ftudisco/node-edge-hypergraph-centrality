using LinearAlgebra
using SparseArrays
# using UnicodePlots
using PyPlot
using StatsBase

include("centrality_tools.jl")
include("sunflower_tools.jl")



number_of_nodes_per_petal = [3 3 3 3 3 3 3 3] 
number_of_nodes_per_petal = [3 4 5 6 7 8 9 10]

B, size_petals = incidence_sunflower_hypergraph(number_of_nodes_per_petal); 
n,m = size(B)

mappings = Dict("linear"    => (x -> x, x -> x, x -> x, x -> x), 
            "log-exp"       => (x -> x, x -> x.^(1/10), x -> log.(x), x -> exp.(x)),
            "max"   => (x -> x, x -> x, x -> x.^15, x -> x.^(1/15))
            )

centralities = Dict()
for maps in mappings
    @show maps[1]
    x,y = compute_centrality(B,maps[2],maxiter=500, tol=1e-6);
    centralities[maps[1]] = Dict("x"=>x, "y"=>y)
end

figure()
for (i,maps) in enumerate(keys(centralities))
    x = centralities[maps]["x"]
    subplot(1,length(keys(centralities)), i, aspect="equal")
    sizes = plot_sunflower(number_of_nodes_per_petal, sizes = x./maximum(x), mode="intvals", max_radious=maximum(number_of_nodes_per_petal))
    PyPlot.title("$maps node centrality")
    # axis("off")
    xticks([])
    yticks([])
end


figure(frameon=false); subplot(111,aspect="equal")
plot_sunflower(number_of_nodes_per_petal)
axis("off")


# compute average centralitiy over over petal sizes
h = fit(Histogram, size_petals[2:end], nbins = length(unique(number_of_nodes_per_petal)))
ind = StatsBase.binindex.(Ref(h),size_petals[2:end])

bins = unique(ind)
average_node_centrality  = Dict()
for maps in mappings 
    x = centralities[maps[1]]["x"][2:end]
    x = x./maximum(x)
    average_node_centrality[maps[1]] = zeros(size(bins))
    # @show x
    for (i,b) in enumerate(bins)
        average_node_centrality[maps[1]][i] = mean(x[ind.==b])
    end
end

figure()
for maps in mappings
    plot(bins,average_node_centrality[maps[1]], "-o")
end

    

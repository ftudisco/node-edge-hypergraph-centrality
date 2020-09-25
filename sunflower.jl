using LinearAlgebra
using SparseArrays
using UnicodePlots
using PyPlot

input("centrality_tools.jl")

# B = incidence matrix ~ n x m (nodes x edges)
function incidence_sunflower_hypergraph(r,k)
    number_of_nodes_per_petal = [k for _ in 1:r]
    
    n = sum(number_of_nodes_per_petal)
    m = length(number_of_nodes_per_petal)
    
    Is = [1 for _ in 1:m]
    Js = [i for i in 1:m]
    n = 1
    for (e,size_e) in enumerate(number_of_nodes_per_petal.-1)
        for i in 1:size_e
            push!(Is,i+n)
            push!(Js,e)
        end 
        n = Is[end]
    end
        
    return sparse(Is, Js,1), Is
end

function incidence_sunflower_hypergraph(number_of_nodes_per_petal)
    n = sum(number_of_nodes_per_petal)
    m = length(number_of_nodes_per_petal)
    
    Is = [1 for _ in 1:m]
    Js = [i for i in 1:m]
    n = 1
    for (e,size_e) in enumerate(number_of_nodes_per_petal.-1)
        for i in 1:size_e
            push!(Is,i+n)
            push!(Js,e)
        end 
        n = Is[end]
    end
        
    return sparse(Is, Js,1)
end





a = 2
b = 3
f = x -> x.^a
g = x -> x.^(1/a)
ϕ = x -> x.^b
ψ = x -> x.^(1/b)

number_of_nodes_per_petal = [3 6 9 12 15];
number_of_nodes_per_petal = [number_of_nodes_per_petal number_of_nodes_per_petal]
B = incidence_sunflower_hypergraph(number_of_nodes_per_petal);
x,y = compute_centrality(B,f,g,ϕ,ψ);
subplot(121)
PyPlot.plot(x,"o")
subplot(122)
PyPlot.plot(y,"x")

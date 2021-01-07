using LinearAlgebra
using SparseArrays
# using UnicodePlots
using PyPlot
using StatsBase

input("centrality_tools.jl")

function compute_centrality(B,(f,g,ϕ,ψ); maxiter=100,tol=1e-3 )
    n,m = size(B)
        
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    # x0 = rand(n,1)./n
    # y0 = rand(m,1)./m
    
     
    for it in 1:maxiter
        if (it<10) || (it%30==0) println("$it ..."); end
        u = g(B*f(y0))
        v = ψ(B'*ϕ(x0))
        x = u./norm(u,1)
        y = v./norm(v,1)
        
        check = norm(x-x0,1) + norm(y-y0)
        if check < tol
            println("$it ===")
            return vec(x),vec(y)
        else
            x0 = copy(x)
            y0 = copy(y)
        end
    end
    println("Warning: Centrality did not converge")
    return vec(x0),vec(y0)
end


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
    
    size_petals = []
    append!(size_petals,n)
    
    Is = [1 for _ in 1:m]
    Js = [i for i in 1:m]
    n = 1
    for (e,size_e) in enumerate(number_of_nodes_per_petal.-1)
        for i in 1:size_e
            push!(Is,i+n)
            push!(Js,e)
            append!(size_petals, size_e)
        end 
        n = Is[end]
    end
        
    return sparse(Is, Js,1), Int64.(size_petals)
end

function coordinates_sunflower(number_of_nodes_per_petal; max_radious = 10)
    m = length(number_of_nodes_per_petal) #number of petals
    θs = range(0,2*π,length=m+1)
    n = sum(number_of_nodes_per_petal) #number of nodes
    XY = [0 0];
    for (i,n_in_petal) in enumerate(number_of_nodes_per_petal)
        θ = θs[i]
        r = range(0,max_radious,length=n_in_petal)
        for j in 1:n_in_petal
            if j > 1
                XY = [XY ;r[j].*[cos(θ) sin(θ)] ]
            else 
                # if i == 1
                #     XY = [XY ;r[j].*[cos(θ) sin(θ)] ]
                # end
            end
        end
    end
    return XY[1:end,:]
end

function plot_sunflower(number_of_nodes_per_petal; max_dot_size = 15, min_dot_size = 1, max_radious = 10, sizes = "equal", mode="centrality")
    n = sum(number_of_nodes_per_petal) #number of nodes
    if sizes == "equal"
        sizes = [1 for _ in 1:n]
    end
    # sizes = linearize(sizes,0.1,1);
    sizes = sizes./maximum(sizes)
    sizes = round.(sizes, digits=5)
    intvals = zeros(size(sizes))
    u = unique(sizes); p = sortperm(u,rev=true);
    # sarray = range(max_dot_size,min_dot_size,length=length(u))
    for (i,s) in enumerate(u)
        intvals[sizes.==u[p[i]]].= max_dot_size/i
    end
    
    XY = coordinates_sunflower(number_of_nodes_per_petal, max_radious=max_radious)
    for i = 1 : size(XY)[1]
        if mode == "centrality"
            plot(XY[i,1],XY[i,2],"o",markersize=max_dot_size*sizes[i], color="#209cee")
        elseif mode == "intvals"
            plot(XY[i,1],XY[i,2],"o",markersize=intvals[i], color="#209cee")
        end
    end
    return intvals
end


function linearize(x,a,b)
    xx = x;
    alpha = (b-a)./(maximum(x)-minimum(x[x.>0]))
    beta = alpha*minimum(x[x.>0])-a;
    xx[x.>0] .= (alpha.*x[x.>0]).-beta;
    return xx
end





number_of_nodes_per_petal = [3 3 3 3 3 3 3 3] 
number_of_nodes_per_petal = [3 6 12 25 50 100]
number_of_nodes_per_petal = [3 3 3 3 10 10 10]
number_of_nodes_per_petal = [3 4 5 6 7 8 9 10]

number_of_nodes_per_petal = rand(3:6, 15)

B, size_petals = incidence_sunflower_hypergraph(number_of_nodes_per_petal); 
n,m = size(B)


mappings = Dict("linear"    => (x -> x, x -> x, x -> x, x -> x), 
            "log-exp"       => (x -> x, x -> x.^(1/10), x -> log.(x), x -> exp.(x)),
            # "max-lin"   => (x -> x.^10, x -> x.^(1/10), x -> x, x -> x),
            # "max-max"   => (x -> x.^10, x -> x.^(1/10), x -> x.^10, x -> x.^(1/10)),
            "max"   => (x -> x, x -> x, x -> x.^15, x -> x.^(1/15))
            # "lin-max"   => (x -> x, x -> x, x -> x.^(10), x -> x.^(1/10))   
            )

centralities = Dict()
for maps in mappings
    @show maps[1]
    x,y = compute_centrality(B,maps[2],maxiter=500);
    centralities[maps[1]] = Dict("x"=>x, "y"=>y)
end

figure()
for (i,maps) in enumerate(keys(centralities))
    x = centralities[maps]["x"]
    subplot(1,length(keys(centralities)), i, aspect="equal")
    sizes = plot_sunflower(number_of_nodes_per_petal, sizes = x./maximum(x), mode="centrality", max_radious=maximum(number_of_nodes_per_petal))
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

    



number_of_nodes_per_petal = [number_of_nodes_per_petal number_of_nodes_per_petal]
B = incidence_sunflower_hypergraph(number_of_nodes_per_petal);
x,y = compute_centrality(B,f,g,ϕ,ψ);
figure(frameon=false); subplot(111,aspect="equal")
sizes = plot_sunflower(number_of_nodes_per_petal, sizes = x./maximum(x), max_radious = 20)
axis("off")

subplot(121)
PyPlot.plot(x,"o")
subplot(122)
PyPlot.plot(y,"x")


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

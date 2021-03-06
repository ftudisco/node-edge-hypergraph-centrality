
using Combinatorics
using DelimitedFiles

function read_walmart_data()
    edges = Dict()
    open("data/walmart-trips/hyperedges-walmart-trips.txt") do f
        for (e,line) in enumerate(eachline(f))
            nodes = [parse(Int64, v) for v in split(line, ',')]
            edge = Set(nodes)
            if (edge in keys(edges))
                edges[edge] = edges[edge]+1
            else
                edges[edge] = 1
            end
        end
    end
    
    Is = []
    Js = []
    w = []
    for (e,edge) in enumerate(keys(edges))
        append!(Is,edge)
        append!(Js,[e for _ in edge])
        append!(w,edges[edge])
    end
    m = length(w)

    labels = Int64[]
    open("data/walmart-trips/node-labels-walmart-trips.txt") do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    
    label_name = Dict()
    open("data/walmart-trips/label-names-walmart-trips.txt") do f
        for (idx,name) in enumerate(eachline(f))
            label_name[idx] = name
        end
    end
    

    return sparse(Is,Js,1),  Float64.(w), labels, label_name, edges
end

function read_tags_data(dataset::String)
    
    read(filename::String) = convert(Vector{Int64}, readdlm(filename, Int64)[:, 1])
    simplices = read("data/$(dataset)/$(dataset)-simplices.txt")
    nverts = read("data/$(dataset)/$(dataset)-nverts.txt")
    
    edges = Dict()
    Is = []
    Js = []
    w = []
    curr_ind = 1
    
    
    for nvert in nverts
        edge = Set(simplices[curr_ind:(curr_ind + nvert - 1)])
        curr_ind += nvert
        
        if (edge in keys(edges))
            edges[edge] = edges[edge]+1
        else
            edges[edge] = 1
        end
        
    end
    
    for (e,edge) in enumerate(keys(edges))
        append!(Is,edge)
        append!(Js,[e for _ in edge])
        append!(w,edges[edge])
    end
    
    m = length(w)
    
    return sparse(Is,Js,1), Float64.(w), edges
end


function read_tags_data_node_labels(dataset::String)
    labels = String[]
    open("data/$(dataset)/$(dataset)-node-labels.txt") do f
        for line in eachline(f)
            data = split(strip(line))
            push!(labels, join(data[2:end], " "))
        end
    end
    return labels
end

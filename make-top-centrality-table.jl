
function make_centrality_table(centralities,labels; number_of_top_nodes = 10)
    
    top_labels_by_centrality = Dict()
    for (i,maps) in enumerate(keys(centralities))
        x = centralities[maps]["x"]
        xid = sortperm(vec(x),rev=true)
        xx = x[xid]
        
        top_nodes_by_centrality = xid[1:number_of_top_nodes]
        top_labels_by_centrality[maps] = labels[top_nodes_by_centrality]
    end

    cc = String(['c' for _ in 1 : length(keys(centralities))])
    
    println("""\\begin{tabular}{$cc}""")
    println(raw"""\toprule""")
    for (j,map) in enumerate(keys(centralities))
        mapname = uppercasefirst(map)
        if j < length(keys(centralities))
            print("""\\textit{$mapname} & """)
        else
            println("""\\textit{$mapname} \\\\""")
        end
    end
    println(raw"""\midrule""")    
    for i in 1:number_of_top_nodes
        for (j,map) in enumerate(keys(centralities))
            nodename = uppercasefirst(replace(top_labels_by_centrality[map][i], "-"=>" "))
            if j < length(keys(centralities))
                print("""$nodename & """)
            else
                println("""$nodename \\\\""")
            end
        end
    end
    println(raw"""\bottomrule""")
    println(raw"""\end{tabular}""")    
    
    return top_labels_by_centrality
end

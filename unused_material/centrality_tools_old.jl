using Distances


function compute_centrality(B,(f,g,ϕ,ψ); maxiter=100,tol=1e-6,edge_weights = ones(size(B)[1],1), node_weights = ones(size(B)[2],1)  )
    n,m = size(B)
        
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    # x0 = rand(n,1)./n
    # y0 = rand(m,1)./m
    
    W = Diagonal(edge_weights)
    N = Diagonal(node_weights)
     
    for it in 1:maxiter
        if (it<10) || (it%30==0) println("$it ..."); end
        u = g(B*W*f(y0))
        v = ψ(B'*N*ϕ(x0))
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

function compute_centrality(B,W,(f,g,ϕ,ψ); maxiter=100,tol=1e-6  )
    n,m = size(B)
        
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    # x0 = rand(n,1)./n
    # y0 = rand(m,1)./m

     
    for it in 1:maxiter
        if (it<10) || (it%30==0) println("$it ..."); end
        u = g(B*W*f(y0))
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

function compute_centrality(B,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
    n,m = size(B)
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    
    for it in 1:maxiter
        if (it%30==0) println("$it ..."); end
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
# function compute_centrality(B,W,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
#     n,m = size(B)
# 
#     x0 = ones(n,1)./n
#     y0 = ones(m,1)./m
# 
#     for it in 1:maxiter
#         if (it%30==0) println("$it ..."); end
#         u = g(B*W*f(y0))
#         v = ψ(B'*ϕ(x0))
#         x = u./norm(u,1)
#         y = v./norm(v,1)
# 
#         check = norm(x-x0,1) + norm(y-y0)
#         if check < tol
#             println("$it ===")
#             return vec(x),vec(y)
#         else
#             x0 = copy(x)
#             y0 = copy(y)
#         end
#     end
#     println("Warning: Centrality did not converge")
#     return vec(x0),vec(y0)
# end


function compute_centrality(B,W,N,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
    n,m = size(B)
    @assert size(B)[2] == size(W)[1]
    @assert size(B)[1] == size(N)[1]
    
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    
    for it in 1:maxiter
        if (it%30==0) println("$it ..."); end
        u = g(B*W*f(y0))
        v = ψ(B'*N*ϕ(x0))
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


# given two vectors u and v returns the sequence of correlation values between top k ranked entries of u and v
function correlation_map(u,v)
    m = length(u);
    spearman = [];
    kendall = [];
    jacc = [];
    idu = sortperm(u,rev=true)
    idv = sortperm(v,rev=true)
    
    for i in 1:floor(log2(m))
        k = Int64(2^i)
        append!(  spearman , corspearman(vec(idu[1:k]), vec(idv[1:k]))  )
        append!(  kendall , corkendall(vec(idu[1:k]), vec(idv[1:k]))  )
        append!(  jacc , Distances.jaccard(vec(idu[1:k]), vec(idv[1:k]))  )
        println("$k...")
    end
    append!(  spearman , corspearman(vec(idu), vec(idv))  )
    append!(  kendall , corkendall(vec(idu), vec(idv))  )
    append!(  jacc , Distances.jaccard(vec(idu), vec(idv))  )
    
    spearman[isnan.(spearman)].=0
    kendall[isnan.(kendall)].=0
    jacc[isnan.(jacc)].=0
    return spearman, kendall, jacc
end

function correlation_table(x,y,centrality_dict)
    x = x./maximum(x)
    y = y./maximum(y)
    correlations = Dict()
    for key in keys(centrality_dict)
        xx = centrality_dict[key]["x"]
        xx = xx./maximum(xx)
        @show length(x), length(xx)
        sp_x,ke_x,ja_x = correlation_map(x,xx)
        # m_x = mxarray(Float64.(vec(x)))
        # m_xx = mxarray(Float64.(vec(xx)))
        # isim_x = mat"isim_new(m_x, m_xx)"
        yy = centrality_dict[key]["y"]
        yy = yy./maximum(yy)
        sp_y,ke_y,ja_y = correlation_map(y,yy)
        # m_y = mxarray(Float64.(vec(y)))
        # m_yy = mxarray(Float64.(vec(yy)))
        # isim_y = mat"isim_new(m_y, m_yy)"
        println("$key")
        correlations[key] = Dict("sp_x"=>sp_x, "sp_y"=>sp_y, "ke_x"=>ke_x, "ke_y"=>ke_y, "ja_x"=>ja_x, "ja_y"=>ja_y)#, "isim_x"=>isim_x, "isim_y"=>isim_y)
    end
    return correlations
end

function centrality_table(αs,βs,B,W;maxiter=500)
    centralities = Dict()
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
    return centralities
end

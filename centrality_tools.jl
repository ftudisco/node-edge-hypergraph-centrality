using SparseArrays


function compute_centrality(B,(f,g,ϕ,ψ); 
                maxiter=100,
                tol=1e-6,
                edge_weights = ones(size(B)[2],1), 
                node_weights = ones(size(B)[1],1),
                mynorm = (x -> norm(x,1))  )
                
    B = sparse(B)
    n,m = size(B)
        
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    # x0 = rand(n,1)./n
    # y0 = rand(m,1)./m
    
    W = spdiagm(0=>edge_weights[:])
    N = spdiagm(0=>node_weights[:])
     
    for it in 1:maxiter
        if (it<10) || (it%30==0) println("$it ..."); end
        u = sqrt.( x0 .* g(B*W*f(y0)) )
        v = sqrt.( y0 .* ψ(B'*N*ϕ(x0)) )
        x = u./mynorm(u)
        y = v./mynorm(v)
        
        check = mynorm(x-x0) + mynorm(y-y0)
        if check < tol
            println("$it ===")
            return vec(x),vec(y)
        else
            x0 = copy(x)
            y0 = copy(y)
        end
    end
    println("Warning: Centrality did not reach tol = $tol in $maxiter iterations\n-------  Relative error reached so far = $check")
    return vec(x0),vec(y0)
end

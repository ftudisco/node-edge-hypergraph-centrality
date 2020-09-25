
function compute_centrality(B,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
    n,m = size(B)
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    
    for it in 1:maxiter
        println("$it ...")
        u = g(B*f(y0))
        v = ψ(B'*ϕ(x0))
        x = u./norm(u,1)
        y = v./norm(v,1)
        
        check = norm(x-x0,1) + norm(y-y0)
        if check < tol
            return vec(x),vec(y)
        else
            x0 = copy(x)
            y0 = copy(y)
        end
    end
    println("Warning: Centrality did not converge")
    return vec(x),vec(y)
end
function compute_centrality(B,W,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
    n,m = size(B)
    @assert size(B)[2] == size(W)[1]
        
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    
    for it in 1:maxiter
        println("$it ...")
        u = g(B*W*f(y0))
        v = ψ(B'*ϕ(x0))
        x = u./norm(u,1)
        y = v./norm(v,1)
        
        check = norm(x-x0,1) + norm(y-y0)
        if check < tol
            return vec(x),vec(y)
        else
            x0 = copy(x)
            y0 = copy(y)
        end
    end
    println("Warning: Centrality did not converge")
    return vec(x),vec(y)
end
function compute_centrality(B,W,N,f,g,ϕ,ψ; maxiter=100,tol=1e-6)
    n,m = size(B)
    @assert size(B)[2] == size(W)[1]
    @assert size(B)[1] == size(N)[1]
    
    x0 = ones(n,1)./n
    y0 = ones(m,1)./m
    
    for it in 1:maxiter
        println("$it ...")
        u = g(B*W*f(y0))
        v = ψ(B'*N*ϕ(x0))
        x = u./norm(u,1)
        y = v./norm(v,1)
        
        check = norm(x-x0,1) + norm(y-y0)
        if check < tol
            return vec(x),vec(y)
        else
            x0 = copy(x)
            y0 = copy(y)
        end
    end
    println("Warning: Centrality did not converge")
    return vec(x),vec(y)
end    

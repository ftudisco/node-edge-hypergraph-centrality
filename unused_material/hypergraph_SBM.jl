

using Distributions


p = 0.3 # prob to select a node from the same set
q = 0.1 # prob to select a node from a different set

blocks = [50 50 50]
n = sum(blocks)
m = 200
k = length(blocks)
g = .4 # prob to have an edge
dist = Geometric(.4)
edge_lengths = rand(dist, m)

S = [];
for (i,k) in enumerate(blocks)
    push!(S,Set(Int64.((i-1)*k+1:i*k)))
end

for ed in edge_lengths
    block = rand(1:k)
    for i = 1 : ed
        
end

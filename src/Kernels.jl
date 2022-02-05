function undef_copy(v::Array{T, d}) where {T, d}
    Array{T, d}(undef, size(v))
end

function make_volterra_kernel(K::Array{T, 2}) where {T <: Number}
    n = size(K, 1)
    
    V = LowerTriangular(Array{T, 2}(undef, (n, n)))
    
    V[1,1] = 0
    
    for i  = 2:n
        V[i,1] = K[i,1]
        V[i,i] = K[i,i]
        
        for j = 2:(i-1)
            V[i,j] = 2 * K[i,j]
        end
    end
    
    return(V)
end

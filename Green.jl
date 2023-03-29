function green_matrix(N::Integer, Rij::Array; k=1)
    G = Array{Complex{Float64}}(undef, N, N)
    for i in eachindex(G)
        @inbounds G[i] = cis(k*Rij[i])/(1im*k*Rij[i])
    end
    # Acessa e preenche a diagonal principal da matriz de Green
    G[LinearAlgebra.diagind(G)].= 1

    return G

end

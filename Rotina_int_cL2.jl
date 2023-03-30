#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------ Estatisticas da Intensidade do Laser Gaussiano -------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

function Rotina_int_cL2(
    N_rep::Integer,
    N::Number,
    Γ₀::Integer,
    k::Number,
    Δ::Number,
    λ::Number,
    h::Number,
    w₀::Number,
    fator_r::Number,
    Distancia_do_sistema::Number,
    N_sensores::Integer,
    Angulo_de_variação_inicial_1::Integer,
    Angulo_de_variação_final_2::Integer,
    Radius::Float64,
)
    s=1e-5
    E₀= √(s * (1 + 4(Δ / Γ₀)^2)/2)
    ρ = (N / ((pi * (Radius^2))*h))
    rₘ = ρ^(-1 / 3) / 5
    # b₀ = (6 * N) / ((k * Radius)^2)
    p = Progress(length(N_rep); showspeed = true)
    betas = pmap(1:N_rep) do rep
        # r_instantaneo = getSphere(Radius, N, rₘ)
        r_instantaneo = getCilindro(Radius, N, rₘ,h)

        Rij_instantaneo = Distances.pairwise(Euclidean(), r_instantaneo, r_instantaneo, dims = 1)

        G = -(Γ₀/2)*green_matrix(N, Rij_instantaneo)
        G[diagind(G)] .= im*Δ -(Γ₀/2)

        Ωₙ = (-im/2).*get_campo_over_atoms(r_instantaneo, w₀, k, E₀) # feixe gaussiano
        βₛ = -(G\Ωₙ)

        βₛ, r_instantaneo
    end



    todas_matrizes = Array{Complex{Float64}}(undef, N, N_rep)
    r_totais = zeros(N, 3, N_rep)

    for rep in 1:N_rep
        todas_matrizes[:, rep] .= betas[rep][1]
        r_totais[:, :, rep] .= betas[rep][2]
    end
    betas_totais = todas_matrizes



    I_instantanea_total = zeros(((N_sensores ÷ 2)-1), N_rep)
    Intensidade = zeros(((N_sensores ÷ 2)-1), N_rep)


    for rep in 1:N_rep
        Sensores = get_tela(N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, fator_r, Radius, Distancia_do_sistema)



        I_instantanea_total[:, rep] = get_detector_intensities_cl(Sensores, r_totais[:, :, rep], betas_totais[:, rep], E₀, w₀, k)
        Intensidade[:, rep] = I_instantanea_total[:, rep] ./ mean(I_instantanea_total[:, rep])

    end

    return Intensidade

end

function calcularParesIndices(x, Δi)
    totalIndices = length(x)
    indices = 1:totalIndices
    y = mod.(indices .+ Δi, totalIndices)
    y[findall(y .== 0)] .= totalIndices
    paresIndices = [indices y]
    return paresIndices
end

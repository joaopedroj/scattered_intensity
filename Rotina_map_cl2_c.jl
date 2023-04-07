#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------ Estatisticas da Intensidade do Laser Gaussiano -------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

function Rotina_g2_paralelo22(
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
    t_correlacao::Any,
)
    s=1e-5
    E₀= √(s * (1 + 4(Δ / Γ₀)^2)/2)
    ρ = (N / ((pi * (Radius^2))*h))
    rₘ = ρ^(-1 / 3) / 5
    # b₀ = (6 * N) / ((k * Radius)^2)
    p = Progress(length(N_rep); showspeed = true)
    betas = pmap(1:N_rep) do rep


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

    g2_parcial = zeros(t_correlacao, N_rep)
    g2_final = zeros(size(t_correlacao, 1))
    g2 = zeros((N_sensores ÷ 2))
    g22 = zeros((N_sensores ÷ 2))


    # matriz_r_final = zeros(((N_sensores ÷ 2)-1), 3)


    for rep in 1:N_rep
        Sensores = get_tela(N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, fator_r, Radius, Distancia_do_sistema)


        I_instantanea_total = zeros(((N_sensores ÷ 2)-1), N_rep)
        Intensidade = zeros(((N_sensores ÷ 2)-1), N_rep)
        

        I_instantanea_total[:, rep] = get_detector_intensities_cl(Sensores, r_totais[:, :, rep], betas_totais[:, rep], E₀, w₀, k)
        Intensidade[:, rep] = I_instantanea_total[:, rep] ./ mean(I_instantanea_total[:, rep])

        # matriz_r_final = Intensidades_esfericas(Sensores, Intensidade[:, rep])

        # # ang_interesse = 59.5
        # # ang_interesse2 = 60.5
        # # intensidade = matriz_r_final[findall(ang_interesse2 .>= matriz_r_final[:, 3] .>= ang_interesse), 1]
        # intensidade = matriz_r_final[:, 1]
        intensidade = Intensidade
        g2 = zeros(length(intensidade), 1)
        g22 = zeros(length(intensidade), 1)
        parestotais = Array{Float64}(undef, length(intensidade), 2, length(t_correlacao))

        for t in t_correlacao
            parestotais[:, :, t] .= calcularParesIndices(intensidade, t)
            for i in 1:length(intensidade)
                t₁ = parestotais[i, 1, t]
                t₂ = parestotais[i, 2, t]
                autocorrelação = intensidade[round(Int64, t₁)] * intensidade[round(Int64, t₂)]
                autocorrelação2 = intensidade[round(Int64, t₁)]
                g2[i] = autocorrelação
                g22[i] = autocorrelação2
            end
            g2_parcial[t, rep] = mean(g2) / mean(g22)^2
        end
    end
    for t in t_correlacao
        g2_final[t] = mean(g2_parcial[t, :])

    end

    return g2_final

end

function calcularParesIndices(x, Δi)
    totalIndices = length(x)
    indices = 1:totalIndices
    y = mod.(indices .+ Δi, totalIndices)
    y[findall(y .== 0)] .= totalIndices
    paresIndices = [indices y]
    return paresIndices
end

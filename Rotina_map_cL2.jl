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
    w₀::Number,
    Distancia_do_sistema::Number,
    N_sensores::Integer,
    Angulo_de_variação_inicial_1::Integer,
    Angulo_de_variação_final_2::Integer,
    Radius::Float64,
    th::Number
)
    s = 1e-6
    # E₀ = √(s * (1 + 4(Δ / Γ₀)^2) / 2)
    E₀ = Γ₀ * √(s / 2)
    ρ = (N / ((4 / 3) * pi * (Radius^3)))
    rₘ = ρ^(-1 / 3) / 5
    # b₀ = (6 * N) / ((k * Radius)^2)

    p = Progress(length(N_rep); showspeed=true)
    betas = pmap(1:N_rep) do rep
        r_instantaneo = getSphere(Radius, N, rₘ)

        Rij_instantaneo = Distances.pairwise(Euclidean(), r_instantaneo, r_instantaneo, dims=1)

        G = -(Γ₀ / 2) * green_matrix(N, Rij_instantaneo)
        G[diagind(G)] .= im * Δ - (Γ₀ / 2)

        Ωₙ = (-im / 2) .* get_campo_over_atoms(r_instantaneo, w₀, k, E₀) # feixe gaussiano
        βₛ = -(G \ Ωₙ)


        βₛ, r_instantaneo
    end




    todas_matrizes = Array{Complex{Float64}}(undef, N, N_rep)
    r_totais = zeros(N, 3, N_rep)

    for rep in 1:N_rep
        todas_matrizes[:, rep] .= betas[rep][1]
        r_totais[:, :, rep] .= betas[rep][2]
    end
    betas_totais = todas_matrizes

    # I_instantanea_total = zeros(((N_sensores ÷ 2) - 1), N_rep)
    I_instantanea_total = zeros(N_sensores, N_rep)
    # Intensidade = zeros(((N_sensores ÷ 2) - 1), N_rep)
    # matriz_r_final = zeros((N_sensores), 3)

    for rep in 1:N_rep
        # Sensores = get_tela(N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, Distancia_do_sistema)
        # Sensores = get_sensors_ring(num_pts=100, φ_inicial=0, φ_final=2π, kR=1, θ=π / 2)
        Sensores = get_sensors_ring(num_pts=100, φ_inicial=0, φ_final=2π, kR=1, θ=th)


        I_instantanea_total[:, rep] = get_detector_intensities_cl(Sensores, r_totais[:, :, rep], betas_totais[:, rep], E₀, w₀, k)



    end
    I_instantanea_total = reduce(vcat, I_instantanea_total)
    int_n = I_instantanea_total ./ mean(I_instantanea_total)
    # variancia = (mean(int .^ 2) / (mean(int)^2)) - 1
    variancia = var(int_n)

    return variancia
end

function calcularParesIndices(x, Δi)
    totalIndices = length(x)
    indices = 1:totalIndices
    y = mod.(indices .+ (Δi - 1), totalIndices)
    y[findall(y .== 0)] .= totalIndices
    paresIndices = [indices y]
    return paresIndices
end

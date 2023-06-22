
function get_campo_over_atoms(r::Array, w₀::Number, k::Number, E₀::Number)
    N_pontos = size(r, 1)
    E = ComplexF64[]
    campos = Array{Complex{Float64}}(undef, N_pontos)

    i = 1
    for row in eachrow(r)
        x = row[1]
        y = row[2]
        z = row[3]

        ## This formula is stable for z==0
        # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
        denominator_factor = 1 .+ 2im .* z / (k * w₀^2)
        Eᵢ = E₀ .* exp.(+im * k * z)
        Eᵢ = Eᵢ .* exp.(-(x .^ 2 + y .^ 2) ./ (denominator_factor .* w₀^2))
        Eᵢ = Eᵢ ./ denominator_factor
        campos[i] = Eᵢ
        i += 1
    end

    return campos
end

function get_eletric_field_one_position(position::Array, w₀::Number, k::Number, E₀::Number)
    x = position[1]
    y = position[2]
    z = position[3]

    ## This formula is stable for z==0
    # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
    denominator_factor = 1 .+ 2im .* z / (k * w₀^2)
    Eᵢ = E₀ .* exp.(+im * k * z)
    Eᵢ = Eᵢ .* exp.(-(x .^ 2 + y .^ 2) ./ (denominator_factor .* w₀^2))
    Eᵢ = Eᵢ ./ denominator_factor
    return Eᵢ
end


function get_field_over_one_point(detector::Array, atoms::Array, β::Array, E₀::Number, k::Number, w₀::Number)
    ## Laser Pump
    E_L = im * get_eletric_field_one_position(detector, w₀, k, E₀)
    ## Scattered
    N = length(β)
    E_scatt = zero(ComplexF64)
    n_position = detector ./ norm(detector)
    dot_n_r = zero(ComplexF64)

    for j in 1:N
        r_nx = atoms[j, 1]
        r_ny = atoms[j, 2]
        r_nz = atoms[j, 3]
        # dot_n_r = cis(- (n_position[1] * r_nx + n_position[2] * r_ny + n_position[3] * r_nz))
        # E_scatt += dot_n_r * β[j]
        ikr = im * k * norm(detector - atoms[j, :])
        E_scatt = E_scatt + (β[j] * exp(ikr) / ikr)
    end
    return abs2(E_L + E_scatt)
end

function get_detector_intensities(all_detectors::Array, r::Array, β::Array, Eₒ::Number, w₀::Number, k::Number)
    n_detectors = size(all_detectors, 1)
    intensities = zeros(n_detectors)
    # Threads.@threads 
    for i in 1:n_detectors
        # one_detector = all_detectors[i, :]
        intensities[i] = get_field_over_one_point(all_detectors[i, :], r, β, Eₒ, k, w₀)
    end
    return intensities
end

# function get_field_over_one_point_cl(detector::Array, atoms::Array, β::Array, E₀::Number, k::Number, w₀::Number)
#     ## Laser Pump
#     E_L = im*get_eletric_field_one_position(detector, w₀, k, E₀)
#     ## Scattered
#     N = length(β)
#     E_scatt = zero(ComplexF64)
#     n_position = detector/norm(detector)
#     for j=1:N
#         r_nx = atoms[j,1]
#         r_ny = atoms[j,2]
#         r_nz = atoms[j,3]
#         dot_n_r = exp(-im*( n_position[1]*r_nx
#         + n_position[2]*r_ny
#         + n_position[3]*r_nz) )
#         E_scatt += dot_n_r*β[j]
#     end
#     ikr = im*k*norm(detector)
#     E_scatt = E_scatt*exp(ikr)/ikr
#     # return abs2(E_scatt)
#     return abs2(E_L + E_scatt)
# end

function get_field_over_one_point_cl(detector::Array, atoms::Array, β::Array, E₀::Number, k::Number, w₀::Number)
    ## Laser Pump
    E_L = (-im / 2) * get_eletric_field_one_position(detector, w₀, k, E₀)

    ## Scattered
    N = length(β)
    E_scatt = zero(ComplexF64)
    n_position = detector / norm(detector)
    for j = 1:N
        r_nx = atoms[j, 1]
        r_ny = atoms[j, 2]
        r_nz = atoms[j, 3]
        dot_n_r = exp(-im * (n_position[1] * r_nx + n_position[2] * r_ny + n_position[3] * r_nz))
        E_scatt += dot_n_r * β[j]
    end
    ikr = im * k * norm(detector)
    E_scatt = (-0.5) * E_scatt * exp(ikr) / ikr #'Gamma=1' --> 'Gamma/2' == '0.5'
    return abs2(E_L + E_scatt)
end

function get_detector_intensities_cl(all_detectors::Array, r::Array, β::Array, Eₒ::Number, w₀::Number, k::Number)
    n_detectors = size(all_detectors, 1)
    intensities = zeros(n_detectors)
    # Threads.@threads 
    for i in 1:n_detectors
        intensities[i] = get_field_over_one_point_cl(all_detectors[i, :], r, β, Eₒ, k, w₀)
    end
    return intensities
end




function get_field_over_one_point_cl_EL(detector::Array, E₀::Number, k::Number, w₀::Number)
    ## Laser Pump
    E_L = (-im / 2) * get_eletric_field_one_position(detector, w₀, k, E₀)


    return abs2(E_L)
end
function get_detector_intensities_cl_EL(all_detectors::Array, Eₒ::Number, w₀::Number, k::Number)
    n_detectors = size(all_detectors, 1)
    intensities = zeros(n_detectors)
    # Threads.@threads 
    for i in 1:n_detectors
        # one_detector = all_detectors[i, :]
        intensities[i] = get_field_over_one_point_cl_EL(all_detectors[i, :], Eₒ, k, w₀)
    end
    return intensities
end






# 

function get_degree_azimut(x, y)

    if x >= 0 && y >= 0
        θ = 90 - atand(x / y)
    elseif x <= 0 && y >= 0
        θ = 90 - atand(x / y)
    elseif x <= 0 && y <= 0
        θ = 270 - atand(x / y)
    elseif x >= 0 && y <= 0
        θ = 270 - atand(x / y)
    end

    return θ
end


# function Intensidade_media_tempo(Sensores::Array, r_totais::Array, betas_totais::Array, E₀::Float64, w₀::Float64, k::Number, nTempos::Integer)
#     I_instantanea_total = zeros(N_sensores, nTempos)

#     for i in 1:nTempos   # evolução Evolucao_temporal
#         I_instantanea_media = zeros(N_sensores, N_rep)

#         for j in 1:N_rep
#             # β_instantaneo = betas_totais[:,i,j]
#             I_instantanea = get_detector_intensities(Sensores, r_totais[:, :, j], betas_totais[:, i, j], E₀, w₀, k)
#             I_instantanea_media[:, j] = I_instantanea
#         end
#         I_instantanea_media_pura = zeros(N_sensores, 1)
#         for k in 1:N_sensores
#             I_instantanea_media_pura[k, 1] = mean(I_instantanea_media[k, :])
#         end
#         I_instantanea_total[:, i] = I_instantanea_media_pura
#     end

#     return I_instantanea_total ./ mean(I_instantanea_total)
# end

# function Intensidade_por_rep(Sensores::Array, r_totais::Array, betas_totais::Array, E₀::Float64, w₀::Float64, k::Number, nTempos::Integer)
#     I_instantanea_total = zeros(N_sensores, nTempos, N_rep)
#     for j in 1:N_rep
#         for i in 1:nTempos
#             I_instantanea = get_detector_intensities(Sensores, r_totais[:, :, j], betas_totais[:, i, j], E₀, w₀, k)
#             I_instantanea_total[:, i, j] = I_instantanea
#         end
#     end
#     return I_instantanea_total ./ mean(I_instantanea_total)
# end

@views function Intensidades_esfericas(Sensores::Array, Intensidade::Array)
    N_sensores = size(Sensores, 1)
    matriz_r = Array{Float64}(undef, N_sensores, 3)
    for i in 1:N_sensores
        matriz_r[i, 1] = Intensidade[i]
        matriz_r[i, 2] = get_degree_azimut(Sensores[i, 1], Sensores[i, 2])
        matriz_r[i, 3] = get_degree_azimut(Sensores[i, 3], sqrt(Sensores[i, 1]^2 + Sensores[i, 2]^2))
    end
    return matriz_r
end

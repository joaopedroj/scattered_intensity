#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------ Estatisticas da Intensidade do Laser Gaussiano -------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

function Rotina1(
    R::Number,
    N::Vector,
    Δ::Number,
    s::Number,
    w₀::Number,
    th::Any)
    ### ------------ SIMULATION SPECS ---------------------
    sensors = get_sensors_ring(; num_pts=30, kR=300, θ=th)
    maxRep = 15

    ### -------- PRODUCE INTENSITIES -----------------
    all_intensities = map(N) do N
        many_intensities = @showprogress pmap(1:maxRep) do rep
            Random.seed!(1134 + rep)

            atoms = Atom(CoupledDipoles.Sphere(), N, R)
            laser = Laser(Gaussian3D(w₀), s, Δ)
            simulation = LinearOptics(Scalar(), atoms, laser)

            βₙ = steady_state(simulation)
            intensities = scattered_intensity(simulation, βₙ, sensors; regime=:far_field)

            intensities
        end

        many_intensities = reduce(vcat, many_intensities)
        all_intensities_over_mean = many_intensities ./ mean(many_intensities)

        all_intensities_over_mean
    end

    variancia = zeros(length(N))
    for i in 1:length(N)
        variancia[i] = var(all_intensities[i])
    end

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

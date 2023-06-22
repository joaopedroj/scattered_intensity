# using Distributed, MPIClusterManagers
# @time manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
# @time addprocs(5, exeflags=`--threads 4`)

using Distributed
@time addprocs(5, exeflags=`--threads 6`)

using Plots
using FileIO
using JLD2
using Dates

@everywhere begin
    using Distances
    using Plots
    using LinearAlgebra
    using Statistics: mean, var
    using StatsBase
    using Plots
    using DifferentialEquations
    using Random
    using SpecialFunctions
    using QuadGK
    using LsqFit
    using Optim
    using ProgressMeter

    include("Pontos.jl")
    include("Green.jl")
    include("Intensidade.jl")
    include("Sensors_positions.jl")
    include("Evolucao_temporal.jl")
    include("Rotina_map_cL.jl")
    include("Rotina_map_cL2.jl")
    include("doZero_atoms.jl")


    Γ₀ = 1
    k = 1
    λ = (2 * π) / k

    N_rep = 75
    N_sensores = 100

    Angulo_de_variação_inicial_1 = 59
    Angulo_de_variação_final_2 = 61




    Radius = 4 * λ
    # w₀ = 2 * λ
    t_w = [Radius * 0.5, Radius * 0.8]
    Distancia_do_sistema = 2 * λ



    # ρ_range = range(5, 30; length=30) ./ λ^3
    # # ./ λ^3
    # t_N = zeros(length(ρ_range))
    # @. t_N = ρ_range * ((4 / 3) * pi * (Radius^3))

    ρ_range = 27 / λ^3
    t_N = ρ_range * ((4 / 3) * pi * (Radius^3))

    th = range(1, 180; length=(180 * 2)) .* (pi / 180)

    # t_correlacao = 1:1
    # # t_d = -0.5:0.1:2.1
    # t_d = range(-1.5, 1.5; length=30)
    t_d = 0.5
    # length(t_d)

end

@time Dados_var = Rotina_g2_paralelo23(N_rep, t_N, Γ₀, k, t_d, λ, t_w, Distancia_do_sistema, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, th)


darOI() = println("Começou o código")
darOI()

println("workers(): $(workers())")




save("map_cl.jld2", Dict("Saida" => Dados_var))
# save("DATA={$ano-$mes-$dia}_r.jld2", Dict("Saida" => Dados_r))
println(" - ")



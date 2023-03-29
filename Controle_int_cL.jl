# using Distributed, MPIClusterManagers
# @time manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
# @time addprocs(5, exeflags=`--threads 4`)

using Distributed
@time addprocs(5, exeflags=`--threads 4`)

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

    include("Green.jl")
    include("Intensidade.jl")
    include("Sensors_positions.jl")
    include("Rotina_int_cL2.jl")
    # include("Rotina_int_cL.jl")
    include("doZero_atoms.jl")


    Γ₀ = 1  
    k = 1
    λ = (2*π)/k


    N_rep = 100
    N_sensores = 36000*2
    angle_sensores = 5 * π / 12

    Angulo_de_variação_inicial_1 = 59
    Angulo_de_variação_final_2 = 61


    

    h=6*λ
    w₀ = 2.5*λ
    Distancia_do_sistema = 0
    # Distancia_do_sistema = h+(2*λ)

    Radius=3*λ
    fator_r = 300*λ
    # fator_r = h+(2*λ)
    # fator_r = -(Radius/2)







    # t_N = 1800:200:7800
    t_N = 6800

    # length(t_N)
    # ρ=zeros(length(t_N))
    # for i in 1:length(t_N)
    #     ρ[i] = (t_N[i] /((pi * (Radius^2))*h))
    # end
    # ρ[26]*(λ^3)
    # ρ.*(λ^3)

    t_d = 0.75
    
end

# @time Dados_var = Rotina_g2_paralelo23(N_rep, Γ₀, k, E₀, Δ, λ, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, t_N, Radius, t_correlacao);
@time Dados_var = Rotina_int_cL2(N_rep,t_N,Γ₀,k,t_d,λ,h,w₀,fator_r,Distancia_do_sistema,N_sensores,Angulo_de_variação_inicial_1,Angulo_de_variação_final_2,Radius)


darOI() = println("Começou o código")
darOI()

println( "workers(): $(workers())" )



DATA = Dates.now()
ano = Dates.year(DATA)
mes = Dates.monthname(DATA)
dia = Dates.day(DATA)

save("int_cL_c.jld2", Dict("Saida" => Dados_var))
# save("DATA={$ano-$mes-$dia}_r.jld2", Dict("Saida" => Dados_r))
println(" - ")


using Distributed
@time addprocs(5, exeflags=`--threads 6`)
# import Pkg
# Pkg.add(url="https://github.com/NoelAraujo/CoupledDipoles.jl")
# Pkg.add("ProgressMeter")
# Pkg.add("Random")
# Pkg.add("CairoMakie")

@everywhere begin
    using CairoMakie, LinearAlgebra
    using Statistics: mean, var
    using CoupledDipoles
    using ProgressMeter, Random

    include("Rotina3.jl")
    include("Rotina4.jl")

    ### ------------ ATOMS SPECS ---------------------
    k = 1
    λ = (2 * π) / k
    R = 4 * λ
    ρ_range = [0.01 / λ^3, 28 / λ^3]
    N = [Int(round(ρ_range[1] * ((4 / 3) * pi * (R^3)))), Int(round(ρ_range[2] * ((4 / 3) * pi * (R^3))))]
    # N = [20, 21]
    ### ------------ LASER SPECS ---------------------
    Δ = 1.0
    s = 1e-6
    w₀ = R / 4

    ### ------------ SENSORS SPECS ---------------------
    th = range(0, pi; length=90)
    # th = zeros(10)
    # for f in eachindex(range(0, 180; length=2))
    #     th[f] = th2[f]
    # end


end

@time Dados_var = Rotina2(R, N, Δ, s, w₀, th)


darOI() = println("Começou o código")
darOI()

println("workers(): $(workers())")




save("variancia_2.jld2", Dict("Saida" => Dados_var))
# save("DATA={$ano-$mes-$dia}_r.jld2", Dict("Saida" => Dados_r))
println(" - ")



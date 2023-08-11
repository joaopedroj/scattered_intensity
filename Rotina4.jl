function Rotina2(
    R::Number,
    N::Vector,
    Δ::Number,
    s::Number,
    w₀::Number,
    th::Any)


    var_final = @showprogress pmap(1:length(th)) do ind
        var = Rotina1(R, N, Δ, s, w₀, th[ind])
        var
    end
    return var_final

end
# function Rotina_g2_paralelo23(N_rep, t_N, Γ₀, k, t_d, λ, w₀, Distancia_do_sistema, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, t_correlacao)

#     v_m = zeros(1, length(t_N), length(t_d))
#     # g2f=zeros(size(t_correlacao, 1), length(t_N), length(t_d))
#     for q in 1:(length(t_N))
#         for d in 1:(length(t_d))
#             v_m[:, q, d] .= Rotina_g2_paralelo22(N_rep, Int(round(t_N[q])), Γ₀, k, t_d[d], λ, w₀, Distancia_do_sistema, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, t_correlacao)

#             # g2f[:,q,d]=Rotina_g2_paralelo22(N_rep, t_N[q], Γ₀, k, E₀, t_d[d], λ, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, t_correlacao)
#         end
#     end

#     return v_m
# end

function Rotina_g2_paralelo23(N_rep, t_N, Γ₀, k, t_d, λ, t_w, Distancia_do_sistema, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, th)

    v_m = zeros(1, length(t_w), length(th))
    # g2f=zeros(size(t_correlacao, 1), length(t_N), length(t_d))
    for w in 1:(length(t_w))
        for t in 1:(length(th))
            v_m[:, w, t] .= Rotina_g2_paralelo22(N_rep, Int(round(t_N)), Γ₀, k, t_d, λ, t_w[w], Distancia_do_sistema, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, th[t])

            # g2f[:,q,d]=Rotina_g2_paralelo22(N_rep, t_N[q], Γ₀, k, E₀, t_d[d], λ, N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Radius, t_correlacao)
        end
    end

    return v_m
end
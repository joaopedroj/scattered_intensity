
function sensors_positions(num_pts::Integer, R::Number, angle::Number)
    ϕ = range(deg2rad(0), deg2rad(360), length=num_pts)
    θ = ones(length(ϕ)).*angle
    x, y, z = R.*cos.(ϕ).*sin.(θ), R.*sin.(ϕ).*sin.(θ), R.*cos.(θ)
    return [x y z]
end


function get_tela(N_sensores::Integer, Angulo_de_variação_inicial_1::Integer, Angulo_de_variação_final_2::Integer, fator_r::Number, Radius::Number, Distancia_do_sistema::Number)


    Sensores = zeros(N_sensores,3)
    α_angulo_vertical = range(Angulo_de_variação_inicial_1+1,Angulo_de_variação_final_2, step=1)

    # α_angulo_vertical = (Angulo_de_variação_final_2 - Angulo_de_variação_inicial_1)
    β_angulo_horizontal = (size(α_angulo_vertical,1)*360)/(N_sensores-1)
    contador_sensores = 0
    Posição_Angular = 0:β_angulo_horizontal:(size(α_angulo_vertical,1)*180)


    for i in α_angulo_vertical.-Angulo_de_variação_inicial_1
        angulos_por_cota = Posição_Angular[findall( (Posição_Angular.>=(i-1)*360).*(Posição_Angular.<=i*360))]
        sensores_na_cota = size(findall( (Posição_Angular.>=(i-1)*360).*(Posição_Angular.<=i*360)),1)
        for j in 1:sensores_na_cota
            Sensores[j+contador_sensores,1] = (Distancia_do_sistema + Radius)*sind(Angulo_de_variação_inicial_1+i)*cosd(angulos_por_cota[j])
            Sensores[j+contador_sensores,2] = (Distancia_do_sistema + Radius)*sind(Angulo_de_variação_inicial_1+i)*sind(angulos_por_cota[j])
            Sensores[j+contador_sensores,3] = (Distancia_do_sistema + Radius)*cosd(Angulo_de_variação_inicial_1+i))+(fator_r)
        end
        contador_sensores += sensores_na_cota
    end
    Sensores=Sensores[findall(Sensores[:,2] .!= 0),:]
    return Sensores
end
function get_tela_c(N_sensores::Integer, Angulo_de_variação_inicial_1::Integer, Angulo_de_variação_final_2::Integer, fator_r::Number, Radius::Number, Distancia_do_sistema::Number)


    Sensores = zeros(N_sensores,3)
    α_angulo_vertical = range(Angulo_de_variação_inicial_1+1,Angulo_de_variação_final_2, step=1)

    # α_angulo_vertical = (Angulo_de_variação_final_2 - Angulo_de_variação_inicial_1)
    β_angulo_horizontal = (size(α_angulo_vertical,1)*360)/(N_sensores-1)
    contador_sensores = 0
    Posição_Angular = 0:β_angulo_horizontal:(size(α_angulo_vertical,1)*180)


    for i in α_angulo_vertical.-Angulo_de_variação_inicial_1
        angulos_por_cota = Posição_Angular[findall( (Posição_Angular.>=(i-1)*360).*(Posição_Angular.<=i*360))]
        sensores_na_cota = size(findall( (Posição_Angular.>=(i-1)*360).*(Posição_Angular.<=i*360)),1)
        for j in 1:sensores_na_cota
            Sensores[j+contador_sensores,1] = (fator_r + Radius)*sind(Angulo_de_variação_inicial_1+i)*cosd(angulos_por_cota[j])
            Sensores[j+contador_sensores,2] = (fator_r + Radius)*sind(Angulo_de_variação_inicial_1+i)*sind(angulos_por_cota[j])
            Sensores[j+contador_sensores,3] = (Distancia_do_sistema)
        end
        contador_sensores += sensores_na_cota
    end
    Sensores=Sensores[findall(Sensores[:,2] .!= 0),:]
    return Sensores
end




function Tela_sensores(N_sensores::Integer, Angulo_de_variação_inicial_1::Integer, Angulo_de_variação_final_2::Integer, Distancia_do_sistema::Number, Radius::Number, Coordenadas::String)

    Sensores_cartesianos = get_tela(N_sensores, Angulo_de_variação_inicial_1, Angulo_de_variação_final_2, Distancia_do_sistema, Radius)
    Sensores = Array{Float64}(undef, N_sensores, 3)

    if Coordenadas == "Esfericas"

        for k in 1:N_sensores

            x1 = Sensores_cartesianos[k,1]
            y1 = Sensores_cartesianos[k,2]
            z1 = Sensores_cartesianos[k,3]

            Sensores[k,1] = atan((((x1^2) + (y1^2))^(1/2))/z1)
            Sensores[k,2] = atan(y1/x1)
            Sensores[k,3] = sqrt((x1^2)+(y1^2)+(z1^2))

        end

    elseif Coordenadas == "Cilindricas" 

        for k in 1:N_sensores

            x1 = Sensores_cartesianos[k,1]
            y1 = Sensores_cartesianos[k,2]
            z1 = Sensores_cartesianos[k,3]

            Sensores[k,1] = (x1^2 + y1^2)^(1/2)
            Sensores[k,2] = atan(y1/x1)
            Sensores[k,3] = z1
        end

    elseif Coordenadas == "Cartesianas"
        Sensores =  Sensores_cartesianos
    end

    return Sensores
end

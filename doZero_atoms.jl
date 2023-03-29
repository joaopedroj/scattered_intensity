function getSphere(Radius, N, rₘ)
    r_new_atom, A = empty_arrays(N) 
    nValid = 0

    for iLoop = 1:10^8
        get_one_random_atom_inside_esfera!(Radius, r_new_atom) #aleatorização

        if is_valid_position(r_new_atom, A[1:nValid, : ], rₘ)  ## AQUI ##
            nValid += 1 
            A[nValid, : ] = r_new_atom
        end

        if nValid == N                                                                                                   # Aqui existe a contagem de atomos validos, dentro do disco
            break
        end
    end
    return A
end

function getCilindro(Radius, N, rₘ, h)
    r_new_atom, A = empty_arrays(N) 
    nValid = 0

    for iLoop = 1:10^8
        r_new_atom = ftn_AtomsOnCylinder(Radius, h, r_new_atom) #aleatorização

        if is_valid_position(r_new_atom, A[1:nValid, : ], rₘ)  ## AQUI ##
            nValid += 1 
            A[nValid, : ] = r_new_atom
        end

        if nValid == N                                                                                                   # Aqui existe a contagem de atomos validos, dentro do disco
            break
        end
    end
    return A
end

function empty_arrays(N::Integer) #gera as matrizes vazias

    r_new_atom = zeros(3)
    A = zeros(N,3)

    return r_new_atom, A
end

function get_one_random_atom_inside_esfera!(Radius::Number, r_new_atom::Array)
    kR = Radius
    U = rand()^(1 ./ 3.0)
    x = (2rand() - 1)
    y = (2rand() - 1)
    z = (2rand() - 1)
    mag = sqrt(x^2 + y^2 + z^2)

    r_new_atom[1] = kR*U*x/mag
    r_new_atom[2] = kR*U*y/mag
    r_new_atom[3] = kR*U*z/mag
    nothing
end

function get_one_random_atom_inside_cilindro!(Radius::Number, r_new_atom::Array)
    k=1
    λ = (2*π)/k
    h = 3*λ
    kR = Radius
    U = rand()^(1 ./ 3.0)
    x = (2rand() - 1)
    y = (2rand() - 1)
    z = (2rand() - 1)*h
    mag = sqrt(x^2 + y^2)

    r_new_atom[1] = kR*U*x/mag
    r_new_atom[2] = kR*U*y/mag
    r_new_atom[3] = z
    nothing
end

function ftn_AtomsOnCylinder(R::Number, h::Number, r_new_atom::Array)
    # Thanks to chatGPT for 'r = R*sqrt(rand())' command
    r, θ = R * sqrt(rand()), 2π * rand()

    x = r * cos(θ)
    y = r * sin(θ)
    z = -h * rand() + h / 2

    r_new_atom = [x, y, z]

    return r_new_atom 
end


function esferic_for_cart(esferic_coordinate::Array)

    azimuth = esferic_coordinate[1]
    alt = esferic_coordinate[2]
    radius = esferic_coordinate[3]


    x = radius * sin.(alt) * cos.(azimuth)
    y = radius * sin.(alt) * sin.(azimuth)
    z = radius * cos.(alt)

    position_cartesian = [x, y, z]
    # Esta função retorna as coordenadas dos pontos que antes eram polares para coordenadas cartesianas
    return position_cartesian
end



function is_valid_position(r_new_atom, A::Array, rₘ::Number) #Valida as posiçoes dos átomos
    return all(get_distance_C_to_D(A, r_new_atom) .≥ rₘ)
end

function get_distance_C_to_D(C::Array, D::Array) #Determina a distância entre os pontos anteriores e os recém criados
    n_rows = size(C, 1)                                    #Aqui temos a contagem dos N atomos novamente
    distanceCD = zeros(n_rows)                             #Define uma matriz que vai guardar a distancia entre os atomos
    @inbounds for i = 1:n_rows                             #Fazer o calculo em loop para cada atomo
        distanceCD[i] = Distances.evaluate(Euclidean(), C[i, :], D)
    end
    return distanceCD
end
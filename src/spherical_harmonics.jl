#
# This code was developed by Dr. Rodolfo Batista Negri at the National Institute for Space Research, Brazil.
#
# Description: 
# This function computes the potential and acceleration for a spacecraft's position using the 
# spherical harmonics method. The expansion is based on the coefficients of the 
# gravitational field of the central body (C and S coefficients), and the satellite's position 
# in space (r_vec). 
# Reference: as described in "Satellite Orbits" by Montenbruck (Page 66).
# Arguments:
# N      -> Order of the multipolar expansion (an integer).
# r_vec  -> A vector containing the Cartesian coordinates of the satellite (x, y, z).
# R0     -> The reference radius of the central body.
# C, S   -> Coefficients associated with the cosine and sine components of the expansion.
# mu     -> Gravitational parameter of the central body (GM, where G is the gravitational constant and M is the mass).

# Returns:
# F -> The total acceleration vector experienced by the satellite due to gravitational forces.
# U -> The gravitational potential at the satellite's location.



function spherical_harmonics(N,r_vec,R0,C,S,mu)
    # Pag. 66 - Montenbruck - Satellite orbits
    # N -> order of the expansion
    # θ -> latitude
    # λ -> longitude
    # R0-> reference radius
    # C and S -> coefficients associated with cosine and sine, respectively
    # mu-> gravitational parameter

    r = norm(r_vec)
    λ = atand( r_vec[2] , r_vec[1] )
    θ = asind( r_vec[3] / r )
    x = r_vec[1]
    y = r_vec[2]
    z = r_vec[3]
    # @show λ,θ

    u = sind(θ)
    P = zeros(N+2,N+2)
    V = zeros(N+2,N+2)
    V[1,1] = R0 / r
    W = zeros(N+2,N+2)
    P[1,1] = 1

    U = 0.0

    for mi = 1:N+2
        m = mi - 1
        for ni = maximum([2 mi]):N+2
            n = ni - 1
            if n == m # calculating sectorial terms
                P[ni,mi] = (2*m-1) * sqrt(1-u^2) * P[ni-1,mi-1]
            elseif n == m + 1
                P[ni,mi] = (2*m+1) * u * P[ni-1,mi]
            else
                P[ni,mi] = 1/(n-m) * ( (2*n-1)*u*P[ni-1,mi] - (n+m-1)*P[ni-2,mi] )
            end
            V[ni,mi] = (R0/r)^(n+1) * P[ni,mi] * cosd(m*λ)
            W[ni,mi] = (R0/r)^(n+1) * P[ni,mi] * sind(m*λ)
            if mi <= N + 1 && ni <= N + 1
                U = U + mu/R0 * ( C[ni,mi] * V[ni,mi] + S[ni,mi] * W[ni,mi] )
            end
            # @show n,m
            # @show P[ni,mi]
        end
    end
    U = U + mu / r
    # Nao excluir o que esta comentado abaixo (bom para checar e codigo ja pronto de recorrencia caso um dia precise)
    # P40 = 1/8*(35*u^4-30*u^2+3)
    # P41 = 5/2*(7*u^3-3*u)*sqrt(1-u^2)
    # P42 = 15/2*(7*u^2-1)*(1-u^2)
    # P43 = 105*u*(1-u^2)^(3/2)
    # P44 = 105*(1-u^2)^2
    # @show P40, P41, P42, P43, P44
    # @show W
    # Vtest = zeros(N+2,N+2)
    # Wtest = zeros(N+2,N+2)
    # Vtest[1,1] = R0/r
    # for mi = 1:N+2 # this calculates de V and W as recurrance
    #     m = mi - 1
    #     for ni = maximum([2 mi]):N+2
    #         n = ni - 1
    #         if n == m
    #             Vtest[ni,mi] = (2*m-1) * (x*R0/r^2*Vtest[ni-1,mi-1]-y*R0/r^2*Wtest[ni-1,mi-1])
    #             Wtest[ni,mi] = (2*m-1) * (x*R0/r^2*Wtest[ni-1,mi-1]+y*R0/r^2*Vtest[ni-1,mi-1])
    #         elseif n == m + 1
    #             Vtest[ni,mi] = (2*n-1)/(n-m) * z * R0/r^2 * V[ni-1,mi]
    #             Wtest[ni,mi] = (2*n-1)/(n-m) * z * R0/r^2 * W[ni-1,mi]
    #         else
    #             Vtest[ni,mi] = (2*n-1)/(n-m) * z * R0/r^2 * V[ni-1,mi] - (n+m-1)/(n-m) * R0^2/r^2 * V[ni-2,mi]
    #             Wtest[ni,mi] = (2*n-1)/(n-m) * z * R0/r^2 * W[ni-1,mi] - (n+m-1)/(n-m) * R0^2/r^2 * W[ni-2,mi]
    #         end
    #     end
    # end
    # @show Wtest
    # sleep(31231)


    xdd = zeros(N+1,N+1)
    ydd = zeros(N+1,N+1)
    zdd = zeros(N+1,N+1)
    for mi = 1:N+1
        m = mi - 1
        for ni = mi:N+1
            n = ni - 1
            if m == 0
                xdd[ni,mi] = mu / R0^2 * ( - C[ni,mi] * V[ni+1,mi+1] )
                ydd[ni,mi] = mu / R0^2 * ( - C[ni,mi] * W[ni+1,mi+1] )
            else
                xdd[ni,mi] = mu / R0^2 / 2 * ( ( - C[ni,mi] * V[ni+1,mi+1] - S[ni,mi] * W[ni+1,mi+1]  )
                + factorial(n-m+2)/factorial(n-m) * ( C[ni,mi] * V[ni+1,mi-1] + S[ni,mi] * W[ni+1,mi-1]  ) )
                ydd[ni,mi] = mu / R0^2 / 2 * ( ( - C[ni,mi] * W[ni+1,mi+1] + S[ni,mi] * V[ni+1,mi+1]  )
                + factorial(n-m+2)/factorial(n-m) * ( - C[ni,mi] * W[ni+1,mi-1] + S[ni,mi] * V[ni+1,mi-1]  ) )
            end
            zdd[ni,mi] = mu / R0^2 * ( (n-m+1) * ( - C[ni,mi] * V[ni+1,mi] - S[ni,mi] * W[ni+1,mi] ) )
        end
    end

    F = [sum(xdd); sum(ydd); sum(zdd)] # acceleration

    return F, U
end


module SmallBodyGravity

# Write your package code here.

function Multipolar_expansion(N,r_vec,R0,C,S,mu)
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


function Polyhedron_Model(p, r_vec)

    centroid_edges, centroid_facets, e_e, edges, facets, vertex, n_f, n_f_e, n_fp_e, r_e_1, r_e_2, r_f_1, r_f_2, r_f_3, G, sigma = p

    # Calculate normal versors for each facet used in the Edge summation
    n_f_Ee = n_f[Int.(edges[:,3]),:]
    n_fp_Ee = n_f[Int.(edges[:,4]),:]
    # Calculate normal versors for each facet used in the Facets summation
    n_f_Ff = n_f;
    # Distance from each edge's/facet's centroid to field point
    re = centroid_edges - repeat(r_vec,size(centroid_edges,1),1);
    rf = centroid_facets - repeat(r_vec,size(centroid_facets,1),1);
    # Distance from each edge's vertex to field point
    re1 = sqrt.(sum((r_e_1 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2));
    re2 = sqrt.(sum((r_e_2 - repeat(r_vec,size(centroid_edges,1),1)).^2,dims=2));
    # Distance from each facet's vertex to field point
    rf1 = r_f_1 - repeat(r_vec,size(centroid_facets,1),1);
    rf2 = r_f_2 - repeat(r_vec,size(centroid_facets,1),1);
    rf3 = r_f_3 - repeat(r_vec,size(centroid_facets,1),1);
    rf1_norm = sqrt.( sum( rf1.^2 , dims=2 ) );
    rf2_norm = sqrt.( sum( rf2.^2 , dims=2 ) );
    rf3_norm = sqrt.( sum( rf3.^2 , dims=2 ) );

    sum_e_U = 0.;
    sum_e_F = 0.;
    sum_f_U = 0.;
    sum_f_F = 0.;
    Lapl = 0.;
    for m = 1:size(centroid_edges,1) # Sum over edges

        Ee = n_f_Ee[m,:]*n_f_e[m,:]' + n_fp_Ee[m,:]*n_fp_e[m,:]'

        L_e = log.( ( re1[m] + re2[m] + e_e[m] ) / ( re1[m] + re2[m] - e_e[m] ));

        sum_e_U = sum_e_U + re[m,:]'*Ee*re[m,:]*L_e;

        sum_e_F = sum_e_F .+ Ee*re[m,:]*L_e;


        if m <= size(centroid_facets,1) # Sum over facets (inside edges' sum to save computational time)

            Ff = n_f_Ff[m,:]*n_f_Ff[m,:]'

            omega_f =  2 * atan( dot( rf1[m,:] , cross( rf2[m,:] , rf3[m,:] ) ) ,
            rf1_norm[m] * rf2_norm[m] * rf3_norm[m] + rf1_norm[m] * dot( rf2[m,:] , rf3[m,:] ) +
            rf2_norm[m] * dot( rf3[m,:] , rf1[m,:] ) + rf3_norm[m] * dot( rf1[m,:] , rf2[m,:] ));

            sum_f_U = sum_f_U + rf[m,:]'*Ff*rf[m,:]*omega_f;

            sum_f_F = sum_f_F .+ Ff*rf[m,:]*omega_f;

            Lapl = Lapl - G * sigma * omega_f; # Laplaciano

        end

    end
    U = G * sigma / 2 * ( sum_e_U - sum_f_U ); # Potential

    F = - G * sigma * ( sum_e_F - sum_f_F ); # Force

    return U, F, Lapl


end


end

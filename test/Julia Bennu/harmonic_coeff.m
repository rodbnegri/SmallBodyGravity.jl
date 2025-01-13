clc
clear
close all

% Obtem os coeficientes da expansão multipolar do potencial gravitacional,
% de ordem n, de um polihedro
% Método: Werner, RA - Spherical harmonic coefficients for the potential of a constant-density polyhedron

%%% TODAS AS VARIAVEIS COM "_var" está adicionado +1 nelas, pois não é
%%% possível criar uma matriz que acesse elemento 0

%%% Se o C00 não der exatamente 1, quer dizer que a massa não está
%%% totalmente coerente com o volume e densidade

%%% ATUALIZAR PARA FICAR MAIS RAPIDO USANDO >>>MULTINOMIAL THEOREM<<<

addpath('C:\Users\rodol\Documents\Programas MATLAB')

load('facets.dat')
load('vertex.dat')

vertex = vertex * 1e3;
radius_each_vertex = sqrt(sum(vertex.^2,2));
a = max( radius_each_vertex );
syms X Y Z
M = 6.68700000e+15  ;
sigma = 2.68176926e+03  ;
n_max = 5;
disp('pode rodar? Já atualizou a massa e densidade do asteroide? E a qtd de harmonicos?')
pause

% calcular volume e obter densidade daqui (massa dada)
% colocar no centroide e alinha com eixos principais

c_symb = sym('C',[n_max+1,n_max+1]);
s_symb = sym('S',[n_max+1,n_max+1]);
c_symb(1,1) = 1;
s_symb(1,1) = 0;
C = zeros(n_max+1,n_max+1);
S = zeros(n_max+1,n_max+1);
C_sum_simplices = C;
S_sum_simplices = S;
C(1,1) = 1;
S(1,1) = 0;
for simplice = 1:size(facets,1)
    
    vertex_1 = vertex( facets( simplice, 1 ), : ) ;
    vertex_2 = vertex( facets( simplice, 2 ), : ) ;
    vertex_3 = vertex( facets( simplice, 3 ), : ) ;

    x1=vertex_1(1);
    x2=vertex_2(1);
    x3=vertex_3(1);
    y1=vertex_1(2);
    y2=vertex_2(2);
    y3=vertex_3(2);
    z1=vertex_1(3);
    z2=vertex_2(3);
    z3=vertex_3(3);
    
    xp = x1*X+x2*Y+x3*Z;
    yp = y1*X+y2*Y+y3*Z;
    zp = z1*X+z2*Y+z3*Z;
    
    rp = sqrt( xp^2 + yp^2 + zp^2 );
    
    for n_var = 1:n_max+1
        for m_var = 1:n_var
            
            
            n = n_var - 1;
            m = m_var - 1;
            
            % Eqs. 7-9 in Werner
            if n >= m
                % sectorial
                if n == m
                    if n == 1
                        c_symb(2,2) = 1 / sqrt(3) * xp / a; 
                        s_symb(2,2) = 1 / sqrt(3) * yp / a; 
                    elseif n > 1
                        c_symb(n_var,m_var) = (2*n-1) / sqrt(2*n*(2*n+1)) * ( xp/a * c_symb(n_var-1,n_var-1) - yp/a * s_symb(n_var-1,n_var-1) );
                        s_symb(n_var,m_var) = (2*n-1) / sqrt(2*n*(2*n+1)) * ( yp/a * c_symb(n_var-1,n_var-1) + xp/a * s_symb(n_var-1,n_var-1) );
                    end
                end
                
                % subdiagonal
                if m == n-1
                    c_symb(n_var,m_var) = (2*n-1) / sqrt(2*n+1) * zp / a * c_symb(n_var-1,m_var);
                    s_symb(n_var,m_var) = (2*n-1) / sqrt(2*n+1) * zp / a * s_symb(n_var-1,m_var);
                end
                
                % vertical
                if n ~= m && m ~= n-1
                    c_symb(n_var,m_var) = (2*n-1) * sqrt( (2*n-1) / ((2*n+1)*(n+m)*(n-m)) ) * zp/a * c_symb(n_var-1,m_var) - ...
                        sqrt( ((2*n-3)*(n+m-1)*(n-m-1)) / ((2*n+1)*(n+m)*(n-m)) ) * (rp/a)^2 * c_symb(n_var-2,m_var) ;
                    s_symb(n_var,m_var) = (2*n-1) / sqrt(2*n+1) * zp / a * s_symb(n_var-1,m_var);
                end
                
                coef_factorial = zeros(n_max+1,n_max+1,n_max+1);
                alpha = zeros(n_max+1,n_max+1,n_max+1);
                beta = zeros(n_max+1,n_max+1,n_max+1);
                for i_var = 1:n_var % finding the coeffs of the trinomials
                    i = i_var - 1;
                    for j_var = 1:n_var-i
                        j = j_var - 1;
                        k = n - i - j;
                        k_var = k + 1;
                        
                        A = c_symb(n_var,m_var);
                        B = s_symb(n_var,m_var);
                        
                        [CoefsA, TrinomialA] = coeffs(A, [X Y Z]);
                        [CoefsB, TrinomialB] = coeffs(B, [X Y Z]);
                        
                        alpha(i_var,j_var,k_var) = CoefsA( TrinomialA==X^i*Y^j*Z^k  );
                        if isempty(TrinomialB)
                            beta(i_var,j_var,k_var) = 0;
                        else
                            beta(i_var,j_var,k_var) = CoefsB( TrinomialB==X^i*Y^j*Z^k  );
                        end
                        
                        coef_factorial(i_var,j_var,k_var) = factorial(i) * factorial(j) * factorial(k) ;
                    end
                end
                
                J = [x1 x2 x3; y1 y2 y3; z1 z2 z3];
                C(n_var,m_var) = det(J)/factorial(n+3) * sum(sum(sum( coef_factorial .* alpha )));
                S(n_var,m_var) = det(J)/factorial(n+3) * sum(sum(sum( coef_factorial .* beta )));
            end
            
        end
    end
    
    clc
    disp(simplice/size(facets,1)*100)
    C_sum_simplices = C_sum_simplices + sigma / M * C % these are the bar C and S (normalized C and S) 
    S_sum_simplices = S_sum_simplices + sigma / M * S
    
end

n_matrix = repmat([0:n]',1,n+1);
m_matrix = n_matrix';
n_minus_m = tril(n_matrix-m_matrix);
n_plus_m = tril(n_matrix+m_matrix);
N_nm = sqrt( 2 * (2*n_matrix+1) .* factorial(n_minus_m) ./ factorial(n_plus_m) );
N_nm(1,1) = sqrt( 1 * (2*0+1) * factorial(0-0) / factorial(0+0) );


C_non_norml = C_sum_simplices.*N_nm;
S_non_norml = S_sum_simplices.*N_nm;

fID = fopen('spherical_harmonics.txt','w');
fprintf(fID,'Mass [kg]: %12.8e \n\n',M);
fprintf(fID,'Density [kg/m3]: %12.8e \n\n',sigma);
fprintf(fID,'Reference distance [m]: %12.8e \n\n',a);
fprintf(fID,'Non-normalized coefficients:\n\n');
fprintf(fID,'%6s %6s \t %12s \t %12s\n','n','m','Cnm','Snm');
for n_var = 1:n_max+1
    for m_var = 1:n_var
        fprintf(fID,'%6d %6d \t %12.8e \t %12.8e\n',n_var-1,m_var-1,C_non_norml(n_var,m_var),S_non_norml(n_var,m_var));
    end
end
fprintf(fID,'\n\n\n\nNormalized coefficients:\n\n');
fprintf(fID,'%6s %6s \t %12s \t %12s\n','n','m','Cnm','Snm');
for n_var = 1:n_max+1
    for m_var = 1:n_var
        fprintf(fID,'%6d %6d \t %12.8e \t %12.8e\n',n_var-1,m_var-1,C_sum_simplices(n_var,m_var),S_sum_simplices(n_var,m_var));
    end
end
fprintf(fID,'\n\n\n\n-------------- Normalized coefficients: --------------\n\n');
for n_var = 1:n_max+1
    for m_var = 1:n_var
        fprintf(fID,'C[%2d,%2d] = %12.8e \n',n_var,m_var,C_sum_simplices(n_var,m_var));
        fprintf(fID,'S[%2d,%2d] = %12.8e \n',n_var,m_var,S_sum_simplices(n_var,m_var));
    end
end
fprintf(fID,'\n\n\n\n-------------- Non-normalized coefficients: --------------\n\n');
for n_var = 1:n_max+1
    for m_var = 1:n_var
        fprintf(fID,'C[%2d,%2d] = %12.8e \n',n_var,m_var,C_non_norml(n_var,m_var));
        fprintf(fID,'S[%2d,%2d] = %12.8e \n',n_var,m_var,S_non_norml(n_var,m_var));
    end
end
fclose(fID);
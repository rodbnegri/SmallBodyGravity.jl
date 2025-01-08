function vertex = polyhedron_correction(vertex_old, faces, M, filename_IN)
% Models from: Dobrovolskis - Inertia of Any Polyhedron


% Finding volume and corresponding density
Volume = 0;
for simplice = 1:size(faces,1)
    vertex_1 = vertex_old( faces( simplice, 1 ), : ) ;
    vertex_2 = vertex_old( faces( simplice, 2 ), : ) ;
    vertex_3 = vertex_old( faces( simplice, 3 ), : ) ;
    Volume = dot(cross(vertex_2,vertex_3),vertex_1) / 6 + Volume ;
end
sigma = M / Volume;

% Fiding the center of mass
R = zeros(1,3);
for simplice = 1:size(faces,1)
    vertex_1 = vertex_old( faces( simplice, 1 ), : ) ;
    vertex_2 = vertex_old( faces( simplice, 2 ), : ) ;
    vertex_3 = vertex_old( faces( simplice, 3 ), : ) ;
    dV = dot(cross(vertex_2,vertex_3),vertex_1) / 6; % volume of the current simplex
    dR = (vertex_2+vertex_3+vertex_1) / 4; % center of mass of the current simplex
    R = dV * dR / Volume + R ;
end

% Finding the products of inertia and the inertia matrix
P = zeros(3,3);
for simplice = 1:size(faces,1)
    for j = 1:3
        for k = 1:3
            D = vertex_old( faces( simplice, 1 ), : ) ;
            E = vertex_old( faces( simplice, 2 ), : ) ;
            F = vertex_old( faces( simplice, 3 ), : ) ;
            dV = dot(cross(E,F),D) / 6; % volume of the current simplex
            P(j,k) = sigma * dV / 20 * ( 2*D(j)*D(k) + 2*E(j)*E(k) + 2*F(j)*F(k) + ...
                D(j)*E(k) + D(k)*E(j) + D(j)*F(k) + D(k)*F(j) + E(j)*F(k) + E(k)*F(j) ) + P(j,k);
        end
    end
end
I = [P(2,2)+P(3,3) -P(1,2) -P(1,3)
    -P(1,2) P(1,1)+P(3,3) -P(2,3)
    -P(1,3) -P(2,3) P(1,1)+P(2,2)];

% Teorema do eixo paralelo, finding the inertia tensor relative to the
% center of mass:
ICM = I - M * [R(2)^2+R(3)^2 -R(1)*R(2) -R(1)*R(3);
                -R(1)*R(2) R(1)^2+R(3)^2 -R(2)*R(3)
                -R(1)*R(3) -R(2)*R(3) R(1)^2+R(2)^2];
            
% Eigenvalues and eigenvectors
[eig_vec, eig_val] = eig(ICM); % each column is an eigen vector, s.t.: ICM*eig_vec = eig_vec*eig_val
eig_vec(:,diag(eig_vec)<0) = -eig_vec(:,diag(eig_vec)<0); 

% Aligning with the principal axes
vertex = zeros(size(vertex_old));
for simplice = 1:size(vertex_old,1)
    vertex(simplice,:) = ( vertex_old(simplice,:) - R ) * eig_vec ;
end

% Calculating equivalent ellipsoid
A = eig_val(1,1);
B = eig_val(2,2);
C = eig_val(3,3);
a = sqrt( 5 * (B+C-A) / 2 / M );
b = sqrt( 5 * (A+C-B) / 2 / M );
c = sqrt( 5 * (B+A-C) / 2 / M );
%% Finding new properties

% Fiding the center of mass
Rnew = zeros(1,3);
for simplice = 1:size(faces,1)
    vertex_1 = vertex( faces( simplice, 1 ), : ) ;
    vertex_2 = vertex( faces( simplice, 2 ), : ) ;
    vertex_3 = vertex( faces( simplice, 3 ), : ) ;
    dV = dot(cross(vertex_2,vertex_3),vertex_1) / 6; % volume of the current simplex
    dR = (vertex_2+vertex_3+vertex_1) / 4; % center of mass of the current simplex
    Rnew = dV * dR / Volume + Rnew ;
end

% Finding the products of inertia and the inertia matrix
P = zeros(3,3);
for simplice = 1:size(faces,1)
    for j = 1:3
        for k = 1:3
            D = vertex( faces( simplice, 1 ), : ) ;
            E = vertex( faces( simplice, 2 ), : ) ;
            F = vertex( faces( simplice, 3 ), : ) ;
            dV = dot(cross(E,F),D) / 6; % volume of the current simplex
            P(j,k) = sigma * dV / 20 * ( 2*D(j)*D(k) + 2*E(j)*E(k) + 2*F(j)*F(k) + ...
                D(j)*E(k) + D(k)*E(j) + D(j)*F(k) + D(k)*F(j) + E(j)*F(k) + E(k)*F(j) ) + P(j,k);
        end
    end
end
Inew = [P(2,2)+P(3,3) -P(1,2) -P(1,3)
    -P(1,2) P(1,1)+P(3,3) -P(2,3)
    -P(1,3) -P(2,3) P(1,1)+P(2,2)];

%% Writing file
fID = fopen('polyhedron_properties.txt','w');
fprintf(fID,'Number of faces: %d \n\n',size(faces,1));
fprintf(fID,'Mass [kg]: %12.8e \n\n',M);
fprintf(fID,'Volume [km3]: %12.8e \n\n',Volume);
fprintf(fID,'Density [kg/m3]: %12.8e \n\n',sigma*1e-9);
fprintf(fID,'Inertia matrix [kg.km2]:\n');
fprintf(fID,'%12.8e %12.8e %12.8e\n',Inew');
fprintf(fID,'\n Center of mass [km]:\n');
fprintf(fID,'%12.8e %12.8e %12.8e\n',Rnew');
fprintf(fID,'\n Equivalent ellipsoid:\n');
fprintf(fID,'a x b x c [km]: %6.4f x %6.4f x %6.4f\n',a,b,c);
fprintf(fID,'Equiv. volume [km3]: %6.4f\n',4/3*pi*a*b*c);
fprintf(fID,'Equiv. density [kg/m3]: %6.4f\n',3*M/4/pi/a/b/c*1e-9);


fprintf(fID,'\n\n\n---------- OLD PROPERTIES AS FOUND IN "%s" ----------\n', filename_IN);
fprintf(fID,'Inertia matrix [kg.km^2]:\n');
fprintf(fID,'%12.8e %12.8e %12.8e\n',I');
fprintf(fID,'\n Center of mass [km]:\n');
fprintf(fID,'%12.8e %12.8e %12.8e\n',R');

% fprintf(fID,'Non-normalized coefficients:\n\n');
% fprintf(fID,'%6s %6s \t %12s \t %12s\n','n','m','Cnm','Snm');
% for n_var = 1:n_max+1
%     for m_var = 1:n_var
%         fprintf(fID,'%6d %6d \t %12.8e \t %12.8e\n',n_var-1,m_var-1,C_non_norml(n_var,m_var),S_non_norml(n_var,m_var));
%     end
% end
fclose(fID);

end
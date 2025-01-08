
% ------------------------------------------------------------------------
% Code Developed by Dr. Rodolfo Batista Negri
%
% Based on the work of  Werner and Scheeres [1]. This code should not be shared
% or distributed without prior consent from Dr. Negri (rodolfobnegri@yahoo.com.br).
%
% If you use this code in any manner, whether integrally, partially, or as
% inspiration, please cite references [2] and [3]. For instance:
% "... used the code developed by Negri [2,3], which applies the asteroid model 
% of Werner and Scheeres [1]."
%
% Please note that Dr. Negri does not assume any responsibility for the misuse 
% of this code. He is unable to provide assistance or support for its implementation.
%
% [1] Werner, R.A. and Scheeres, D.J., 1996. Exterior gravitation of a polyhedron 
% derived and compared with harmonic and mascon gravitation representations of 
% asteroid 4769 Castalia. Celestial Mechanics and Dynamical Astronomy, 65, pp.313-344.
%
% [2] Negri, R. B., “A Study of Dynamics, Guidance, Navigation, and Control 
% Applied to Asteroid Deflection," Ph.D. thesis, National Institute for Space Research,
% São José dos Campos, 2022.
%
% [3] Negri, R. B. and Prado, A.F. B. A., 2022. Autonomous and Robust Orbit-Keeping 
% for Small-Body Missions. Journal of Guidance, Control, and Dynamics, 45(3), pp.587-598.
% ------------------------------------------------------------------------
function [faces, vertex, edges, r_f_1, r_f_2, r_f_3, r_e_1, r_e_2, e_e,...
    n_f, centroid_faces, centroid_edges, n_f_e, n_fp_e] ...
    = Polyhedron_data(filename_IN, folder_path, Edge_Coord_file,M)

file_edge = [folder_path '\' Edge_Coord_file];
file_in = [folder_path '\' filename_IN];

% coords da edge, duas primeiras colunas dão os vertex de cada edge,
% duas colunas restantes dão as faces que formam a edge
edges = load(file_edge); 

faces = load('faces.dat');
vertex = load('vertex.dat');

% Coordinates of each vertex of each face
r_f_1 =  vertex(faces(:,1),:);
r_f_2 =  vertex(faces(:,2),:);
r_f_3 =  vertex(faces(:,3),:);

% Coordinates of each vertex of each edge
r_e_1 =  vertex(edges(:,1),:);
r_e_2 =  vertex(edges(:,2),:);

% Size of each edge
e_e_vec = r_e_2 - r_e_1;
e_e = sqrt( sum( (r_e_2 - r_e_1).^2, 2) );

% Corresponding faces for each edge
f_e = faces(edges(:,3),:);
fp_e = faces(edges(:,4),:);

% Normal vector to each face
norma = sqrt(sum(cross(r_f_2 - r_f_1,r_f_3-r_f_2).^2,2));
n_f = cross(r_f_2 - r_f_1,r_f_3-r_f_2) ./ [norma norma norma];

% Find centroid of each face and edge (it could be any point on both)
centroid_faces = ( r_f_1 + r_f_2 + r_f_3 ) / 3;
centroid_edges = ( r_e_1 + r_e_2 ) / 2;

% Find normals to each edge
norma = sqrt(sum( (centroid_edges - centroid_faces(edges(:,3),:)).^2,2));
nn_paralel = (centroid_edges - centroid_faces(edges(:,3),:)) ./ [norma norma norma]; % versor from centroid of the face to the edge's centroid
nn_perpendicular = cross(nn_paralel,e_e_vec); % vector perpendicular to edge and face
norma = sqrt( sum(cross(nn_perpendicular, e_e_vec).^2,2) );
n_f_e = cross(nn_perpendicular, e_e_vec) ./ [norma norma norma]; % versor parallel to face and perpendicular to edge
theta = acosd( sum(nn_paralel.*n_f_e,2) );
n_f_e(  theta > 90 , : ) = - n_f_e(  theta > 90 , : ); % making sure n_f_e is pointing outside its face
norma = sqrt(sum( (centroid_edges - centroid_faces(edges(:,4),:)).^2,2));
nn_paralel = (centroid_edges - centroid_faces(edges(:,4),:)) ./ [norma norma norma]; % versor from centroid of the face to the edge's centroid
nn_perpendicular = cross(nn_paralel,e_e_vec); % vector perpendicular to edge and face
norma = sqrt( sum(cross(nn_perpendicular, e_e_vec).^2,2) );
n_fp_e = cross(nn_perpendicular, e_e_vec) ./ [norma norma norma]; % versor no parallel to face and perpendicular to edge
theta = acosd( sum(nn_paralel.*n_fp_e,2) );
n_fp_e(  theta > 90 , : ) = - n_fp_e(  theta > 90 , : ); % making sure n_f_e is pointing outside its face

% Create data for other programs
f_edges = fopen('edges.dat','w');
f_r_f_1 = fopen('r_f_1.dat','w');
f_r_f_2 = fopen('r_f_2.dat','w');
f_r_f_3 = fopen('r_f_3.dat','w');
f_r_e_1 = fopen('r_e_1.dat','w');
f_r_e_2 = fopen('r_e_2.dat','w');
f_e_e = fopen('e_e.dat','w');
f_n_f = fopen('n_f.dat','w');
f_centroid_faces = fopen('centroid_faces.dat','w');
f_centroid_edges = fopen('centroid_edges.dat','w');
f_n_f_e = fopen('n_f_e.dat','w');
f_n_fp_e = fopen('n_fp_e.dat','w');

fprintf(f_edges,'%.20e %.20e %.20e %.20e\n',edges');
fprintf(f_r_f_1,'%.20e %.20e %.20e\n',r_f_1');
fprintf(f_r_f_2,'%.20e %.20e %.20e\n',r_f_2');
fprintf(f_r_f_3,'%.20e %.20e %.20e\n',r_f_3');
fprintf(f_r_e_1,'%.20e %.20e %.20e\n',r_e_1');
fprintf(f_r_e_2,'%.20e %.20e %.20e\n',r_e_2');
fprintf(f_e_e,'%.20e\n',e_e);
fprintf(f_n_f,'%.20e %.20e %.20e\n',n_f');
fprintf(f_centroid_faces,'%.20e %.20e %.20e\n',centroid_faces');
fprintf(f_centroid_edges,'%.20e %.20e %.20e\n',centroid_edges');
fprintf(f_n_f_e,'%.20e %.20e %.20e\n',n_f_e');
fprintf(f_n_fp_e,'%.20e %.20e %.20e\n',n_fp_e');
fclose all;

% % check if n_f is correct
% quiver3(centroid_faces(:,1),centroid_faces(:,2),centroid_faces(:,3),...
%     n_f(:,1),n_f(:,2),n_f(:,3))
% axis equal
% grid on
% return

% % check if n_f_e and n_fp_e are correct
% index = round( rand * size(edges,1) ); 
% POS_TRI = [r_f_1(edges(index,3),:);r_f_2(edges(index,3),:);r_f_3(edges(index,3),:)];
% patch(POS_TRI(:,1),POS_TRI(:,2),POS_TRI(:,3),'b')
% hold on
% quiver3(centroid_edges(index,1),centroid_edges(index,2),centroid_edges(index,3),...
%     n_f_e(index,1),n_f_e(index,2),n_f_e(index,3),'b')
% axis equal
% POS_TRI = [r_f_1(edges(index,4),:);r_f_2(edges(index,4),:);r_f_3(edges(index,4),:)];
% patch(POS_TRI(:,1),POS_TRI(:,2),POS_TRI(:,3),'r')
% quiver3(centroid_edges(index,1),centroid_edges(index,2),centroid_edges(index,3),...
%     n_fp_e(index,1),n_fp_e(index,2),n_fp_e(index,3),'r') 
% grid on
% return
end
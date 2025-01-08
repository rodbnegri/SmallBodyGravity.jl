clc
clear
close all


M = 7.329e10;

file_obj = 'bennu.obj';
pasta = pwd;

Create_EdgeCoordinates(file_obj, 'Edge_Coords.dat', pasta, M);

[faces, vertex, edges, r_f_1, r_f_2, r_f_3, r_e_1, r_e_2, e_e,...
    n_f, centroid_faces, centroid_edges, n_f_e, n_fp_e] ...
    = Polyhedron_data(file_obj, pasta, 'Edge_Coords.dat');
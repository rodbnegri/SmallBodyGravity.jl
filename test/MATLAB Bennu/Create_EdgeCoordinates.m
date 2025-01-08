function [] = Create_EdgeCoordinates(filename_IN, filename_OUT, folder_path,M)

% Esse programa gera arquivos para representar os edges do polihedro.obj, as
% duas primeiras componentes mostram os vertex que formam o edge e as duas
% últimas mostram os faces

file_in = [folder_path '\' filename_IN];
file_out =  [folder_path '\' filename_OUT];

fid  = fopen(file_in, 'r') ;
data = cell(1e6, 4) ;                    % Prealloc.
rCnt = 0 ;                               % Row counter.
while ~feof(fid)
    rCnt = rCnt + 1 ;
    data{rCnt,1} = fscanf(fid, '%s', 1) ;
    data{rCnt,2} = fscanf(fid, '%f', 1) ;
    data{rCnt,3} = fscanf(fid, '%f', 1) ;
    data{rCnt,4} = fscanf(fid, '%f', 1) ;
end
fclose(fid) ;
data = data(1:rCnt,:) ;                  % Truncate.

m = size(data,1);
n = size(data,2);
vert_or_fac = cell2mat(data(1:m,1));
v = cell2mat( data( vert_or_fac == 'v' , 2:n ) );
f = cell2mat( data( vert_or_fac == 'f' , 2:n ) );

% Function to position the reference system in its center of mass and, align
% with the principal axis of inertia and find the density
v = polyhedron_correction(v, f, M, filename_IN); % new vertex

qtd_edges = (size(v,1) + size(f,1))*2 ;% size(v,1) + size(f,1) - 2 ; Euler's formula only works for convex polyhedron
e_alfa = zeros( qtd_edges , 4 );
h = waitbar(0,'Espera aí, tio! Tô achando as coordenadas das Edges do poliedro!');
s_antes = 0;
for i = 1 : size(f,1) % running in each face
    
    waitbar(i/size(f,1))
    
    AA = [];
    
    nff = size( find( e_alfa( : , 4 ) == i) , 1 ); % number of contact faces already found for the current face
    
    for j = i : size(f,1) % from the current face to the end
        
        int = intersect( f(i,:) , f(j,:) ); % find contact faces
        
        if size(int,2) == 2
            
            AA = [AA; int i j]; % armazenar as contact faces
            
        end
        
        if size( AA , 1 ) + nff == 3
            break % terminar loop se já achou as 3 faces que estao em contato
        end
        
    end
    
    s_depois = size(AA,1) + s_antes; % isso é só para ir atualizando as rows em e_alfa
    if size(AA,1) ~= 0
        e_alfa( s_antes+1:s_depois , : ) = AA;
    end
    s_antes = s_depois;
    
end
close(h)
e_alfa( ~any(e_alfa,2), : ) =  []; % deletando as linhas com zeros (como a formula 
% de euler nao funciona pra polihedros nao convexos, tem que fazer um 
% e_alfa grande e deletar o que sobrou de 0)

fileID = fopen(file_out,'w');
fprintf(fileID,'%d %d %d %d\n',e_alfa');

f_faces = fopen('faces.dat','w');
f_vertex = fopen('vertex.dat','w');
fprintf(f_faces,'%.20e %.20e %.20e\n',f');
fprintf(f_vertex,'%.20e %.20e %.20e\n',v');
fclose all;

end
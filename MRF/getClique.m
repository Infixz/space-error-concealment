%Function that gets the clique and the complementary clique of a pixel.

function cliques = getClique(bloque, r, c)

% 3x5 masks for clique and complementary clique as illustrated in Fig.1 %%%
maskC = [0 5 0 3 0;
         7 6 4 2 1;
         0 8 0 0 0];
     
maskCP = [0 0 0 8 0;
          1 2 4 6 7;
          0 3 0 5 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cliques = [c; c']
cliques = -ones(2,8);

%Clique
for i = -2:0
    for j = -2:2
        if bloque(r+i,c+j) >= 0 && ~(i == 0 && j == 0) && maskC(i+3,j+3) > 0
            cliques(1,maskC(i+3,j+3)) = bloque(r+i,c+j);
        end
    end
end
    
%Clique'
for i = 0:2
    for j = -2:2
        if bloque(r+i,c+j) >= 0 && ~(i == 0 && j == 0) && maskCP(i+1,j+3) > 0
            cliques(2,maskCP(i+1,j+3)) = bloque(r+i,c+j);
        end
    end
end



end

function [block yAux] = biLinInterpol_v5(mask, center, angle, mb_size)

r = center(1);
c = center(2);
y = mask;

[d1 d2] = size(mask);
yAux = zeros(d1,d2,3);
yAux(:,:,1) = mask;
yAux(:,:,2) = mask;
yAux(:,:,3) = mask;
block = -ones(mb_size);

if angle < 90
    
    for i = 1:mb_size
        for j = 1:mb_size
            slope = tan(2*pi/360*(90 - angle));
            above = i;
            below = mb_size - i + 1;
            
            left = round(above*slope);
            right = round(below*slope);
            
            %FROM TOP
            if left <= j
                %TO BOTTOM
                if right <= mb_size - j + 1
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);                   
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r-1,c+j-1-left) + d1/(d1 + d2)*y(r+mb_size, c+j-1+right);
                    block(i,j) = d2/(d1 + d2)*y(r-1,c+j-1-left) + d1/(d1 + d2)*y(r+mb_size, c+j-1+right);
%                     yAux(r-1,c+j-1-left,1) = 255;
%                     yAux(r-1,c+j-1-left,2) = 0;
%                     yAux(r-1,c+j-1-left,3) = 0;
%                     yAux(r+mb_size, c+j-1+right,1) = 0;
%                     yAux(r+mb_size, c+j-1+right,2) = 255;
%                     yAux(r+mb_size, c+j-1+right,3) = 0;
                    
                
                %TO RIGHT    
                else
                    slope = tan(2*pi/360*angle);
                    aux = (mb_size - j + 1)*slope;
                    below = round(aux);
                    right = mb_size - j + 1;
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r-1,c+j-1-left) + d1/(d1 + d2)*y(r+i-1+below, c+mb_size);
                    block(i,j) = d2/(d1 + d2)*y(r-1,c+j-1-left) + d1/(d1 + d2)*y(r+i-1+below, c+mb_size);
%                     yAux(r-1,c+j-1-left,1) = 255;
%                     yAux(r-1,c+j-1-left,2) = 0;
%                     yAux(r-1,c+j-1-left,3) = 0;
%                     yAux(r+i-1+below, c+mb_size,1) = 0;
%                     yAux(r+i-1+below, c+mb_size,2) = 255;
%                     yAux(r+i-1+below, c+mb_size,3) = 0;
                end
                
            %FROM LEFT    
            else
                slope = tan(2*pi/360*angle);
                aux = (slope*j);
                above = round(aux);
                left = j;
                %TO BOTTOM
                if right <= mb_size - j + 1
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r+i-1-above,c-1) + d1/(d1 + d2)*y(r+mb_size, c+j-1+right);
                    block(i,j) = d2/(d1 + d2)*y(r+i-1-above,c-1) + d1/(d1 + d2)*y(r+mb_size, c+j-1+right);
                    yAux(r+i-1-above,c-1,1) = 255;
                    yAux(r+i-1-above,c-1,2) = 0;
                    yAux(r+i-1-above,c-1,3) = 0;
                    yAux(r+mb_size, c+j-1+right,1) = 0;
                    yAux(r+mb_size, c+j-1+right,2) = 255;
                    yAux(r+mb_size, c+j-1+right,3) = 0;

                %TO RIGHT
                else
                    slope = tan(2*pi/360*angle);
                    aux = (mb_size - j + 1)*slope;
                    below = round(aux);
                    right = mb_size - j + 1;
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r+i-1-above,c-1) + d1/(d1 + d2)*y(r+i-1+below, c+mb_size);
                    block(i,j) = d2/(d1 + d2)*y(r+i-1-above,c-1) + d1/(d1 + d2)*y(r+i-1+below, c+mb_size);
%                     yAux(r+i-1-above,c-1,1) = 255;
%                     yAux(r+i-1-above,c-1,2) = 0;
%                     yAux(r+i-1-above,c-1,3) = 0;
%                     yAux(r+i-1+below, c+mb_size,1) = 0;
%                     yAux(r+i-1+below, c+mb_size,2) = 255;
%                     yAux(r+i-1+below, c+mb_size,3) = 0;
                end                  
            end
        end
    end
    
else
    for i = 1:mb_size
        for j = 1:mb_size
            slope = tan(2*pi/360*(angle - 90));
            above = i;
            below = mb_size - i + 1;
            
            right = round(above*slope);
            left = round(below*slope);
            
            %FROM TOP
            if right <= mb_size - j + 1
                %TO BOTTOM
                if left <= j
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r-1,c+j-1+right) + d1/(d1 + d2)*y(r+mb_size, c+j-1-left);
                    block(i,j) = d2/(d1 + d2)*y(r-1,c+j-1+right) + d1/(d1 + d2)*y(r+mb_size, c+j-1-left);
%                     yAux(r-1,c+j-1+right,1) = 255;
%                     yAux(r-1,c+j-1+right,2) = 0;
%                     yAux(r-1,c+j-1+right,3) = 0;
%                     yAux(r+mb_size, c+j-1-left,1) = 0;
%                     yAux(r+mb_size, c+j-1-left,2) = 255;
%                     yAux(r+mb_size, c+j-1-left,3) = 0;
                
                %TO LEFT
                else
                    slope = tan(2*pi/360*(90-(angle-90)));
                    aux = j*slope;
                    below = round(aux);
                    left = j;
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r-1,c+j-1+right) + d1/(d1 + d2)*y(r+i-1+below, c-1);
                    block(i,j) = d2/(d1 + d2)*y(r-1,c+j-1+right) + d1/(d1 + d2)*y(r+i-1+below, c-1);
%                     yAux(r-1,c+j-1+right,1) = 255;
%                     yAux(r-1,c+j-1+right,2) = 0;
%                     yAux(r-1,c+j-1+right,3) = 0;
%                     yAux(r+i-1+below, c-1,1) = 0;
%                     yAux(r+i-1+below, c-1,2) = 255;
%                     yAux(r+i-1+below, c-1,3) = 0;
                end
                
            %FROM RIGHT    
            else
                slope = tan(2*pi/360*(90-(angle-90)));
                aux = (mb_size - j + 1)*slope;
                above = round(aux);
                right = mb_size - j + 1;
                %TO BOTTOM
                if left <= j
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r+i-1-above,c+mb_size) + d1/(d1 + d2)*y(r+mb_size, c+j-1-left);
                    block(i,j) = d2/(d1 + d2)*y(r+i-1-above,c+mb_size) + d1/(d1 + d2)*y(r+mb_size, c+j-1-left);
%                     yAux(r+i-1-above,c+mb_size,1) = 255;
%                     yAux(r+i-1-above,c+mb_size,2) = 0;
%                     yAux(r+i-1-above,c+mb_size,3) = 0;
%                     yAux(r+mb_size, c+j-1-left,1) = 0;
%                     yAux(r+mb_size, c+j-1-left,2) = 255;
%                     yAux(r+mb_size, c+j-1-left,3) = 0;
                
                %TO LEFT
                else
                    slope = tan(2*pi/360*(90-(angle-90)));
                    aux = j*slope;
                    below = round(aux);
                    left = j;
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    y(r+i-1, c+j-1) = d2/(d1 + d2)*y(r+i-1-above,c+mb_size) + d1/(d1 + d2)*y(r+i-1+below, c-1);
                    block(i,j) = d2/(d1 + d2)*y(r+i-1-above,c+mb_size) + d1/(d1 + d2)*y(r+i-1+below, c-1);
%                     yAux(r+i-1-above,c+mb_size,1) = 255;
%                     yAux(r+i-1-above,c+mb_size,2) = 0;
%                     yAux(r+i-1-above,c+mb_size,3) = 0;
%                     yAux(r+i-1+below, c-1,1) = 0;
%                     yAux(r+i-1+below, c-1,2) = 255;
%                     yAux(r+i-1+below, c-1,3) = 0;
                end
                    
            end
            
        end
    end
end



end

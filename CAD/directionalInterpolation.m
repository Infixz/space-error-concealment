%This funtion conceals the missing block by applying directional
%interpolation.

function block_conc = directionalInterpolation(support_area, angle, mb_size)

% INPUT:
%   support_area - area around the missing macroblock (included), that also
%                contains correctly received and decoded pixels
%   angle        - angle (in degrees) that indicates the direction of
%                interpolation, NOTE: angles are counted clockwise (i.e.
%                45deg is located in the second quadrant and 135deg in the
%                first quadrant)
%   mb_size      - macroblock dimensions (in pixels)
%
% OUTPUT:
%   block_conc   - mb_size x mb_size matrix containing the concealed block
%

r = mb_size+1;
c = mb_size+1;
block_conc = -ones(mb_size);

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
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r-1,c+j-1-left) + d1/(d1 + d2)*support_area(r+mb_size, c+j-1+right);
                    
                
                %TO RIGHT    
                else
                    slope = tan(2*pi/360*angle);
                    aux = (mb_size - j + 1)*slope;
                    below = round(aux);
                    right = mb_size - j + 1;
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r-1,c+j-1-left) + d1/(d1 + d2)*support_area(r+i-1+below, c+mb_size);
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
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r+i-1-above,c-1) + d1/(d1 + d2)*support_area(r+mb_size, c+j-1+right);

                %TO RIGHT
                else
                    slope = tan(2*pi/360*angle);
                    aux = (mb_size - j + 1)*slope;
                    below = round(aux);
                    right = mb_size - j + 1;
                    d1 = sqrt(above^2 + left^2);
                    d2 = sqrt(below^2 + right^2);
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r+i-1-above,c-1) + d1/(d1 + d2)*support_area(r+i-1+below, c+mb_size);
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
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r-1,c+j-1+right) + d1/(d1 + d2)*support_area(r+mb_size, c+j-1-left);
                
                %TO LEFT
                else
                    slope = tan(2*pi/360*(90-(angle-90)));
                    aux = j*slope;
                    below = round(aux);
                    left = j;
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r-1,c+j-1+right) + d1/(d1 + d2)*support_area(r+i-1+below, c-1);
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
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r+i-1-above,c+mb_size) + d1/(d1 + d2)*support_area(r+mb_size, c+j-1-left);
                
                %TO LEFT
                else
                    slope = tan(2*pi/360*(90-(angle-90)));
                    aux = j*slope;
                    below = round(aux);
                    left = j;
                    d1 = sqrt(above^2 + right^2);
                    d2 = sqrt(below^2 + left^2);
                    block_conc(i,j) = d2/(d1 + d2)*support_area(r+i-1-above,c+mb_size) + d1/(d1 + d2)*support_area(r+i-1+below, c-1);
                end
                    
            end
            
        end
    end
end



end

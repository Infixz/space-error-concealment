%This function returns the most voted edge angle in the neighbourhood of
%the missing macroblock.

function theta = getAngle(neighbourhood_N, T, mb_size)

% INPUT:
%   neighbourhood_N - area around the missing macroblock (included), that also
%                   contains correctly received and decoded pixels
%   T               - threshold used in Eq.(10)
%   mb_size         - macroblock dimensions (in pixels)
%
% OUTPUT:
%   theta           - the angle (in degrees) of the edges with the
%                   strongest gradient, set to -1 if the gradient is not
%                   strong enough (monotone area, Eq.(10))
%

%We assume that the missing macroblock is located exactly at the center of
%the support area (which is squared). Thus, the top-left corner of the
%missing macroblock within the support area is located at
%[delimiter, delimiter].
delimiter = (length(neighbourhood_N(1,:)) - mb_size)/2;

D = zeros(1,8);

%Eight directional edge categories (Fig.4)
angle_range = 0:22.5:22.5*7;

% Sobel operators %%%%%%
Sy = [-1 -2 -1;
       0  0  0;
       1  2  1];
%NOTE: There is probably a typo in the paper (Eq.(6)) since Eq.(7) holds
%only if the first rows of Sy is negative and the last one positive, not
%viceversa as written in Eq.(6)
   
Sx = [-1 0 1;
      -2 0 2;
      -1 0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%

%Exploring the support area...
for i = 2:length(neighbourhood_N(1,:))-1
    for j = 2:length(neighbourhood_N(1,:))-1
        
        %We make sure no pixel in 3x3 area needed for Sobel operator is
        %missing
        if min(min(neighbourhood_N(i-1:i+1,j-1:j+1))) >= 0         
            
            %Applying the Sobel operator, Eq. (4) and Eq. (5)
            gx = sum(sum(Sx.*neighbourhood_N(i-1:i+1,j-1:j+1)));
            gy = sum(sum(Sy.*neighbourhood_N(i-1:i+1,j-1:j+1)));
                        
            slope = (gx/gy);
            offset = -(i) - ((j)*slope);
            
            %Checking whether the (prolongation of the) edge crosses the missing macroblock        
            for p = delimiter+1:delimiter+delimiter
                if (slope*p+offset >= -(delimiter+delimiter) && slope*p+offset <= -(delimiter+1))  || (slope>delimiter && j>delimiter && j<= delimiter+delimiter)
                    
                    %NOTE: the angle "theta" from Eq. (7) is the angle of
                    %the normal to the edge. In order to evaluate Eq.(8)
                    %we must compute the angles of the edges
                    %themselves.
                    theta = atan(slope)/pi*180;                    
                    if theta < 0
                        theta = theta + 180;
                    end
                    D(mod(round(theta/22.5)+8,8)+1) = D(mod(round(theta/22.5)+8,8)+1) + sqrt(gx^2 + gy^2);
                    
                    break
                end
            end
        end
    end
end

%Applying Eq.(9)
[value index] = max(D);

%Applying Eq.(10)
if value >= T
    %Edge area
    theta = angle_range(index);
else
    %Monotone area
    theta = -1;
end

end

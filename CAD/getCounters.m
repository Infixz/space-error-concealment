%This functions gets the eight counters "d" described in Section 3.

function d = getCounters(support_area, mb_size)

% INPUT:
%   support_area - area around the missing macroblock (included), that also
%                contains correctly received and decoded pixels
%   mb_size      - macroblock dimensions (in pixels)
%
% OUTPUT:
%   d           - gradient counters
%

d = zeros(1,8);

%We assume that the missing macroblock is located exactly at the center of
%the support area (which is squared). Thus, the top-left corner of the
%missing macroblock within the support area is located at
%[delimiter, delimiter].
delimiter = (length(support_area(1,:)) - mb_size)/2;

% Sobel operators %%%%%%
Sy = [-1 -2 -1;
       0  0  0;
       1  2  1];
   
Sx = [-1 0 1;
      -2 0 2;
      -1 0 1];
%%%%%%%%%%%%%%%%%%%%%%%%

%Exploring the support area...
for i = 2:length(support_area(1,:))-1
    for j = 2:length(support_area(1,:))-1
        
        %We make sure no pixel in 3x3 area needed for Sobel operator is
        %missing
        if min(min(support_area(i-1:i+1,j-1:j+1))) >= 0         
            
            %Applying the Sobel operator
            gx = sum(sum(Sx.*support_area(i-1:i+1,j-1:j+1)));
            gy = sum(sum(Sy.*support_area(i-1:i+1,j-1:j+1)));
            
            slope = (gx/gy);   
            offset = -(i) - ((j)*slope);             
            
            %Checking whether the (prolongation of the) edge crosses the missing macroblock        
            for p = delimiter+1:delimiter+delimiter
                if (slope*p+offset >= -(delimiter+delimiter) && slope*p+offset <= -(delimiter+1))  || (slope>delimiter && j>delimiter && j<= delimiter+delimiter)                                      
                    
                    theta = atan(slope)/pi*180;  
                    if theta < 0
                        theta = theta + 180;
                    end
                    
                    %Computing the counters "d" (Section 3)
                    d(mod(round(theta/22.5)+8,8)+1) = d(mod(round(theta/22.5)+8,8)+1) + sqrt(gx^2 + gy^2);                
                    break
                end
            end
        end
    end
end
  
end

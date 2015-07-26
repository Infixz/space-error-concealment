%This functions gets the eight counters "cm" described in Section 3.

function cm = getCounters(support_area, mb_size)

% INPUT:
%   support_area - area around the missing macroblock (included), that also
%                contains correctly received and decoded pixels
%   mb_size      - macroblock dimensions (in pixels)
%
% OUTPUT:
%   cm           - gradient counters
%

cm = zeros(1,8);

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

delta = 22.5/2;
angle_range = delta:22.5:180;

%Exploring the support area...
for i = 2:length(support_area(1,:))-1
    for j = 2:length(support_area(1,:))-1
        
        %We make sure no pixel in 3x3 area needed for Sobel operator is
        %missing
        if min(min(support_area(i-1:i+1,j-1:j+1))) >= 0         
            
            %Applying the Sobel operator, Eq. (9) and Eq. (10)
            %NOTE: in the paper, the aforementioned equations are erroneous
            %and do not correspond to the Sobel operator, although the
            %corresponding text says so. In this implementation, we apply
            %the correct operator (see "Sy" and "Sx").
            gx = sum(sum(Sx.*support_area(i-1:i+1,j-1:j+1)));
            gy = sum(sum(Sy.*support_area(i-1:i+1,j-1:j+1)));
            
            slope = (gx/gy);   
            offset = -(i) - ((j)*slope);             
            
            %Checking whether the (prolongation of the) edge crosses the missing macroblock        
            for p = delimiter+1:delimiter+delimiter
                if (slope*p+offset >= -(delimiter+delimiter) && slope*p+offset <= -(delimiter+1))  || (slope>delimiter && j>delimiter && j<= delimiter+delimiter)                                      
                    
                    %NOTE: the angle "theta" from Eq. (12) is the angle of
                    %the normal to the edge. In order to get the counters
                    %"cm" we must compute the angles of the edges
                    %themselves.
                    theta = atan(slope)/pi*180;  
                    if theta < 0
                        theta = theta + 180;
                    end
                    
                    %Computing the counters "cm" according to Fig. 1.
                    if theta  < angle_range(1) || theta > angle_range(8)
                        cm(8) = cm(8) + sqrt(gx^2 + gy^2);
                    else
                        for q = 2:8
                            if theta>angle_range(q-1) && theta <=angle_range(q)
                                cm(q-1) = cm(q-1) + sqrt(gx^2 + gy^2);
                                break
                            end
                        end
                    end                    
                    break
                end
            end
        end
    end
end
  
%Although it is not specified in the paper, in the unlikely case of totally
%homogeneous support area or if none of the computed edges crosses the lost
%macroblock, we set the counters "cm" to 1. If we did not do that, the
%counters would remain set to 0 which leads to a numerical instability
%(0/0, see Eq. (14)).
if sum(cm) == 0
    cm = ones(size(cm));
end

end

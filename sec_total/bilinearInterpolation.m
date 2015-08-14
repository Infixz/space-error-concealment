%This function performs the concealment for one macroblock according to
%Eq.(1)

function y = bilinearInterpolation(support_area, lambda)

% INPUT:
%   support_area - contains all the available pixels (1 pixel thick layer)
%                around the missing macroblock
%   lambda       - parameter that determines the combination of vertical
%                and horizontal interpolations, see Eq.(1)
% OUTPUT
%   y            - concealed macroblock
%

mb_size = (length(support_area(1,:))-2);
y = zeros(mb_size);

for i = 1+1:length(support_area(1,:))-1
    if support_area(1,i) < 0
        support_area(1,i) = support_area(end,i);
    end
    if support_area(end,i) < 0
        support_area(end,i) = support_area(1,i);
    end
    if support_area(i,1) < 0
        support_area(i,1) = support_area(i,end);
    end
    if support_area(i,end) < 0
        support_area(i,end) = support_area(i,1);
    end
end


for i = 1:mb_size
    for j = 1:mb_size
        
        mu1 = j/(mb_size+1);    %horizontal distance
        mu2 = i/(mb_size+1);    %vertical distance
        
        y(i,j) = lambda*(mu1*support_area(i+1,end) + (1-mu1)*support_area(i+1,1)) + (1-lambda)*(mu2*support_area(end,j+1) + (1-mu2)*support_area(1,j+1));
        
    end
end
        



end

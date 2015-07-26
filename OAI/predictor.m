%Function that estimates the value of Xn, Eq.(2)

function Xn = predictor(support_area, r, c, T)

% INPUT:
%   support_area - area around the missing macroblock (included), that also
%                contains correctly received and decoded pixels
%   [r, c]       - coordenates (row and column) within the support_area of
%                the pixel to be concealed
%   T            - T parameter, see Fig. 2
%
% OUTPUT:
%   Xn           - estimated pixel value
%

[rows cols] = size(support_area);

%Allocating memory for "y" (Mx1) and "C" (MxN). Note that these two
%variables will need to be resized afterwards as Nn > Nn_valid
%and Mn > Mn_valid
y = zeros((2*T+1)*(2*T+1)-1,1);
C = zeros((2*T+1)*(2*T+1)-1,8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mask_Nn_valid = zeros(3);
Nn_valid = zeros(1,8);
counter_Nn_valid = 1;
for i = -1:1
    for j = -1:1
        %Checking validity
        if r+i>0 && c+j>0 && r+i<=rows && c+j<=cols && support_area(r+i,c+j) >= 0
            mask_Nn_valid(i+2,j+2) = 1;
            Nn_valid(counter_Nn_valid) = support_area(r+i,c+j);
            counter_Nn_valid = counter_Nn_valid + 1;
        end
    end
end
Nn_valid = Nn_valid(1:counter_Nn_valid-1);


counterY = 1;
for i = -T:T
    for j = -T:T
        %For every pixel in the valid set "Mn"...
        if r+i>0 && c+j>0 && r+i<=rows && c+j<=cols && support_area(r+i,c+j) >= 0
            
            %...we extract the valid set "Nn" (the available pixels from the
            %8 closest ones) and update the matrix "C"
            counterC = 1;
            for p = -1:1
                for q = -1:1
                    if mask_Nn_valid(p+2,q+2) > 0 && support_area(r+i+p,c+j+q) >= 0
                        C(counterY,counterC) = support_area(r+i+p,c+j+q);
                        counterC = counterC + 1;
                    end
                end
            end
            
            %Updating "y"
            if counterC-1 == sum(mask_Nn_valid(:))
                y(counterY) = support_area(r+i,c+j);
                counterY = counterY + 1;           
            end
        end   
    end
end

%Dropping the "overallocated" elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = y(1:counterY-1);
C = C(1:counterY-1,1:sum(mask_Nn_valid(:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Applying Eq.(5)
a = (C'*C)\(C'*y);

%Applying Eq.(2)
Xn = Nn_valid*a;


%Saturating the output, if necessary (although it is not specified in the
%paper)
if Xn < 0
    Xn = 0;
elseif Xn > 255
    Xn = 255;
    
%In same cases, the vector "a" (Eq.(5)) cannot be computed since the
%matrices are badly scaled. Although this situation is not contemplated in
%the paper, if it does happen, we set the pixel to the mean value of its
%valid set Nn.
elseif isnan(Xn)
    Xn = mean(Nn_valid);    
end

end

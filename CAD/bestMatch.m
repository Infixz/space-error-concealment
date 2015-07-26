%This function finds the best matching macroblock by minimizing the squared
%error, see Eq.(4)

function block_conc = bestMatch(center, mask, mb_size)

% INPUT:
%   center     - 2x1 vector containing the row and the column of the top-left
%              pixel of the missing macroblock within the image
%   mask       - received (lossy) frame
%   mb_size    - macroblock dimensions (in pixels)
%
% OUTPUT:
%   block_conc - mb_size x mb_size matrix containing the concealed block
%

r = center(1);
c = center(2);

%Getting the outer layer of the missing macroblock %%%%%%%%%%%%%%%%%%%%%%%%
neighbourhood = zeros(1,4*mb_size + 4);
neighbourhood(1:mb_size+2) = mask(r-1,c-1:c+mb_size);
neighbourhood(mb_size+3:mb_size+3+mb_size) = mask(r:r+mb_size,c+mb_size);
neighbourhood(mb_size+mb_size+4:mb_size+mb_size+4+mb_size) = mask(r+mb_size,c+mb_size-1:-1:c-1);
neighbourhood(mb_size+mb_size+mb_size+5:mb_size+mb_size+mb_size+5+mb_size-1) = mask(r+mb_size-1:-1:r,c-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Search area as specified in Section 4, see also Eq.(4) %%%%%%%%%%%%%%%%%%%
search_area = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

block_conc = zeros(mb_size);
[rows cols] = size(mask);
TH = 10^6;

%Exploring the search area
for i = max(2,r-search_area*mb_size):min(rows-mb_size,r+mb_size+search_area*mb_size)
    for j = max(2,c-search_area*mb_size):min(cols-mb_size,c+mb_size+search_area*mb_size)
        
        %If the candidate macroblock is complete...
        mb = mask(i:i+mb_size-1,j:j+mb_size-1);
        if min(mb(:)) < 0
            continue
        
        %...then get its outer layer...    
        else            
            candidate = zeros(1,4*mb_size + 4);
            candidate(1:mb_size+2) = mask(i-1,j-1:j+mb_size);
            candidate(mb_size+3:mb_size+3+mb_size) = mask(i:i+mb_size,j+mb_size);
            candidate(mb_size+mb_size+4:mb_size+mb_size+4+mb_size) = mask(i+mb_size,j+mb_size-1:-1:j-1);
            candidate(mb_size+mb_size+mb_size+5:mb_size+mb_size+mb_size+5+mb_size-1) = mask(i+mb_size-1:-1:i,j-1);
            
            %...and compute the squared error
            aux = (candidate >= 0);
            sq_error = sum((((neighbourhood - candidate).^2).*aux))/sum(aux);
            if sq_error < TH
                block_conc = mb;
                TH = sq_error;
            end
        end
    end
end   

end

% This function performs the content estimation as described in Section 3.

function [content angSelector] = decideMode(d)

% INPUT:
%   d           - gradient counters, see Section 3.
% 
% OUTPUT:
%   content     - content estimate
%               = 1 -> smooth
%               = 2 -> edges
%               = 3 -> texture
%   angSelector - 1x8 boolean vector, set to zero unless the content
%               estimate is "edges", in that case is set to 1 for every
%               relevant direction (edge)
%

angSelector = zeros(1,8);

%Thresholds as specified in Section 5 %%%%%%%%%%%%%%%%%
T = 0.55;
Td = 3000;
TN = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Filtering weak edges and getting the strongest edge %%
d = d.*(d>T*max(d));
dmax = max(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth content
if dmax < Td
    content = 1;

%edges
elseif dmax >= Td && sum(d > 0) <= TN
    content = 2;
    angSelector = d;
    
%texture
else
    content = 3;
end

end

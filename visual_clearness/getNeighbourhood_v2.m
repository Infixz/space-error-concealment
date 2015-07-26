function [set1 set2 set3] = getNeighbourhood_v2(mb, linea, theta, mb_size, side)

%tres conjuntos: 1. linea, 2. izda/top y 3. dcha/bottom, sin ordenar

mb = double(mb);
[dim dimm] = size(mb);
set1 = -ones(1,dim*dimm);
set2 = -ones(1,dim*dimm);
set3 = -ones(1,dim*dimm);

thickness = 3;
counter2 = 1;
counter3 = 1;

for k = 1:length(linea(1,:))
    set1(k) = mb(linea(1,k),linea(2,k));
    mb(linea(1,k),linea(2,k)) = -1;
end

%TOP AND BOTTOM
if theta < 45 || theta > 135
    c = 0;
    r = 0;
    for i = 1:length(linea(1,:))
        if c ~= linea(2,i) || r ~= linea(1,i)
            r = linea(1,i);
            c = linea(2,i);
            %top
            iteratorTop = 1;
            while true
                if r-iteratorTop > 0 && iteratorTop <= thickness && mb(r-iteratorTop,c) >= 0
                    set2(counter2) = mb(r-iteratorTop,c);
                    counter2 = counter2 + 1;
                    iteratorTop = iteratorTop + 1;
                else
                    break;
                end
            end
            %bottom
            iteratorBottom = 1;
            while true
                if r+iteratorBottom <= mb_size && iteratorBottom <= thickness && mb(r+iteratorBottom,c) >= 0
                    set3(counter3) = mb(r+iteratorBottom,c);
                    counter3 = counter3 + 1;
                    iteratorBottom = iteratorBottom + 1;
                else
                    break;
                end
            end          
        end
    end
    
    
%LEFT AND RIGHT    
else
    r = 0;
    c = 0;
    for i = 1:length(linea(1,:))
        if c ~= linea(2,i) || r ~= linea(1,i)
            r = linea(1,i);
            c = linea(2,i);
            %left
            iteratorLeft = 1;
            while true
                if c-iteratorLeft > 0 && iteratorLeft <= thickness && mb(r,c-iteratorLeft) >= 0
                    set2(counter2) = mb(r,c-iteratorLeft);
                    counter2 = counter2 + 1;
                    iteratorLeft = iteratorLeft + 1;
                else
                    break;
                end
            end
            %right
            iteratorRight = 1;
            while true
                if c+iteratorRight <= mb_size && iteratorRight <= thickness && mb(r,c+iteratorRight) >= 0
                    set3(counter3) = mb(r,c+iteratorRight);
                    counter3 = counter3 + 1;
                    iteratorRight = iteratorRight + 1;
                else
                    break;
                end
            end          
        end
    end 
end
   
set1 = set1(1:length(linea(1,:)));
set2 = set2(1:counter2-1);
set3 = set3(1:counter3-1);

%FOR BOTTOM & TOP
if strcmp(side, 'BT') && theta >= 45
    aux2 = set2;
    aux3 = set3;
    clear set2
    clear set3
    set2 = aux3;
    set3 = aux2;
end

%FOR LEFT & RIGHT
if strcmp(side, 'LR') && theta >= 135
    aux2 = set2;
    aux3 = set3;
    clear set2
    clear set3
    set2 = aux3;
    set3 = aux2;
end
    
end

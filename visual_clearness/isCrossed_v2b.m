function [yOut yEnter xEnter yExit xExit] = isCrossed_v2b(rf, cf, ri, ci, mb_size)

% ri - row initial
% ci - column initial

% rf - row final
% cf - column final

%mb_size

x1 = ci;
y1 = -ri;
x2 = cf;
y2 = -rf;

mat = [x1 1 y1; x2 1 y2];
mat = rref(mat);
slope = mat(1,end);
offset = mat(2,end);

yOut = false;
xEnter = 0;
yEnter = 0;
xExit = 0;
yExit = 0;
TH = 1;

x = mb_size:0.01:(mb_size + mb_size - 1);
for i = 1:length(x)
    y = slope*x(i) + offset;
    if -y >= mb_size+TH && -y <= mb_size + mb_size-TH
        yOut = true;
        break;
    end
end

%linea vertical
aux = abs(ri - rf)/(abs(ci - cf) + 0.00001);
if aux > mb_size && (ci > mb_size && ci <= 2*mb_size) && (cf > mb_size && cf <= 2*mb_size)
    yOut = true;
    yEnter = mb_size + 1;
    yExit = mb_size*2;
    xEnter = cf;
    xExit = cf;
    return
elseif aux > mb_size
    yOut = false;
    return
end
    
%linea horizontal
if abs(slope) < 1/mb_size && (ri > mb_size && ri <= 2*mb_size) && (rf > mb_size && rf <= 2*mb_size)
    yOut = true;
    xEnter = mb_size + 1;
    xExit = 2*mb_size;
    yEnter = ri;
    yExit = rf;
    return
elseif slope == 0
    yOut = false;
    return
end

%if crossed
if yOut
    counter = 1;
    filas = zeros(1,2);
    columnas = zeros(1,2);
    x = -2*mb_size:0.01:2*mb_size;
      
    %aaa = zeros(1,length(x));
    %bbb = zeros(1,length(x));
    
    for i = 1:length(x)
        if counter == 3
            break
        end
        dX = x(i);
        dY = slope*x(i);
        rAux = round(-rf + dY);
        cAux = round(cf + dX);
        
        %aaa(i) = -rAux;
        %bbb(i) = cAux;
        
        if ((-rAux) == mb_size+1 || (-rAux) == 2*mb_size) && (cAux <= 2*mb_size && cAux > mb_size)
            if counter > 1
                if ((-rAux) ~= filas(counter-1) && cAux ~= columnas(counter-1))
                    filas(counter) = -rAux;
                    columnas(counter) = cAux;
                    counter = counter + 1;
                end
            else
                filas(counter) = -rAux;
                columnas(counter) = cAux;
                counter = counter + 1;
            end        
            
        elseif (cAux == mb_size+1 || cAux == 2*mb_size) && ((-rAux) <= 2*mb_size && (-rAux) > mb_size)
            if counter > 1
                if ((-rAux) ~= filas(counter-1) && cAux ~= columnas(counter-1))
                    filas(counter) = -rAux;
                    columnas(counter) = cAux;
                    counter = counter + 1;
                end
            else
                filas(counter) = -rAux;
                columnas(counter) = cAux;
                counter = counter + 1;
            end
        end
    end
else
    yOut = false;
    return
end

if norm([rf-filas(1) cf-columnas(1)]) < norm([rf-filas(2) cf-columnas(2)]) 
    xEnter = columnas(1);
    yEnter = filas(1);
    xExit = columnas(2);
    yExit = filas(2);
else
    xEnter = columnas(2);
    yEnter = filas(2);
    xExit = columnas(1);
    yExit = filas(1);
end

if 0 & columnas(1) == 0 && columnas(2) == 0 && filas(1) == 0 && filas(2) == 0
    slope
    ri
    ci
    rf
    cf
end

%plot(x,aaa,'r',x,bbb,'b')
    
end








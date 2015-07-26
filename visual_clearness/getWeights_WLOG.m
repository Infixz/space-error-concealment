function y = getWeights_WLOG(ang, puntos, mb_size, N)

%devuelve pesos, si la interpolacion no cuenta, peso puesto a cero
%puntos = [xEnter1 yEnter1 xExit1 yExit1; ...]

dim = length(ang);
y = ones(mb_size, mb_size, dim);
if sum(ang) == -dim
    return
end

flags = zeros(1,dim);
slope = zeros(1,dim);
offset = zeros(1,dim);
epsilon = 0.00001;

for i = 1:N
    if ang(i) < 0
        continue
    end
    %linea vertical
    aux = abs(puntos(i,2) - puntos(i,4))/(abs(puntos(i,1) - puntos(i,3)) + epsilon);
    if aux > mb_size && (puntos(i,1) > mb_size && puntos(i,1) <= 2*mb_size) && (puntos(i,3) > mb_size && puntos(i,3) <= 2*mb_size)
        flags(i) = true;
    else
        %slope & offset
        mat = [puntos(i,1) 1 -puntos(i,2); puntos(i,3) 1 -puntos(i,4)];
        mat = rref(mat);
        slope(i) = mat(1,end);
        offset(i) = mat(2,end);
    end
end


%para cada pixel del MB
for w = 1:N
    if ang(w) < 0
        continue
    else
        for i = mb_size+1:mb_size+mb_size
            for j = mb_size+1:mb_size+mb_size

                if flags(w)
                    y(j-mb_size,i-mb_size,w) = abs(i - puntos(w,1)) + epsilon;
                else
                    y(j-mb_size,i-mb_size,w) = abs(slope(w)*i + j + offset(w))/sqrt(slope(w)^2 + 1) + epsilon;
                end
            end
        end
    end
end
      
%normalizamos
for i = mb_size+1:mb_size+mb_size
    for j = mb_size+1:mb_size+mb_size
        total = sum(y(j-mb_size,i-mb_size,1:N));
        for w = 1:dim
            if ang(w) < 0
                continue
            else
                y(j-mb_size,i-mb_size,w) = y(j-mb_size,i-mb_size,w)/total;
            end
        end
    end
end
    
    

end
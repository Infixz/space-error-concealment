function y = variance(linea, block, mb_size, theta)

if 0
gx = 0;
gy = 0;

totalX = 0;
totalY = 0;

for i = 1:length(linea(1,:))
    r = linea(1,i);
    c = linea(2,i);
    m = block(r,c);
    
    %ROW GRAD
    if r > 1 && r < mb_size
        
        gx = gx + (abs(m - block(r-1,c)) + abs(m - block(r+1,c)))/2;
        totalX = totalX + 1;
    end
    
    %COL GRAD
    if c > 1 && c < mb_size
        gy = gy + (abs(m - block(r,c-1)) + abs(m - block(r,c+1)))/2;
        totalY = totalY + 1;
    end
end

gx = gx/totalX;
gy = gy/totalY;
y = sqrt(gx^2 +gy^2);



elseif 0
[side1 side2] = getDescriptors(block, linea, theta, mb_size);
total = 0;
y = 0;
for i = 1:length(side1(1,:))
    if side1(1,i) >= 0 && side2(1,i) >= 0
        y = y + (side1(1,i)-side2(1,i));
        total = total + 1;
    end
end

y = abs(y)/total;




else
mb = block;
for k = 1:length(linea(1,:))
    mb(linea(1,k),linea(2,k)) = -1;
end
LEN_MAX = 3;

auxY = -ones(size(linea));
auxX = -ones(size(linea));

for i = 1:length(linea(1,:))
    %Y GRAD
    c = linea(2,i);
    for k = 1:LEN_MAX
        r = linea(1,i)+k;
        if r <= mb_size && mb(r,c) >= 0
            auxY(1,i) = mb(r,c);
            break
        end
    end
    for k = 1:LEN_MAX
        r = linea(1,i)-k;
        if r > 0 && mb(r,c) >= 0
            auxY(2,i) = mb(r,c);
            break
        end
    end
    
    %X GRAD
    r = linea(1,i);
    for k = 1:LEN_MAX
        c = linea(2,i)+k;
        if c <= mb_size && mb(r,c) >= 0
            auxX(1,i) = mb(r,c);
            break
        end
    end
    for k = 1:LEN_MAX
        c = linea(2,i)-k;
        if c > 0 && mb(r,c) >= 0
            auxX(2,i) = mb(r,c);
            break
        end
    end
end

for i = 1:length(auxY(1,:))
    if auxY(1,i) < 0 || auxY(2,i) < 0
        auxY(1,i) = 0;
        auxY(2,i) = 0;
    end
    if auxX(1,i) < 0 || auxX(2,i) < 0
        auxX(1,i) = 0;
        auxX(2,i) = 0;
    end
end

yGrad = auxY(1,:) - auxY(2,:);
xGrad = auxX(1,:) - auxX(2,:);

y = mean(sqrt(xGrad.^2 + yGrad.^2));

end

end



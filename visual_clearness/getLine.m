function yOut = getLine(x1, y1, x2, y2, mb_size)

selector = 0;   %0 => from point to point
                %1 => across the whole MB, passing through the points

if selector
    equ = [x1 1 y1; x2 1 y2];
    equ = rref(equ);
    slope = equ(1,end);
    offset = equ(2,end);
    %beginning
    x1 = 1;
    y1 = round(slope*x1 + offset);
    if y1 < 1
        y1 = 1;
        x1 = round((y1 - offset)/slope);
    end
    %end
    x2 = mb_size;
    y2 = round(slope*x2 + offset);
    if y2 > mb_size
        y2 = mb_size;
        x2 = round((y2 - offset)/slope);
    end
end
                
                
len = 2*mb_size;
x = linspace(x1,x2, len);
y = linspace(y1,y2, len);

yOut = -ones(2,len);
yOut(:,1) = [x1; y1];
counter = 1;

for i = 2:len
    if round(x(i)) ~= yOut(1,counter) || round(y(i)) ~= yOut(2,counter)
        counter = counter + 1;
        yOut(:,counter) = [round(x(i)); round(y(i))];
    end
end

[v i] = min(yOut(1,:));
yOut = yOut(:,1:i-1);

end

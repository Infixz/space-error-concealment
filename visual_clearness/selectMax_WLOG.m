function [ang w puntos] = selectMax_WLOG(left, right, top, bottom)

%X = [theta; H; diff; mu2; mu3; var2; var3; xEnter; yEnter; xExit; yExit];
%        1   2    3    4    5    6      7      8      9       10     11

TH = 3.5;
dim = length(left(1,:));
puntos = zeros(4,dim,4);
puntos(:,:,1) = [left(8,:); right(8,:); top(8,:); bottom(8,:)];
puntos(:,:,2) = [left(9,:); right(9,:); top(9,:); bottom(9,:)];
puntos(:,:,3) = [left(10,:); right(10,:); top(10,:); bottom(10,:)];
puntos(:,:,4) = [left(11,:); right(11,:); top(11,:); bottom(11,:)];

index = zeros(1,4);
sigma = [left(2,:).*left(3,:); right(2,:).*right(3,:); top(2,:).*top(3,:); bottom(2,:).*bottom(3,:)];

% POINT 1
[vL index(1)] = max(sigma(1,:));
[vR index(2)] = max(sigma(2,:));
[vT index(3)] = max(sigma(3,:));
[vB index(4)] = max(sigma(4,:));

angVec = [left(1,index(1)) right(1,index(2)) top(1,index(3)) bottom(1,index(4))];
strVec = [vL vR vT vB];
[w1 i] = max(strVec);
ang1 = angVec(i);
punto1 = puntos(i,index(i),:);

for l = 1:4
    for k = 1:dim
        xe = puntos(l,k,1);
        ye = puntos(l,k,2);
        xx = puntos(l,k,3);
        yx = puntos(l,k,4);
        dd = [sqrt((xe-punto1(1))^2 +  (ye-punto1(2))^2) sqrt((xx-punto1(3))^2 +  (yx-punto1(4))^2)];
        if dd(1) <= TH && dd(2) <= TH
            sigma(l,k) = -1;
        end
    end
end

%POINT 2
[vL index(1)] = max(sigma(1,:));
[vR index(2)] = max(sigma(2,:));
[vT index(3)] = max(sigma(3,:));
[vB index(4)] = max(sigma(4,:));

angVec = [left(1,index(1)) right(1,index(2)) top(1,index(3)) bottom(1,index(4))];
strVec = [vL vR vT vB];
[w2 i] = max(strVec);
ang2 = angVec(i);
punto2 = puntos(i,index(i),:);

for l = 1:4
    for k = 1:dim
        xe = puntos(l,k,1);
        ye = puntos(l,k,2);
        xx = puntos(l,k,3);
        yx = puntos(l,k,4);
        dd = [sqrt((xe-punto2(1))^2 +  (ye-punto2(2))^2) sqrt((xx-punto2(3))^2 +  (yx-punto2(4))^2)];
        if dd(1) <= TH && dd(2) <= TH
            sigma(l,k) = -1;
        end
    end
end

%POINT 3
[vL index(1)] = max(sigma(1,:));
[vR index(2)] = max(sigma(2,:));
[vT index(3)] = max(sigma(3,:));
[vB index(4)] = max(sigma(4,:));

angVec = [left(1,index(1)) right(1,index(2)) top(1,index(3)) bottom(1,index(4))];
strVec = [vL vR vT vB];
[w3 i] = max(strVec);
ang3 = angVec(i);
punto3 = puntos(i,index(i),:);

for l = 1:4
    for k = 1:dim
        xe = puntos(l,k,1);
        ye = puntos(l,k,2);
        xx = puntos(l,k,3);
        yx = puntos(l,k,4);
        dd = [sqrt((xe-punto3(1))^2 +  (ye-punto3(2))^2) sqrt((xx-punto3(3))^2 +  (yx-punto3(4))^2)];
        if dd(1) <= TH && dd(2) <= TH
            sigma(l,k) = -1;
        end
    end
end


%POINT 4
[vL index(1)] = max(sigma(1,:));
[vR index(2)] = max(sigma(2,:));
[vT index(3)] = max(sigma(3,:));
[vB index(4)] = max(sigma(4,:));

angVec = [left(1,index(1)) right(1,index(2)) top(1,index(3)) bottom(1,index(4))];
strVec = [vL vR vT vB];
[w4 i] = max(strVec);
ang4 = angVec(i);
punto4 = puntos(i,index(i),:);

for l = 1:4
    for k = 1:dim
        xe = puntos(l,k,1);
        ye = puntos(l,k,2);
        xx = puntos(l,k,3);
        yx = puntos(l,k,4);
        dd = [sqrt((xe-punto4(1))^2 +  (ye-punto4(2))^2) sqrt((xx-punto4(3))^2 +  (yx-punto4(4))^2)];
        if dd(1) <= TH && dd(2) <= TH
            sigma(l,k) = -1;
        end
    end
end


%POINT 5
[vL index(1)] = max(sigma(1,:));
[vR index(2)] = max(sigma(2,:));
[vT index(3)] = max(sigma(3,:));
[vB index(4)] = max(sigma(4,:));

angVec = [left(1,index(1)) right(1,index(2)) top(1,index(3)) bottom(1,index(4))];
strVec = [vL vR vT vB];
[w5 i] = max(strVec);
ang5 = angVec(i);
punto5 = puntos(i,index(i),:);

for l = 1:4
    for k = 1:dim
        xe = puntos(l,k,1);
        ye = puntos(l,k,2);
        xx = puntos(l,k,3);
        yx = puntos(l,k,4);
        dd = [sqrt((xe-punto5(1))^2 +  (ye-punto5(2))^2) sqrt((xx-punto5(3))^2 +  (yx-punto5(4))^2)];
        if dd(1) <= TH && dd(2) <= TH
            sigma(l,k) = -1;
        end
    end
end


if w2 < 0.5*w1
    ang2 = -1;
end
if w3 < 0.5*w1
    ang3 = -1;
end
if w4 < 0.5*w1
    ang4 = -1;
end
if w5 < 0.5*w1
    ang5 = -1;
end

ang = [ang1 ang2 ang3 ang4 ang5];
w = [w1 w2 w3 w4 w5];
puntos = [  punto1(1) punto1(2) punto1(3) punto1(4);
            punto2(1) punto2(2) punto2(3) punto2(4);
            punto3(1) punto3(2) punto3(3) punto3(4);
            punto4(1) punto4(2) punto4(3) punto4(4);
            punto5(1) punto5(2) punto5(3) punto5(4)];


end






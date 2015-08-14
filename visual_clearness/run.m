%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2011 by                                          %
%   Jan Koloda                                                     %
%   Universidad de Granada, Granada, Spain                         %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   Koloda, J., Sanchez, V. and Peinado, A.M., "Spatial Error      % 
%   Concealment Based on Edge Visual Clearness for Image/Video     %
%   Communication", Circuits, Systems and Signal Processing        %
%                                                                  %
%   This program is free of charge for personal and scientific     %
%   use (with proper citation). The author does NOT give up his    %
%   copyright. Any commercial use is prohibited.                   %
%   YOU ARE USING THIS PROGRAM AT YOUR OWN RISK! THE AUTHOR        %
%   IS NOT RESPONSIBLE FOR ANY DAMAGE OR DATA-LOSS CAUSED BY THE   %
%   USE OF THIS PROGRAM.                                           %
%                                                                  %
%   If you have any questions please contact:                      %
%                                                                  %
%   Jan Koloda                                                     %
%   Dpt. Signal Theory, Networking and Communication               %
%   Universidad de Granada                                         %
%   Granada, Spain                                                 %
%                                                                  %
%                                                                  %
%   email: janko @ ugr.es                                          %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% init
clear;
clc;

%% Loading the image and parsing it to grey level, if needed
[filename,pathname]=uigetfile({'*.tiff';'*.tif';'*.*';'*.bmp'},'图像文件');
file=[pathname,filename];
img = imread(file);
[r c d] = size(img);
if d > 1
    img = double(rgb2gray(img));
else
    img = double(img);
end

%% Generating losses
tic    %开始计时

slice_to_be_lost = 1;
nSlices = 2;
mb_size = 16;
mode = 'default';

%Cropping the image so it is made of an integer number of macroblocks %%%%%
img = img(1:floor(r/mb_size)*mb_size,1:floor(c/mb_size)*mb_size);
[r c] = size(img);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mask centers] = simuLoss(img, mb_size, nSlices, slice_to_be_lost, mode);
received_frame = mask;

%% Concealment

STEP_SIZE = 1;


N = 5;
img = double(img);
resRho = 0.9;
resThe = 2;
minLen = 1;
BORDER = 4;
 

maskN5 = mask;
maskN2 = mask;
maskN3 = mask;
maskN4 = mask;

len = length(centers(1,:));
[dim1 dim2] = size(mask);


h = waitbar(0);
set(h,'Name','Processing')

for w = 1:len
    r = centers(1,w);
    c = centers(2,w);    

    if c > 1 && r > 1 && r < dim1 - mb_size && c < dim2 - mb_size

        varL = -ones(1,2*mb_size+1);
        varR = -ones(1,2*mb_size+1);
        varT = -ones(1,2*mb_size+1);
        varB = -ones(1,2*mb_size+1);

        left = -ones(11,2*mb_size+1);
        right = -ones(11,2*mb_size+1);
        top = -ones(11,2*mb_size+1);
        bottom = -ones(11,2*mb_size+1);

        leftPts = -ones(4,2*mb_size+1);
        rightPts = -ones(4,2*mb_size+1);
        topPts = -ones(4,2*mb_size+1);
        bottomPts = -ones(4,2*mb_size+1);

        for i = -mb_size:STEP_SIZE:mb_size
            blockL = img(r+i:r+i+mb_size-1,c-mb_size:c-mb_size+mb_size-1);   %left
            blockR = img(r+i:r+i+mb_size-1,c+mb_size:c+mb_size+mb_size-1);   %right
            blockT = img(r-mb_size:r-mb_size+mb_size-1,c+i:c+i+mb_size-1);   %top
            blockB = img(r+mb_size:r+mb_size+mb_size-1,c+i:c+i+mb_size-1);    %bottom

            th = [0.4 0.7];
            mbL = edge(blockL, 'canny',th);
            mbR = edge(blockR, 'canny',th);
            mbT = edge(blockT, 'canny',th);
            mbB = edge(blockB, 'canny',th);

            %TOP
            [H,T,R] = hough(mbT, 'RhoResolution', resRho, 'ThetaResolution', resThe);         
            P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
            lines = houghlines(mbT,T,R,P,'FillGap',20,'MinLength',minLen);
            [check d2] = size(fieldnames(lines));
            if check > 0
                [d1 d2] = size(lines);
                for k = 1:d2
                    ri = lines(k).point1(2);
                    ci = lines(k).point1(1);
                    rf = lines(k).point2(2);
                    cf = lines(k).point2(1);
                    col = P(k,2);
                    [yOut yEnter xEnter yExit xExit] = isCrossed_v2b(rf, cf+(i+mb_size), ri, ci+(i+mb_size), mb_size);
                    yOut = yOut && (ri > mb_size - BORDER || rf > mb_size - BORDER);
                    if yOut
                        topPts(:,i+mb_size+1) = [xEnter; yEnter; xExit; yExit];
                        linea = getLine(ri, ci, rf, cf, mb_size);
                        varT(1,i+mb_size+1) = variance(linea, blockT, mb_size, T(col)+90);
                        break;
                    end
                end
            end
            if check > 0 && yOut
                linea = getLine(ri, ci, rf, cf, mb_size);
                [set1 set2 set3] = getNeighbourhood_v2(blockT, linea, T(col)+90, mb_size, 'BT');
                top(:,i+mb_size+1) = getDescriptors(linea, blockT, mb_size, T(col)+90, set2, set3, H(P(k,1),P(k,2)), xEnter, yEnter, xExit, yExit);
            end

            %BOTTOM
            [H,T,R] = hough(mbB, 'RhoResolution', resRho, 'ThetaResolution', resThe);         
            P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
            lines = houghlines(mbB,T,R,P,'FillGap',20,'MinLength',minLen);
            [check d2] = size(fieldnames(lines));
            if check > 0
                [d1 d2] = size(lines);
                for k = 1:d1
                    ri = lines(k).point1(2);
                    ci = lines(k).point1(1);
                    rf = lines(k).point2(2);
                    cf = lines(k).point2(1);
                    col = P(k,2);
                    [yOut yEnter xEnter yExit xExit] = isCrossed_v2b(rf+mb_size+mb_size, cf+(i+mb_size), ri+mb_size+mb_size, ci+(i+mb_size), mb_size);
                    yOut = yOut && (ri < BORDER || rf < BORDER);
                    if yOut
                        bottomPts(:,i+mb_size+1) = [xEnter; yEnter; xExit; yExit];
                        linea = getLine(ri, ci, rf, cf, mb_size);
                        varB(1,i+mb_size+1) = variance(linea, blockB, mb_size, T(col)+90);
                        break;
                    end
                end
            end
            if check > 0 && yOut
                linea = getLine(ri, ci, rf, cf, mb_size);
                [set1 set2 set3] = getNeighbourhood_v2(blockB, linea, T(col)+90, mb_size, 'BT');
                bottom(:,i+mb_size+1) = getDescriptors(linea, blockB, mb_size, T(col)+90, set2, set3, H(P(k,1),P(k,2)), xEnter, yEnter, xExit, yExit);
            end

            %LEFT
            [H,T,R] = hough(mbL, 'RhoResolution', resRho, 'ThetaResolution', resThe);         
            P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
            lines = houghlines(mbL,T,R,P,'FillGap',20,'MinLength',minLen);
            [check d2] = size(fieldnames(lines));
            if check > 0
                [d1 d2] = size(lines);
                for k = 1:d1
                    ri = lines(k).point1(2);
                    ci = lines(k).point1(1);
                    rf = lines(k).point2(2);
                    cf = lines(k).point2(1);
                    col = P(k,2);
                    [yOut yEnter xEnter yExit xExit] = isCrossed_v2b(rf+(i+mb_size), cf, ri+(i+mb_size), ci, mb_size);
                    yOut = yOut && (ci > mb_size - BORDER || cf > mb_size - BORDER);
                    if yOut
                        leftPts(:,i+mb_size+1) = [xEnter; yEnter; xExit; yExit];
                        linea = getLine(ri, ci, rf, cf, mb_size);
                        varL(1,i+mb_size+1) = variance(linea, blockL, mb_size, T(col)+90);
                        break;
                    end
                end
            end
            if check > 0 && yOut
                linea = getLine(ri, ci, rf, cf, mb_size);
                [set1 set2 set3] = getNeighbourhood_v2(blockL, linea, T(col)+90, mb_size, 'LR');
                left(:,i+mb_size+1) = getDescriptors(linea, blockL, mb_size, T(col)+90, set2, set3, H(P(k,1),P(k,2)), xEnter, yEnter, xExit, yExit);
            end

            %RIGHT
            [H,T,R] = hough(mbR, 'RhoResolution', resRho, 'ThetaResolution', resThe);         
            P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
            lines = houghlines(mbR,T,R,P,'FillGap',20,'MinLength',minLen);
            [check d2] = size(fieldnames(lines));
            if check > 0
                [d1 d2] = size(lines);
                for k = 1:d1
                    ri = lines(k).point1(2);
                    ci = lines(k).point1(1);
                    rf = lines(k).point2(2);
                    cf = lines(k).point2(1);
                    col = P(k,2);
                    [yOut yEnter xEnter yExit xExit] = isCrossed_v2b(rf+(i+mb_size), cf+mb_size+mb_size, ri+(i+mb_size), ci+mb_size+mb_size, mb_size);
                    yOut = yOut && (ci < BORDER || cf < BORDER);
                    if yOut
                        rightPts(:,i+mb_size+1) = [xEnter; yEnter; xExit; yExit];
                        linea = getLine(ri, ci, rf, cf, mb_size);
                        varR(1,i+mb_size+1) = variance(linea, blockR, mb_size, T(col)+90);
                        break;
                    end
                end
            end
            if check > 0 && yOut
                linea = getLine(ri, ci, rf, cf, mb_size);
                [set1 set2 set3] = getNeighbourhood_v2(blockR, linea, T(col)+90, mb_size, 'LR');
                right(:,i+mb_size+1) = getDescriptors(linea, blockR, mb_size, T(col)+90, set2, set3, H(P(k,1),P(k,2)), xEnter, yEnter, xExit, yExit);
            end
        end

        [vL iL] = max(left(2,:).*varL);
        [vR iR] = max(right(2,:).*right(3,:));
        [vT iT] = max(top(2,:).*varT);
        [vB iB] = max(bottom(2,:).*varB);
        puntos = zeros(4,4);
        puntos(:,1) = leftPts(:,iL);
        puntos(:,2) = rightPts(:,iR);
        puntos(:,3) = topPts(:,iT);
        puntos(:,4) = bottomPts(:,iB);


        angVec = [left(1,iL) right(1,iR) top(1,iT) bottom(1,iB)];
        strVec = [vL vR vT vB];
        [w1 i] = max(strVec);
        ang1 = angVec(i);
        punto1 = puntos(:,i);
        strVec(i) = -1;
        [w2 i] = max(strVec);
        ang2 = angVec(i);
        punto2 = puntos(:,i);
        strVec(i) = -1;
        [w3 i] = max(strVec);
        if w3 > 0.5*w1
            ang3 = angVec(i);
        else
            ang3 = -1;
        end
        punto3 = puntos(:,i);
      
        N = 5;
        [ang wei puntos] = selectMax_WLOG(left, right, top, bottom);
        weights = getWeights_WLOG(ang, puntos, mb_size, N);

        if sum(ang(1:N)) == -N
            [block yaux] = biLinInterpol_v5(mask, [r c], 0, mb_size);
        else
            expo = 2;
            bloque = zeros(mb_size,mb_size,N);
            block = zeros(mb_size);
            for q = 1:N
                if ang(q) >= 0
                    [bloque(:,:,q) yaux] = biLinInterpol_v5(mask, [r c], ang(q), mb_size);
                end
            end      
            for l = 1:mb_size
                for m = 1:mb_size
                    p1 = wei/sum(wei);
                    p1 = p1(1:N);
                    p2 = zeros(1,N);
                    for q = 1:N
                        if ang(q) < 0
                            continue
                        else
                            p2(q) = (1-weights(l,m,q)^expo);
                        end
                    end

                    for q = 1:N
                        block(l,m) = block(l,m) + p1(q)*p2(q)/sum(p1.*p2) * bloque(l,m,q);
                    end

                end
            end
        end

        maskN5(r:r+mb_size-1,c:c+mb_size-1) = block;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = 4;
        [ang wei puntos] = selectMax_WLOG(left, right, top, bottom);
        weights = getWeights_WLOG(ang, puntos, mb_size, N);

        if sum(ang(1:N)) == -N
            [block yaux] = biLinInterpol_v5(mask, [r c], 0, mb_size);
        else
            expo = 2;
            bloque = zeros(mb_size,mb_size,N);
            block = zeros(mb_size);
            for q = 1:N
                if ang(q) >= 0
                    [bloque(:,:,q) yaux] = biLinInterpol_v5(mask, [r c], ang(q), mb_size);
                end
            end      
            for l = 1:mb_size
                for m = 1:mb_size
                    p1 = wei/sum(wei);
                    p1 = p1(1:N);
                    p2 = zeros(1,N);
                    for q = 1:N
                        if ang(q) < 0
                            continue
                        else
                            p2(q) = (1-weights(l,m,q)^expo);
                        end
                    end

                    for q = 1:N
                        block(l,m) = block(l,m) + p1(q)*p2(q)/sum(p1.*p2) * bloque(l,m,q);
                    end

                end
            end
        end
        maskN4(r:r+mb_size-1,c:c+mb_size-1) = block;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = 3;
        [ang wei puntos] = selectMax_WLOG(left, right, top, bottom);
        weights = getWeights_WLOG(ang, puntos, mb_size, N);

        if sum(ang(1:N)) == -N
            [block yaux] = biLinInterpol_v5(mask, [r c], 0, mb_size);
        else
            expo = 2;
            bloque = zeros(mb_size,mb_size,N);
            block = zeros(mb_size);
            for q = 1:N
                if ang(q) >= 0
                    [bloque(:,:,q) yaux] = biLinInterpol_v5(mask, [r c], ang(q), mb_size);
                end
            end      
            for l = 1:mb_size
                for m = 1:mb_size
                    p1 = wei/sum(wei);
                    p1 = p1(1:N);
                    p2 = zeros(1,N);
                    for q = 1:N
                        if ang(q) < 0
                            continue
                        else
                            p2(q) = (1-weights(l,m,q)^expo);
                        end
                    end

                    for q = 1:N
                        block(l,m) = block(l,m) + p1(q)*p2(q)/sum(p1.*p2) * bloque(l,m,q);
                    end

                end
            end
        end
        maskN3(r:r+mb_size-1,c:c+mb_size-1) = block;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N = 2;
        [ang wei puntos] = selectMax_WLOG(left, right, top, bottom);
        weights = getWeights_WLOG(ang, puntos, mb_size, N);

        if sum(ang(1:N)) == -N
            [block yaux] = biLinInterpol_v5(mask, [r c], 0, mb_size);
        else
            expo = 2;
            bloque = zeros(mb_size,mb_size,N);
            block = zeros(mb_size);
            for q = 1:N
                if ang(q) >= 0
                    [bloque(:,:,q) yaux] = biLinInterpol_v5(mask, [r c], ang(q), mb_size);
                end
            end      
            for l = 1:mb_size
                for m = 1:mb_size
                    p1 = wei/sum(wei);
                    p1 = p1(1:N);
                    p2 = zeros(1,N);
                    for q = 1:N
                        if ang(q) < 0
                            continue
                        else
                            p2(q) = (1-weights(l,m,q)^expo);
                        end
                    end

                    for q = 1:N
                        block(l,m) = block(l,m) + p1(q)*p2(q)/sum(p1.*p2) * bloque(l,m,q);
                    end

                end
            end
        end
        maskN2(r:r+mb_size-1,c:c+mb_size-1) = block;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end    

    waitbar(w/length(centers(1,:)),h,[num2str(round(100*w/length(centers(1,:)))) '% is done'])
    
end
close(h)
            
%% Showing results

% MS-SSIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = [0.01 0.03];
winsize = 8;
sigma = 0.5;
window = fspecial('gaussian', winsize, sigma);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'wtd_sum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
imshow(received_frame,[0 255])
subplot(1,2,2)
imshow(maskN5,[0 255])
title(['PSNR = ' num2str(psnr(img,maskN5)) 'dB   |   MS-SSIM = ' num2str(ssim_mscale_new(img, maskN5, K, window, level, weight, method)*100)])
toc
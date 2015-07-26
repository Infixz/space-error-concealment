%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2011 by                                          %
%   Jan Koloda                                                     %
%   Universidad de Granada, Granada, Spain                         %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   Li, X., Orchard, M.T., "Novel Sequential Error-Concealment     %
%   Techniques Using Orientation Adaptive Interpolation",          %
%   IEEE Transactions on Circuits and Systems for Video Technology %
%   October 2002, Vol. 12, No. 10                                  %
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

tic%% Loading the image and parsing it to grey level, if needed

img = imread('D:\MATlab\toolbox\images\imdemos\1.png');
[r c d] = size(img);
if d > 1
    img = double(rgb2gray(img));
else
    img = double(img);
end

%% Generating losses

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

%Distances (Dk) used in Eq.(7) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = zeros(mb_size*mb_size,4);
counter = 1;
for i = 1:mb_size
    for j = 1:mb_size
        D(counter,1) = i;               %D0
        D(counter,2) = j;               %D1
        D(counter,3) = mb_size-i+1;     %D2
        D(counter,4) = mb_size-j+1;     %D3
        counter = counter + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters described in Table 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [0 0 2 2 1 1 3 3];
b = [1 3 1 3 0 2 0 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Computing the weights omega according to Eq.(7) %%%%%%%%%%%%%%%%%%%%%%%%%%
omega = zeros(mb_size,mb_size,8);
for k = 1:8
    counter = 1;
    for i = 1:mb_size
        for j = 1:mb_size
            omega(i,j,k) = 0.25 - (3*D(counter,a(k)+1) + D(counter,b(k)+1))/(8*sum(D(counter,:)));
            counter = counter + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameter T, see Fig. 2
T = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0);
set(h,'Name','Processing')

%Concealment is carried out sequentially (macroblock by macroblock)
for k = 1:length(centers(1,:))        
    
    r = centers(1,k);
    c = centers(2,k);
    scans = zeros(mb_size,mb_size,8);
    
    %Eight different scaning orders (see Fig. 5):
    %D1
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for i = 1:mb_size
        for j = 1:mb_size            
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,1) = x;
        end
    end
    
    %D2
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for i = 1:mb_size
        for j = mb_size:-1:1
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,2) = x;
        end
    end
    
    %D3
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for i = mb_size:-1:1
        for j = 1:mb_size
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,3) = x;
        end
    end
    
    %D4
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for i = mb_size:-1:1
        for j = mb_size:-1:1
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,4) = x;
        end
    end
    
    %D5
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for j = 1:mb_size
        for i = 1:mb_size
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,5) = x;
        end
    end
    
    %D6
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for j = 1:mb_size
        for i = mb_size:-1:1
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,6) = x;
        end
    end
    
    %D7
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for j = mb_size:-1:1
        for i = 1:mb_size
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,7) = x;
        end
    end
    
    %D8
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    for j = mb_size:-1:1
        for i = mb_size:-1:1
            x = predictor(support_area,mb_size+i,mb_size+j,T);
            support_area(mb_size+i,mb_size+j) = x;
            scans(i,j,8) = x;
        end
    end
    
    %Combination of the eight scans
    recovery = zeros(mb_size,mb_size);
    for iter = 1:8
        recovery = recovery + scans(:,:,iter).*omega(:,:,iter);
    end
    
    mask(r:r+mb_size-1,c:c+mb_size-1) = recovery;
    waitbar(k/length(centers(1,:)),h,[num2str(round(100*k/length(centers(1,:)))) '% is done'])        
    
end
close(h)

% MS-SSIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = [0.01 0.03];
winsize = 8;
sigma = 0.25;
window = fspecial('gaussian', winsize, sigma);
level = 5;
weight = [0.0448 0.2856 0.3001 0.2363 0.1333];
method = 'wtd_sum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
imshow(received_frame,[0 255])
subplot(1,2,2)
imshow(mask,[0 255])
title(['PSNR = ' num2str(psnr(img,mask)) 'dB   |   MS-SSIM = ' num2str(ssim_mscale_new(img, mask, K, window, level, weight, method)*100)])




toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2011 by                                          %
%   Jan Koloda                                                     %
%   Universidad de Granada, Granada, Spain                         %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   Sun, H., Kwok, W., "Concealment of damaged block transform     %
%   coded images using projections onto convex sets",              %
%   IEEE Transactions on Image Processing, April 1995, pp.470-477  %
%                                                                  %
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
mode = 'specify';

%Cropping the image so it is made of an integer number of macroblocks %%%%%
img = img(1:floor(r/mb_size)*mb_size,1:floor(c/mb_size)*mb_size);
[r c] = size(img);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mask centers] = simuLoss(img, mb_size, nSlices, slice_to_be_lost, mode);
received_frame = mask;

%% Concealment

%NOTE: neighbourhood_N includes the missing macroblock and its 8 closest
%neighbouring macroblocks

%Parameters initialization (these values are suggested by the authors in
%Section IV (Simulation Results)
T = 5000;
Rth = 3;
Bth = 3;
num_iter = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Constructing the low-pass filter P2smooth, as described in Eq.(17) %%%%%%%
origin = 3*mb_size/2+1;
P2smooth = ones(3*mb_size);
for i = 1:3*mb_size
    for j = 1:3*mb_size
        if sqrt((i-origin)^2 + (j-origin)^2) > Rth
            P2smooth(i,j) = 0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0);
set(h,'Name','Processing')

%Concealment is carried out sequentially (macroblock by macroblock)
for counter = 1:length(centers(1,:))
    r = centers(1,counter);
    c = centers(2,counter);
    neighbourhood_N = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);

    theta = getAngle(neighbourhood_N, T, mb_size);
    
    %Initialization of "f" (f_0, see Section IV)
    f = neighbourhood_N;
    f(neighbourhood_N<0) = mean(mean(neighbourhood_N(neighbourhood_N>=0)));   
    
    %If the area around the missing macroblock is monotone (Eq.(10)), we
    %apply Eq.(20) using...
    if theta < 0
        for i = 1:num_iter         
            
            % ...P2smooth (Eq.(17)) and...
            fff = fftshift(fft2(f));
            fff = fff.*P2smooth;
            fff = ifft2(fftshift(fff));
            
            % ...P1 (Eq.(15))
            f(mb_size+1:2*mb_size,mb_size+1:2*mb_size) = fff(mb_size+1:2*mb_size,mb_size+1:2*mb_size);
            f(f < 0) = 0;
            f(f > 255) = 255;
        end

    %If the area around the missing macroblock is edgy (Eq.(10)), we
    %apply Eq.(20) using...
    else
        for i = 1:num_iter  
            
            %...P2edges (Eq.(19)) and...
            P2edges = ones(3*mb_size);
            for p = 1:3*mb_size
                for q = 1:3*mb_size
                    if abs(-(p-origin) - (q-origin)*tan(deg2rad(theta + 90))) > Bth
                        P2edges(p,q) = 0;
                    end
                end
            end
            fff = fftshift(fft2(f));
            fff = fff.*P2edges;
            fff = ifft2(fftshift(fff));
            
            %...P1 (Eq.(15))
            f(mb_size+1:2*mb_size,mb_size+1:2*mb_size) = real(fff(mb_size+1:2*mb_size,mb_size+1:2*mb_size));
            f(f < 0) = 0;
            f(f > 255) = 255;
        end
    end
    waitbar(counter/length(centers(1,:)),h,[num2str(round(100*counter/length(centers(1,:)))) '% is done'])
    mask(r:r+mb_size-1,c:c+mb_size-1) = f(mb_size:mb_size+mb_size-1,mb_size:mb_size+mb_size-1);    
end
close(h)

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
imshow(mask,[0 255])
title(['PSNR = ' num2str(psnr(img,mask)) 'dB   |   MS-SSIM = ' num2str(ssim_mscale_new(img, mask, K, window, level, weight, method)*100)])

toc






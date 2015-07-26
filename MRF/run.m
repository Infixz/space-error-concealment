%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2011 by                                          %
%   Jan Koloda                                                     %
%   Universidad de Granada, Granada, Spain                         %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   Shirani, S., Kossentini, F., Ward, R., "An adaptive Markov     %
%   random field based error concealment method for video          %
%   communication in an error prone environment", ICASSP 1999      %
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

tic%% Loading the image and parsing it to grey level, if necessary

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

% parameters 
alpha = 1;          %Eq. (13)
max_iter = 10;      %Number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0);
set(h,'Name','Processing')

%Concealment is carried out sequentially (macroblock by macroblock)
for w = 1:length(centers(1,:))    
    r = centers(1,w);
    c = centers(2,w);
    
    %support_area includes the missing macroblock and its 8 closest
    %neighbouring macroblocks
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);
    
    %Getting "cm" parameter, Eq. (13)
    cm = getCounters(support_area, mb_size);
    
    %Iterative reconstruction
    for iter = 1:max_iter
        aux = support_area;
        %For every lost pixel...
        for i = 1:mb_size
            for j = 1:mb_size
                %...we get the two cliques...
                cliques = getClique(support_area, mb_size+i, mb_size+j);
                numerator = 0;
                denominator = 0;
                %...and conceal (Eq. (14))
                for k = 1:8
                    if cliques(1,k) >= 0
                        numerator = numerator + cliques(1,k)*cm(k)*alpha;
                        denominator = denominator + cm(k)*alpha;
                    end
                    if cliques(2,k) >= 0
                        numerator = numerator + cliques(2,k)*cm(k)*alpha;
                        denominator = denominator + cm(k)*alpha;
                    end
                end
                support_area(mb_size+i,mb_size+j) = numerator/denominator;
            end
        end
    end
    mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1) = support_area;
    waitbar(w/length(centers(1,:)),h,[num2str(round(100*w/length(centers(1,:)))) '% is done'])
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
            
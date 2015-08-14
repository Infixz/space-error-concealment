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

angBins = [0 22.5 45 67.5 90 112.5 135 157.5];
h = waitbar(0);
set(h,'Name','Processing')
for w = 1:length(centers(1,:))   
    
    r = centers(1,w);
    c = centers(2,w);      
    
    %support_area includes the missing macroblock and its 8 closest
    %neighbouring macroblocks
    support_area = mask(r-mb_size:r+mb_size+mb_size-1,c-mb_size:c+mb_size+mb_size-1);

    % Content estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d = getCounters(support_area, mb_size);
    [content angSelector] = decideMode(d);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    block_conc = zeros(mb_size);
    %Smooth content, see Eq.(1)
    if content == 1
        block_conc = bilinearInterpolation(support_area(mb_size:2*mb_size+1,mb_size:2*mb_size+1), 0.5);
    
    %Content with edges
    elseif content == 2
        for i = 1:8
            if angSelector(i) > 0
                %See Eq.(2)
                block_aux = directionalInterpolation(support_area, 180-angBins(i), mb_size);
                block_conc = block_conc + block_aux*d(i);
            end
        end
        %See Eq.(3)
        block_conc = block_conc/sum(d.*(angSelector>0));
        
    %Content with texture, see Eq.(4)    
    elseif content == 3
        block_conc = bestMatch([r c], received_frame, mb_size);
    end
    
    mask(r:r+mb_size-1,c:c+mb_size-1) = block_conc;    
    waitbar(w/length(centers(1,:)),h,[num2str(round(100*w/length(centers(1,:)))) '% is done'])
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
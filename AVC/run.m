%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                  %
%   Copyright (c) 2011 by                                          %
%   Jan Koloda                                                     %
%   Universidad de Granada, Granada, Spain                         %
%   - all rights reserved -                                        %
%                                                                  %
%   This is an implementation of the algorithm described in:       %
%   Salama, P., Shroff, N.B., Coyle, E.J., Delp, E.J., "Error      %
%   concealment techniques for encoded video streams", ICIP 1995   %
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

%% Loading the image and parsing it to grey level, if needed

[filename,pathname]=uigetfile({'*.jpg';'*.bmp';'*.tif';'*.*'},'ͼ���ļ�');
file=[pathname,filename];
X=imread(file);
img = imread('D:\MATlab\toolbox\images\imdemos\1.png');
[r c d] = size(img);
if d > 1
    img = double(rgb2gray(img));
else
    img = double(img);
end

%% Generating losses
tic    %��ʼ��ʱ

slice_to_be_lost = 1;
nSlices = 2;
mb_size = 16;
mode = 'default';

%Cropping the image so it is made of an integer number of macroblocks %%%%%
img = img(1:floor(r/mb_size)*mb_size,1:floor(c/mb_size)*mb_size);
[rows, cols] = size(img);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mask centers] = simuLoss(img, mb_size, nSlices, slice_to_be_lost, mode);
received_frame = mask;

%% Concealment

% lambda parameter, see Eq.(1), NOTE: lambda = 0.5 corresponds to the
% non-normative H.264 concealment technique
lambda = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = waitbar(0);
set(h,'Name','Processing')

%Concealment is carried out sequentially (macroblock by macroblock)
for w = 1:length(centers(1,:))
    waitbar(w/length(centers(1,:)),h,[num2str(round(100*w/length(centers(1,:)))) '% is done'])
    r = centers(1,w);
    c = centers(2,w);
    
    %support_area contains the available "ring" of pixels around the
    %missing macroblocks (excluding the four corners)
    support_area = -ones(mb_size+2);
    support_area(2:mb_size+1,2:mb_size+1) = mask(r:r+mb_size-1,c:c+mb_size-1);
    if r > 1
        support_area(1,2:end-1) = mask(r-1,c:c+mb_size-1);
    end
    if r < rows - mb_size + 1
        support_area(end,2:end-1) = mask(r+mb_size,c:c+mb_size-1);
    end
    if c > 1
        support_area(2:end-1,1) = mask(r:r+mb_size-1,c-1);
    end
    if c < cols - mb_size + 1
        support_area(2:end-1,end) = mask(r:r+mb_size-1,c+mb_size);
    end
    
    conc_block = bilinearInterpolation(support_area, lambda);    
    mask(r:r+mb_size-1,c:c+mb_size-1) = conc_block;    
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
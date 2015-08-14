%This function simulates damages in the received frame corresponding to the
%loss of one NALU (i.e. slice) given a slicing mode and the slice to be
%lost

function [y centers] = simuLoss(img, mb_size, nSlices, nLoss, mode)

% INPUT:
%   img     - original frame
%   nSlice  - number of slices the image is made of
%   nLoss   - indicates the slice that is lost
%   mb_size - macroblock dimensions

% OUTPUT:
%   y       - damaged frame, missing pixels are set to -1
%   centers - 2xN vector, where N is the number of missing macroblocks
%           - [row_1 row_2 ... row_N; 
%              col_1 col_2 ... col_N]
%           - row_i (col_i) indicates the row (column) of the top-left
%            corner of the i^th missing macroblock
%

if nargin < 3
    mode = 'default';
end

y = img;
[rows cols] = size(img);

%Dispersed slicing, as described in Recommendation H.264
if strcmp(mode, 'dispersed')
        
    %Extracting variables needed for slicing (see Recommendation H.264)
    PicWidthInMbs = cols/mb_size;
    PicHeightInMbs = rows/mb_size;
    PicSizeInMapUnits = PicWidthInMbs * PicHeightInMbs;
    num_slice_groups_minus1 = nSlices - 1;
    mapUnitToSliceGroupMap = zeros(1,PicSizeInMapUnits);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %see 8.2.2.2 "Specification for dispersed slice group map type" in
    %Recommendation H.264 (p. 116)
    for i = 0:PicSizeInMapUnits-1
        mapUnitToSliceGroupMap(i+1)  = mod((mod(i,PicWidthInMbs) + floor((floor(i/PicWidthInMbs)*(num_slice_groups_minus1+1))/2)),(num_slice_groups_minus1+1));       
    end
    
    %Setting the missing pixels to -1 and computing "centers"
    counter = 1;
    counter_centers = 1;
    centers = zeros(2,sum(mapUnitToSliceGroupMap(counter)+1 == nLoss));
    for i = 1:mb_size:rows-mb_size+1
        for j = 1:mb_size:cols-mb_size+1                        
                        
            if mapUnitToSliceGroupMap(counter)+1 == nLoss                            
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);
                centers(:,counter_centers) = [i;j];                
                counter_centers = counter_centers + 1;
            end
            counter = counter + 1;
            
        end
    end
    
%Slicing for testing purposes, the "famous" 25% loss where every missing
%macroblock has all its 8 closest neighbouring macroblocks available
% (REFERENCE) XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
elseif strcmp(mode, 'default')
    auxX = floor(cols/(2*mb_size));
    auxY = floor(rows/(2*2*mb_size));
    centers = zeros(2,auxX*auxY);
    counter = 1;

    for i = 1+mb_size:2*mb_size:rows
        for j = 1:2*mb_size:cols
            if i > 1 && i < rows - mb_size && j > 1 && j < cols - mb_size
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);                
                centers(:,counter) = [i;j];
                counter = counter + 1;
            end
        end
    end

    
elseif strcmp(mode, 'specify')
    % 本模式基于FMO TYPE2 丢失模式，构造成图像roi区域丢失
    PicWidthInMbs = cols/mb_size;    % 以块大小计量图片：宽度
    PicHeightInMbs = rows/mb_size;   % 高度
    PicSizeInMapUnits = PicWidthInMbs * PicHeightInMbs;    % 用块大小计算图像面积
    num_slice_groups_minus1 = nSlices - 1;    % nSlices:图像由多少slice组成
    mapUnitToSliceGroupMap = zeros(1,PicSizeInMapUnits);    %散布图做编号
    % configure
    mapUnitToSliceGroupMap(461:467) = 1;
    mapUnitToSliceGroupMap(461+32:467+32) = 1;
    mapUnitToSliceGroupMap(461+32*2:467+32*2) = 1;
    mapUnitToSliceGroupMap(461+32*3:467+32*3) = 1;
    
    
    % C实现为：
    % for( i = 0; i < p_Img->PicSizeInMapUnits; i++ ){
    % p_Img->MapUnitToSliceGroupMap [i] = (byte) (((i%p_Img->PicWidthInMbs)+(((i/p_Img->PicWidthInMbs)*(pps->num_slice_groups_minus1+1))>>1))
    %  %(pps->num_slice_groups_minus1+1));}
    
    %for i = 0:PicSizeInMapUnits-1
    %    mapUnitToSliceGroupMap(i+1)  = mod((mod(i,PicWidthInMbs) + floor((floor(i/PicWidthInMbs)*(num_slice_groups_minus1+1))/2)),(num_slice_groups_minus1+1));       
    %end
    
    %Setting the missing pixels to -1 and computing param:centers:
    counter = 1;
    counter_centers = 1;
    
    centers = zeros(2,sum(mapUnitToSliceGroupMap(counter)+1 == nLoss));    
    %centers?什么作用，为什么它是二维矩阵
    
    for i = 1:mb_size:rows-mb_size+1
        for j = 1:mb_size:cols-mb_size+1
                        
            if mapUnitToSliceGroupMap(counter) == nLoss                            
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);% 选定块置 -1
                centers(:,counter_centers) = [i;j];                
                counter_centers = counter_centers + 1;
            end
            counter = counter + 1;
            
        end
    end
        
    
else
    display('Error: Invalid slicing mode')
end


end

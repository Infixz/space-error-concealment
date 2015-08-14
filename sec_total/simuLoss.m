%This function simulates damages in the received frame corresponding to the
%loss of one NALU (i.e. slice) given a slicing mode and the slice to be
%lost

function [y, centers] = simuLoss(img, mb_size, nSlices, nLoss, mode, loss_rate)

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

%% ��ʼ�����ж�ģʽ�����ñ���
if nargin < 3   % NARGIN��Number of function input arguments
    mode = 'default';
elseif nargin < 6
    loss_rate = 0.10;
end

y = img;
[rows, cols] = size(img);%�õ�ͼ��ʵ�ʵĳ���

%% Dispersed slicing, as described in Recommendation H.264
% �����ֲ�
if strcmp(mode, 'dispersed')    %�Ƚ��ַ���ѡ��ģʽ
    %Extracting variables needed for slicing (see Recommendation H.264)
    PicWidthInMbs = cols/mb_size;    %�Կ��С����ͼƬ������
    PicHeightInMbs = rows/mb_size;   % �߶�
    PicSizeInMapUnits = PicWidthInMbs * PicHeightInMbs;    %�ÿ��С����ͼ�����
    num_slice_groups_minus1 = nSlices - 1;    %nSlices:ͼ���ɶ���slice���
    mapUnitToSliceGroupMap = zeros(1,PicSizeInMapUnits);    %ɢ��ͼ�����
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %see 8.2.2.2 "Specification for dispersed slice group map type" in
    %Recommendation H.264 (p. 116)
    
    % Cʵ��Ϊ��
    % for( i = 0; i < p_Img->PicSizeInMapUnits; i++ ){
    % p_Img->MapUnitToSliceGroupMap [i] = (byte) (((i%p_Img->PicWidthInMbs)+(((i/p_Img->PicWidthInMbs)*(pps->num_slice_groups_minus1+1))>>1))
    %  %(pps->num_slice_groups_minus1+1));}
    for i = 0:PicSizeInMapUnits-1
        % ͨ����ż���趨slice group = 0 or 1
        mapUnitToSliceGroupMap(i+1)  = mod((mod(i,PicWidthInMbs) + floor((floor(i/PicWidthInMbs)*(num_slice_groups_minus1+1))/2)),(num_slice_groups_minus1+1));       
    end
    
    %Setting the missing pixels to -1 and computing param:centers:
    counter = 1;
    counter_centers = 1;
    
    centers = zeros(2,sum(mapUnitToSliceGroupMap(counter)+1 == nLoss));    
    %centers?ʲô���ã�Ϊʲô���Ƕ�ά����
    
    for i = 1:mb_size:rows-mb_size+1
        for j = 1:mb_size:cols-mb_size+1
                        
            if mapUnitToSliceGroupMap(counter)+1 == nLoss                            
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);% ѡ������ -1
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
    
    
elseif strcmp(mode, 'random')
    % ���ģʽ������������
    losses = rand(1,rows*cols/mb_size^2);
    [~, indices] = sort(losses, 'descend');
    indices = indices(1:round(length(indices)*loss_rate));    %Ĭ�϶�ʧ��0.1
    centers = zeros(2,length(indices));
    indices = sort(indices,'ascend');
    iterator = 1;
    counter = 1;
    
    
    for i = 1:mb_size:rows - mb_size
        for j = 1:mb_size:cols - mb_size
            if counter <= length(indices) && iterator == indices(counter)
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);                
                centers(:,counter) = [i;j];
                counter = counter + 1;
            end
            iterator = iterator + 1;
        end
    end

elseif strcmp(mode, 'specify')
    % ��ģʽ����FMO TYPE2 ��ʧģʽ�������ͼ��roi����ʧ
    PicWidthInMbs = cols/mb_size;    % �Կ��С����ͼƬ������
    PicHeightInMbs = rows/mb_size;   % �߶�
    PicSizeInMapUnits = PicWidthInMbs * PicHeightInMbs;    % �ÿ��С����ͼ�����
    num_slice_groups_minus1 = nSlices - 1;    % nSlices:ͼ���ɶ���slice���
    mapUnitToSliceGroupMap = zeros(1,PicSizeInMapUnits);    %ɢ��ͼ�����
    % configure
    mapUnitToSliceGroupMap(461:467) = 1;
    mapUnitToSliceGroupMap(461+32:467+32) = 1;
    mapUnitToSliceGroupMap(461+32*2:467+32*2) = 1;
    mapUnitToSliceGroupMap(461+32*3:467+32*3) = 1;
    
    
    % Cʵ��Ϊ��
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
    %centers?ʲô���ã�Ϊʲô���Ƕ�ά����
    
    for i = 1:mb_size:rows-mb_size+1
        for j = 1:mb_size:cols-mb_size+1
                        
            if mapUnitToSliceGroupMap(counter) == nLoss                            
                y(i:i+mb_size-1,j:j+mb_size-1) = -ones(mb_size);% ѡ������ -1
                centers(:,counter_centers) = [i;j];                
                counter_centers = counter_centers + 1;
            end
            counter = counter + 1;
            
        end
    end
    
else    %end of if
    display('Error: Invalid slicing mode')
end

end
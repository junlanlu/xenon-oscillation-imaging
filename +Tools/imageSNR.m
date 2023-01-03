function [imgSNR, imgSNR_Rayleigh, imgSignal, imgNoise] = imageSNR(image, mask, noiseSize, min_numNoise)

% image
% mask: binary
% noiseSize: size of one side of noise ROI
% min_numNoise: min voxels in noise cube to use ROI in SNR calculation 
 
imSize = size(image);   

% set up dilate/erode paramters
SE = 1; % structuring element size
[x,y,z]=meshgrid(-SE:SE,-SE:SE, -SE:SE);
nhood=x.^2+y.^2+z.^2 <=SE^2;                
se1=strel('arbitrary',nhood);
% dilate mask and create masked image for noise calculation
dilated_mask =imdilate(mask,se1);
dilated_masked_image = image;
dilated_masked_image(dilated_mask ~= 0) = NaN;


% erode mask and create masked image for signal calculation
eroded_mask = imerode(mask,se1);

% calculate the index values of the first noise voxel
k = 0:noiseSize-1;
nosie_x = k.*imSize(2)+ones(noiseSize,1);
noise_xy = nosie_x + (k.*ones(noiseSize,1))';
noise_xyz = repmat(noise_xy(:),1,noiseSize) + k*imSize(1)*imSize(2);
noiseVoxel = noise_xyz(:);

% calculate the starting index for all noise voxelsnumNoiseROIs = floor(imSize/noiseSize);
numNoiseROIs = imSize/noiseSize;
startingVoxel_x = permute((0:numNoiseROIs(1)-1).*noiseSize+ones(numNoiseROIs(1),1),[2 1]);
startingVoxel_xy = startingVoxel_x+(0:numNoiseROIs(2)-1).*noiseSize*imSize(1);
startingVoxel_xyz = repmat(startingVoxel_xy(:),1,numNoiseROIs(3)) + (0:numNoiseROIs(3)-1)*imSize(1)*imSize(2)*noiseSize;
startingVoxel = startingVoxel_xyz(:);

% move initial voxel to all starting position
noiseIndex = permute(noiseVoxel' + startingVoxel-1',[2 1]);

% remove noiseIndex with too few voxels
noiseValues = dilated_masked_image(noiseIndex);
count_numNoise = sum(~isnan(noiseValues));
noiseValues(:,count_numNoise < min_numNoise) = [];

imgNoise = median(nanstd(noiseValues));
imgSignal = mean(image(eroded_mask ~= 0));
imgSNR = imgSignal/imgNoise;
imgSNR_Rayleigh = imgSNR * 0.66;

%% This was used for testing the calculated noise index positions

% clear
% 
% noiseSize = 8;               % size of one size of noise cube
% min_numNoise = noiseSize^3;   % min voxels in noise cube 
% 
% % dilate mask
% 
% image = ones(128,128,128);
% image(54:74,54:74,54:74) = 0;
% mask = zeros(128,128,128);
% mask(54:74,54:74,54:74) = 1;
% 
% image(mask ~= 0) = NaN;
% 
    % imSize = size(image);
    % noiseSize = 128/8;
    % 
    % 
    % imageIndex = reshape(1:length(image(:)),[128, 128, 128]);
    % for k = 1:ndims(imageIndex)
    %     c = floor(imSize(k)/noiseSize);
    %     d = rem(imSize(k), noiseSize);
    %     partition(k,:) = ones(1, noiseSize)*c;
    %     partition(k,1:d) = partition(1:d)+1;
    % end 
    % 
    % noiseROIs = mat2cell(imageIndex, partition(1,:), partition(2,:),partition(3,:));
    % % test = cellfun(@std2,noiseROIs);
    % 
    % noiseIndex_matrix = cellfun(@(x) reshape(x,[],1),noiseROIs,'UniformOutput',false);
    % noiseIndex_vector = reshape(noiseIndex_matrix,1,[]);
    % noiseIndex2 = cell2mat(noiseIndex_vector);


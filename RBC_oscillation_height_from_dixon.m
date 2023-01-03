%% get paths and load data
if ~exist('pathbundle','var') || ~isfield(pathbundle,'theDixonFile') 
    disp('Locate the Dixon file');
    pathbundle.theDixonFile = filepath();
    pathbundle.filter = 0;
end 

%% decompose the path and names 

[pathbundle.filepath,pathbundle.filename] = fileparts(pathbundle.theDixonFile);
[upperPath, deepestFolder, ~] = fileparts(pathbundle.filepath);
pathbundle.subname = deepestFolder;
k = strfind(pathbundle.filename,'_FID');
pathbundle.filename = pathbundle.filename(k:k+8); clear k;

%% read in data 

twix_obj = mapVBVD(pathbundle.theDixonFile);
twix_obj.image.flagIgnoreSeg = true; % This line ignore the extra dimension that will explode your memory.
twix_obj.image.flagRemoveOS = false; % This line removes oversampling
data = squeeze(double(twix_obj.image())); %read out image data 
nSpec = 1;

data_dis = data(:,4:2:end-1); % dissolved
k0 = data_dis(1,200:end);
k0_smooth = smooth(k0,5);

phaseShift = -pi:pi/20:pi;

for k = 1:length(phaseShift);
     
    k0_shifted = k0_smooth.*exp(1i*phaseShift(k));
    RBC = real(k0_shifted);
    barrier = imag(k0_shifted);
    
    bar_rmse(k) = std(barrier);
    
    k = k + 1;
end 

[xx yy] = min(bar_rmse);
k0_shifted = k0_smooth.*exp(1i*phaseShift(yy));
RBC = real(k0_shifted);
barrier = imag(k0_shifted);

figure(); plot(RBC), title('RBC')
figure(); plot(barrier), title('Barrier')



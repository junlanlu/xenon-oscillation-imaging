function dyn = fitDynamicSpec(raw_path, dyn_save_name)

% raw_path: full filepath of the raw scanner data
% save_name: name of .mat file for saving. The output will be saved in
%            the same folder ans the raw data

peaks = 3;
nToAvg = 5; % 1 for high flip
skipSize = 1;

[raw_folder, ~, scanner] = fileparts(raw_path);
% Read in twix or P file and define associated variables
[raw_fids, dwell_time, npts, tr, xeGam] = readRawDyn(raw_path);

% For high flip, also comment out SIFT
% raw_fids = raw_fids(10:end,5:end);

% SIFT raw fids
raw_fids = raw_fids(:,1:end-21);
raw_fids = SIFT(raw_fids, dwell_time, tr);
% raw_fids(:,901:end-1) = []; % shorten long aquisitions

% Separate fids from gas frames
gas_fid = mean(raw_fids(:,end-1:end),2);

% Separate dissolved frames
fids = raw_fids(:,1:(end-2));

npts = size(raw_fids,1);            % Number of samples                 
t = dwell_time*(0:(npts-1))';       % FID time vector
nFrames = size(fids,2);             % Number of disolved frames
t_tr = tr*(1:nFrames);              % tr time vector

switch peaks
    case 3 
%         Three peak fit
        area_orig = [1 1 1];
        freq_orig = [0 -20.7 -218.4]*xeGam;
%         freq_orig = [0 -20 -218]*xeGam+7000; For high flip
        fwhmL_orig = [8.7 5.0 1.2]*xeGam;
        fwhmG_orig = [0 6.1 0]*xeGam;
        phase_orig = [0 0 0];
        titles = {'RBC' 'Barrier' 'Gas'};
    case 4
        % Four peak fit
        area_orig = [1 1 1 1];
        freq_orig = [-20 -285 -393 -3600]*2;
        fwhm_orig = [215 200 150 70];
        phase_orig = [0 0 0 0];
        titles = {'RBC' 'Barrier 1' 'Barrier 2' 'Gas'};
    case 5
        % Five peak fit
        area_orig = [1 1 1 1 1];
        freq_orig = [-20 -285 -393 -3465 -3842]*2;
        fwhm_orig = [215 200 150 70 30];
        phase_orig = [0 0 0 0 0 ];
        titles = {'RBC' 'Barrier 1' 'Barrier 2' 'Gas 1' 'Gas 2'};
end 

%% Average data for high SNR and fit for starting guesses
avgdata = mean(fids(:,150:250),2);
nmrFit = NMR_TimeFit_v(avgdata, t, area_orig,freq_orig,fwhmL_orig,fwhmG_orig,phase_orig,[],[]);
nmrFit = nmrFit.fitTimeDomainSignal();

% Create guesses from fitted data
area = nmrFit.area; freq = nmrFit.freq; fwhmL = nmrFit.fwhm;
fwhmG = nmrFit.fwhmG; phase = nmrFit.phase; 

% Define reference frequency for ppm conversion
ref_freq = nmrFit.freq(end);

% Fit spectra
nComp = length(area);
startingTimePoints = 1:skipSize:(nFrames-nToAvg);
nTimePoints = length(startingTimePoints);
           
area_dyn = zeros(nTimePoints,nComp);
freq_dyn = zeros(nTimePoints,nComp);
fwhmL_dyn = zeros(nTimePoints,nComp);
fwhmG_dyn = zeros(nTimePoints,nComp);
phase_dyn = zeros(nTimePoints,nComp);
snrs = zeros(nTimePoints,1);
residualSpectrum = zeros(npts,nTimePoints);
tzero = min(t(:));

t_dyn = repmat(t_tr(startingTimePoints)',[1, nComp]);

for k = 1:peaks
    disp([sprintf('%8.3f', nmrFit.area(k)/sum(nmrFit.area(2))), ' ' ...
        sprintf('%8.1f',(nmrFit.freq(k)-ref_freq)/xeGam),  '  ' ...
        sprintf('%8.1f',nmrFit.fwhm(k)/xeGam), '  ' ...
        sprintf('%8.1f',nmrFit.fwhmG(k)/xeGam), '  ' ...
        sprintf('%9.1f',nmrFit.phase(k))]);

end

avgdata = movmean(fids,nToAvg,2,'Endpoints','discard');

tic
parfor iTimePoint = 1:nTimePoints
    nmrFit = NMR_TimeFit_v(avgdata(:,iTimePoint), t, area,freq,fwhmL,fwhmG,phase,[],[]);
    nmrFit = nmrFit.fitTimeDomainSignal();

    area_dyn(iTimePoint,:) = nmrFit.area(:);
    freq_dyn(iTimePoint,:) = (nmrFit.freq(:)-ref_freq)/xeGam; %ppm units
    fwhmL_dyn(iTimePoint,:) = nmrFit.fwhm(:)/xeGam; %ppm units
    fwhmG_dyn(iTimePoint,:) = nmrFit.fwhmG(:)/xeGam; %ppm units
    phase_dyn(iTimePoint,:) = nmrFit.phase(:);
    
    zeroPaddedTime = tzero + dwell_time*((1:nmrFit.zeroPadSize)-1)';
    fittedSignal = nmrFit.calcTimeDomainSignal(zeroPaddedTime);
    
    fittedSpectrum = dwell_time*fftshift(fft(fittedSignal));
    residualSpectrum(:,iTimePoint) = (nmrFit.spectralDomainSignal - fittedSpectrum);
    snrs(iTimePoint) = snr(fittedSignal,avgdata(:,iTimePoint)-fittedSignal);   
end
toc

% Create Structure for Easy Saving of Data
fnames = {'area','freq','fwhmL','fwhmG','phase','t','resSpec','ref_freq','snrs'};

dyn.(fnames{1}) = area_dyn;
dyn.(fnames{2}) = freq_dyn;
dyn.(fnames{3}) = fwhmL_dyn;
dyn.(fnames{4}) = fwhmG_dyn;
dyn.(fnames{5}) = phase_dyn;
dyn.(fnames{6}) = t_dyn;
dyn.(fnames{7}) = residualSpectrum;
dyn.(fnames{8}) = ref_freq;
dyn.(fnames{9}) = snrs;

dyn.processed_date = datestr(now);
dyn.raw_file = raw_path;

save(fullfile(raw_folder,dyn_save_name),'dyn')

 

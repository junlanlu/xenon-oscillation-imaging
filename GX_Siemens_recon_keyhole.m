function GX_Siemens_recon_keyhole(params)
recon_size = 128 * [1, 1, 1];
overgrid_factor_high = 4;

%% read in data
% Load proton mask
if contains(params.mask_file, '.mat')
    load(params.mask_file, 'mask_reg')
elseif contains(params.mask_file, '.nii')
    mask_reg = double(niftiread(params.mask_file));
else
    ME = MException('myComponent:inputError', 'No valid mask file');
    throw(ME);
end
protonMask = permute((mask_reg), [1, 3, 2]);
% Get rbc2barrier from dynamic spectroscopy
if exist([fileparts(params.spect_file), '\BHs.mat'], 'file')
    load([fileparts(params.spect_file), '\BHs.mat'])
else
    BHs = [2, 8];
end
[nmrFit, ~, fids, spect_tr] = calculateStaticSpectroscopy(params.spect_file, BHs, 'V');

if isfield(params, 'rbc2bar')
    rbc2bar = params.rbc2bar;
elseif spect_tr == 0.020
    rbc2bar = nmrFit.area(1) / nmrFit.area(2) * 0.89; % scale down rbc2bar to account for different TR
elseif spect_tr == 0.015
    rbc2bar = nmrFit.area(1) / nmrFit.area(2);
elseif spect_tr == 0.073
    rbc2bar = nmrFit.area(1) / nmrFit.area(2) * 0.53;
else
    rbc2bar = 0.6002;
end

[subjFolder, Subject] = fileparts(params.mask_file);
if contains(Subject, '_gx')
    Subject = erase(Subject, '_gx');
end
if contains(Subject, '_highBW')
    Subject = erase(Subject, '_highBW');
    Subject = append(Subject, '_s1');
end
if contains(Subject, '_hrep')
    Subject = erase(Subject, '_hrep');
    Subject = append(Subject, '_s2');
end
if isfield(params, 'subject_name')
    Subject = params.subject_name;
end
ReportTitle = [Subject, ' Keyhole Summary Detrend'];

twix_obj = mapVBVD(params.dixon_file);
twix_obj.image.flagIgnoreSeg = true; % This line ignore the extra dimension that will explode your memory.
twix_obj.image.flagRemoveOS = false; % This line removes oversampling
protocol_name = twix_obj.hdr.Config.ProtocolName; % get the sequence name
data = squeeze(double(twix_obj.image.unsorted())); %read out image data
nSpec = 1;

%% Calculate Trajectory and remove bad rays
[data_dis, data_gas, traj_dis, traj_gas, data_type] = calculateTrajDixon(twix_obj, data, nSpec);
[npts, nFrames_raw] = size(data);
% scale trajectory to image output size
scale_size = npts ./ recon_size;
traj_gas_scale = traj_gas ./ scale_size;
traj_dis_scale = traj_dis ./ scale_size;

%  Enforce Nyquist (no trajectory > +/- 0.5)
%[traj_gas_scale, data_gas] = MRI.DataProcessing.enforceNyquistBounds(traj_gas_scale, data_gas);
%[traj_dis_scale, data_dis] = MRI.DataProcessing.enforceNyquistBounds(traj_dis_scale, data_dis);

% scale the proton mask if needed
%if npts == 64
%     protonMask = imresize3(double(protonMask),[npts npts npts],'nearest');
%     protonMask(protonMask <= 0.75) = 0; protonMask(protonMask > 0.75) = 1; protonMask = logical(protonMask);
% end

%% Prep variables for keyhole analysis
disFID = conj(reshape(data_dis, [npts, length(data_dis) / npts])) * 1e5;
gasFID = conj(reshape(data_gas, [npts, length(data_gas) / npts])) * 1e5;
traj = reshape(traj_dis_scale, [3, npts, length(data_dis) / npts]);

% Start by cutting out approach to steady state
N_Begin_Cutoff = params.N_Begin_Cutoff;

%For consistency, cut out points in Gas, Dissolved, and trajectories
disFID(:, 1:N_Begin_Cutoff) = [];
gasFID(:, 1:N_Begin_Cutoff) = [];
traj(:, :, 1:N_Begin_Cutoff) = [];
traj_dis_scale_keyhole = traj_dis_scale(npts*N_Begin_Cutoff+1:end, :);
traj_dis_keyhole = traj_dis(npts*N_Begin_Cutoff+1:end, :);

% Cut out views at the end in case of poor breath-hold
N_End_Cutoff = params.N_End_Cutoff;

%For consistency, cut out points in Gas, Dissolved, and trajectories
if N_End_Cutoff > 0
    disFID(:, size(disFID, 2)-N_End_Cutoff+1:end) = [];
    gasFID(:, size(gasFID, 2)-N_End_Cutoff+1:end) = [];
    traj(:, :, size(traj, 3)-N_End_Cutoff+1:end) = [];
    traj_dis_scale_keyhole = traj_dis_scale_keyhole(1:end-npts*N_End_Cutoff, :);
    traj_dis_keyhole = traj_dis_keyhole(1:end-npts*N_End_Cutoff, :);
end

%% Bin into High and Low keys based on RBC Oscillations

% get the TR
TR = 2 * twix_obj.hdr.Config.TR(1) * 1e-6; % in seconds

% bin the data
[data_dis_highkey, data_dis_lowkey, data_dis_tot, rbc2bar_all, Shift_Raw_FIDs_Fig, Smooth_Detrend_Fig, Key_Radius] = ...
    dissolved_phase_RBC_detrend_bin_peaks(traj, gasFID, disFID, TR, rbc2bar, params);
% remove data points without information
high_zeros = find(data_dis_highkey == 0);
traj_dis_highkey = traj_dis_keyhole;
traj_dis_highkey(high_zeros, :) = [];
data_dis_highkey(high_zeros) = [];

low_zeros = find(data_dis_lowkey == 0);
traj_dis_lowkey = traj_dis_keyhole;
traj_dis_lowkey(low_zeros, :) = [];
data_dis_lowkey(low_zeros) = [];

%% Recon
reconVol_dis_tot = reconDixonLowRes(recon_size(1), data_dis_tot(:), traj_dis_keyhole);
imslice(abs(reconVol_dis_tot));
title('Dissolved Recon Total');
reconVol_dis_highkey = reconDixonLowRes(recon_size(1), data_dis_highkey(:), traj_dis_highkey);
reconVol_dis_lowkey = reconDixonLowRes(recon_size(1), data_dis_lowkey(:), traj_dis_lowkey);
% use conjugate so images have proper orientation
reconVol_dis = reconDixonLowRes(recon_size(1), conj(data_dis), traj_dis);
reconVol_gas_highSNR = reconDixonLowRes(recon_size(1), conj(data_gas)*1e5, traj_gas);
reconVol_gas_highreso = reconDixonHighRes(recon_size(1), conj(data_gas)*1e5, traj_gas);

figure(3);
imslice(abs(reconVol_dis_highkey));
title('Dissolved Recon High Key');
figure(4);
imslice(abs(reconVol_dis_lowkey));
title('Dissolved Recon Low Key');
figure(5);
imslice(abs(reconVol_dis));
title('Dissolved Recon');
figure(6);
imslice(abs(reconVol_gas_highSNR));
title('Gas Recon High SNR');
figure(7);
imslice(abs(reconVol_gas_highreso));
title('Gas Recon High Resolution');

%%

reconVol_dis_tot = permute(reconVol_dis_tot, [1, 3, 2]);
reconVol_dis_highkey = permute(reconVol_dis_highkey, [1, 3, 2]);
reconVol_dis_lowkey = permute(reconVol_dis_lowkey, [1, 3, 2]);

reconVol_dis = permute(reconVol_dis, [1, 3, 2]);
reconVol_gas_highSNR = permute(reconVol_gas_highSNR, [1, 3, 2]);
reconVol_gas_highreso = permute(reconVol_gas_highreso, [1, 3, 2]);

% Inconveniently, these historic sequences are not collected in same
% orientation

[Shifted_Image, Ang, B0PhaseMap_Tot] = SinglePointDixonV2(...
        reconVol_dis, -rbc2bar, reconVol_gas_highSNR, logical(protonMask));
[Shifted_Tot_Image, Tot_Ang, B0PhaseMap] = SinglePointDixonV2(...
        reconVol_dis_tot, -rbc2bar_all.Tot, reconVol_gas_highSNR, logical(protonMask));
[Shifted_High_Image, High_Ang, B0PhaseMap_High] = SinglePointDixonV2(...
        reconVol_dis_highkey, -rbc2bar_all.High, reconVol_gas_highSNR, logical(protonMask));
[Shifted_Low_Image, Low_Ang, B0PhaseMap_Low] = SinglePointDixonV2(...
        reconVol_dis_lowkey, -rbc2bar_all.Low, reconVol_gas_highSNR, logical(protonMask));

[RBC_Osc, oscillation_stats, RBC_PctRBC_Key_Diff_Fig, Binning_Montage_Fig, Histograms] = ...
    analyze_display_detrend_keyhole_dissolvedPhase(Shifted_High_Image, Shifted_Low_Image, ...
        Shifted_Tot_Image, reconVol_gas_highreso, protonMask, rbc2bar_all, 'Siemens', Subject, params.filepath);

%% Get SNR Values
RBC_Tot = abs(real(Shifted_Tot_Image));
RBC_High = abs(real(Shifted_High_Image));
RBC_Low = abs(real(Shifted_Low_Image));
Bar_Tot = abs(imag(Shifted_Tot_Image));
Bar_High = abs(imag(Shifted_High_Image));
Bar_Low = abs(imag(Shifted_Low_Image));
Dis_Tot = abs(Shifted_Tot_Image);
Dis_High = abs(Shifted_High_Image);
Dis_Low = abs(Shifted_Low_Image);
Gas_Image = reconVol_gas_highSNR;

% Get image SNR
[RBC_High_SNR, RBC_High_SNR_Rayleigh, RBC_High_Signal, RBC_High_Noise] = ...
    Tools.imageSNR(RBC_High, protonMask, 8, 0.75*8^3);
[RBC_Low_SNR, RBC_Low_SNR_Rayleigh, RBC_Low_Signal, RBC_Low_Noise] = ...
    Tools.imageSNR(RBC_Low, protonMask, 8, 0.75*8^3);
[RBC_SNR, RBC_SNR_Rayleigh, RBC_Signal, RBC_Noise] = ...
    Tools.imageSNR(RBC_Tot, protonMask, 8, 0.75*8^3);

[Bar_High_SNR, Bar_High_SNR_Rayleigh, Bar_High_Signal, Bar_High_Noise] = ...
    Tools.imageSNR(Bar_High, protonMask, 8, 0.75*8^3);
[Bar_Low_SNR, Bar_Low_SNR_Rayleigh, Bar_Low_Signal, Bar_Low_Noise] = ...
    Tools.imageSNR(Bar_Low, protonMask, 8, 0.75*8^3);
[Bar_SNR, Bar_SNR_Rayleigh, Bar_Signal, Bar_Noise] = ...
    Tools.imageSNR(Bar_Tot, protonMask, 8, 0.75*8^3);

[Dis_High_SNR, Dis_High_SNR_Rayleigh, Dis_High_Signal, Dis_High_Noise] = ...
    Tools.imageSNR(abs(Shifted_High_Image), protonMask, 8, 0.75*8^3);
[Dis_Low_SNR, Dis_Low_SNR_Rayleigh, Dis_Low_Signal, Dis_Low_Noise] = ...
    Tools.imageSNR(abs(Shifted_Low_Image), protonMask, 8, 0.75*8^3);
[Dis_SNR, Dis_SNR_Rayleigh, Dis_Signal, Dis_Noise] = ...
    Tools.imageSNR(abs(Shifted_Tot_Image), protonMask, 8, 0.75*8^3);

[Gas_SNR, Gas_SNR_Rayleigh, Gas_Signal, Gas_Noise] = ...
    Tools.imageSNR(abs(Gas_Image), protonMask, 8, 0.75*8^3);

[recon_dis_SNR, recon_dis_SNR_Rayleigh, recon_dis_Signal, recon_dis_Noise] = ...
    Tools.imageSNR(abs(reconVol_dis), protonMask, 8, 0.75*8^3);

RBC_Tot_Value = RBC_Tot(protonMask == 1);
Bar_Tot_Value = Bar_Tot(protonMask == 1);
RBC_High_Value = RBC_High(protonMask == 1);
Bar_High_Value = Bar_High(protonMask == 1);
RBC_Low_Value = RBC_Low(protonMask == 1);
Bar_Low_Value = Bar_Low(protonMask == 1);

check_tot_shift = mean(RBC_Tot_Value) / mean(Bar_Tot_Value);
check_high_shift = mean(RBC_High_Value) / mean(Bar_High_Value);
check_low_shift = mean(RBC_Low_Value) ./ mean(Bar_Low_Value);

%% Make Gas and Mask FiguresGas
maskFigure = figure('Name', 'maskFigure');
ImSize = recon_size(1);
clf
if ImSize == 128
    start = 32;
    slices = start:4:start + 4 * 15;
else
    start = 16;
    slices = start:2:start + 2 * 15;
end

gas_fig = abs(reconVol_gas_highreso);

[ha, ~] = tight_subplot(2, length(slices)/2, 0, 0, [0.06, 0.06]);
red = cat(3, ones(ImSize, ImSize), ones(ImSize, ImSize)*.2, ones(ImSize, ImSize)*.2);
for idx = 1:length(slices)
    axes(ha(idx));
    imshow(squeeze(gas_fig(:, slices(idx), :)), [0, max(gas_fig(:))])
    hold on
    B = bwboundaries(squeeze(protonMask(:, slices(idx), :)));
    visboundaries(B, 'Linewidth', 0.25)
    % h = imshow(red);
    % hold off
    % set(h,'AlphaData',squeeze(protonMask(:,slices(idx),:))*0.2)
end
set(gcf, 'Position', [233, 553, 1454, 413])

%% RBC Figure
RBC_Tot_Figure = figure('Name', 'RBC_Tot_Figure');
set(RBC_Tot_Figure, 'color', 'white', 'Units', 'inches', 'Position', [0.25, 0.25, 12, 7])
clf

% RBC_Tot = (real(Shifted_Low_Image));%-real(Shifted_Low_Image));
RBC_fig = abs(real(Shifted_Tot_Image));
[ind_start, ind_inter, ind_end] = decideStartInterval(mask_reg);
montage(permute(RBC_fig(:, int16(linspace(ind_start, ind_end, 16)), :), [1, 3, 2]), 'size', [2, 8], 'DisplayRange', [0, max(RBC_fig(:))])

title('RBC Image')

%% Oscillation Figure
RBC_Oscs_Figure = figure('Name', 'RBC_Oscs_Figure');
set(RBC_Oscs_Figure, 'color', 'white', 'Units', 'inches', 'Position', [0.25, 0.25, 12, 7])
clf

RBC_fig = (RBC_High - RBC_Low) ./ mean(RBC_Tot(protonMask == 1)) * 100;
[ind_start, ind_inter, ind_end] = decideStartInterval(mask_reg);
montage(permute(RBC_fig(:, int16(linspace(ind_start, ind_end, 16)), :), [1, 3, 2]), 'size', [2, 8], 'DisplayRange', [-10, max(RBC_fig(:))]);

colorbar
title('RBC (High - Low)/Mean Image')

%% PPT Report

%Start new presentation
% isOpen  = exportToPPTX();
% if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
%     exportToPPTX('close');
% end
pptx = exportToPPTX('', 'Dimensions', [16, 9], ...
    'Title', ReportTitle, ...
    'Author', 'CPIR @ CCHMC');

%Add slides
pptx.addSlide(); %Title Slide
pptx.addTextbox(sprintf(Subject), ...
    'Position', [0, 0, 16, 4.5], ...
    'VerticalAlignment', 'bottom', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 72, ...
    'FontWeight', 'bold');
pptx.addTextbox(sprintf(['Processing Date: ', datestr(date, 29)]), ...
    'Position', [0, 4.5, 16, 4.5], ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 36);

pptx.addSlide(); %FIDs
pptx.addPicture(Shift_Raw_FIDs_Fig);
pptx.addTextbox(sprintf('k0 for FIDs'));

pptx.addSlide(); %Binning Figure
pptx.addPicture(Smooth_Detrend_Fig);
pptx.addTextbox(sprintf('Binning Workflow'));

pptx.addSlide();
pptx.addTextbox(sprintf('Image Characteristics'));

tableData = { ...
    'Gas SNR', Gas_SNR; ...
    'Total Dissolved SNR', Dis_SNR; ...
    'High Dissolved SNR', Dis_High_SNR; ...
    'Low Dissolved SNR', Dis_Low_SNR; ...
    'Total Barrier SNR', Bar_SNR; ...
    'High Barrier SNR', Bar_High_SNR; ...
    'Low Barrier SNR', Bar_Low_SNR; ...
    'Total RBC SNR', RBC_SNR; ...
    'High RBC SNR', RBC_High_SNR; ...
    'Low RBC SNR', RBC_Low_SNR; ...
    'GX Recon Dissolved SNR', recon_dis_SNR; ...
    'Total RBC/Barrier', check_tot_shift; ...
    'High RBC/Barrier', check_high_shift; ...
    'Low RBC/Barrier', check_low_shift; ...
    'Views Removed (Front)', N_Begin_Cutoff; ...
    'Views Removed (Back)', N_End_Cutoff; ...
    'Key Radius (nPts)', Key_Radius; ...
    'GX Sequence', data_type;};

pptx.addTable(tableData, 'Position', [2, 1, 8, 4], ...
    'Vert', 'middle', 'Horiz', 'center', 'FontSize', 13, ...
    'ColumnWidth', [0.5, 0.5], 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'LineWidth', 2);

pptx.addSlide(); %Mask
pptx.addPicture(maskFigure);
pptx.addTextbox(sprintf('Proton Mask with Gas Image'));

pptx.addSlide(); %RBC Montage
pptx.addPicture(RBC_Tot_Figure);
pptx.addTextbox(sprintf('RBC Image'));

pptx.addSlide(); %RBC Montage
pptx.addPicture(RBC_Oscs_Figure);
pptx.addTextbox(sprintf('RBC (High - Low)/Mean Image'))

pptx.addSlide(); %RBC
pptx.addPicture(RBC_PctRBC_Key_Diff_Fig);
pptx.addTextbox(sprintf('RBC Difference'));

pptx.addSlide(); %RBC Binned
pptx.addPicture(Binning_Montage_Fig);
pptx.addTextbox(sprintf(['Binned RBC Oscillations Images for Subject ', Subject]));

pptx.addSlide(); %Histograms
pptx.addPicture(Histograms);
pptx.addTextbox(sprintf(['Whole Lung Histograms for Subject ', Subject]));

pptx.addSlide();
pptx.addTextbox(sprintf('Oscillation Characteristics'));

tableData = { ...
    'Subject', oscillation_stats.Subject; ...
    'Date', oscillation_stats.Date; ...
    'Mean oscillation Percent', oscillation_stats.Pct_Mean_Mean; ...
    'Mean oscillation STD', oscillation_stats.Pct_Mean_STD; ...
    'Mean oscillation CV', oscillation_stats.Pct_Mean_CV; ...
    'MeanRBCMeasure', oscillation_stats.MeanRBCMeasure,; ...
    'Total RBC SNR', RBC_SNR; ...
    '% Low Oscillation (Red)', 100 * oscillation_stats.Bin1; ...
    '% Low Oscillation (Red-Orange)', 100 * (oscillation_stats.Bin1 + oscillation_stats.Bin2); ...
    '% High Oscillation (Purple 1-3)', 100 * (oscillation_stats.Bin6 + oscillation_stats.Bin7 + oscillation_stats.Bin8); ...
    '% Low SNR (analysis ignored)', 100 * (oscillation_stats.Bin9); ...
    'RBC:membrane', rbc2bar, ...
    };

pptx.addTable(tableData, 'Position', [2, 1, 10, 6], ...
    'Vert', 'middle', 'Horiz', 'center', 'FontSize', 9, ...
    'ColumnWidth', [0.5, 0.5], 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'LineWidth', 1);

% summary slide
pptx.addSlide();
pptx.addTextbox(sprintf('RBC Oscillation Analysis - Exploratory'));
pptx.addPicture(Binning_Montage_Fig, 'Position', [1.5, 5, 14, 4]);
pptx.addTable(tableData, 'Position', [.5, .5, 4, 4], ...
    'Vert', 'middle', 'Horiz', 'center', 'FontSize', 11, ...
    'ColumnWidth', [0.5, 0.5], 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'LineWidth', 2);
pptx.addPicture(Histograms, 'Position', [5, .5, 11, 4]);
pptx.addTextbox(sprintf('Binned Oscillations'), 'FontSize', 14, ...
    'Position', [7, 4.75, 2, 1]);

cbar_Fig = figure('Name', 'Colorbar');
[cbar, map, alphachannel] = imread('img/colorbar_binnedosc_white.png');
h = imshow(cbar);
colormap(map);
set(cbar_Fig, 'Color', 'w');
pptx.addPicture(tightfig(cbar_Fig), 'Position', [.5, 5, 1.25, 3]);
%Save presentation and close presentation -- overwrite file if it already exists
newFile = pptx.save(fullfile('Summaries_Long/', ReportTitle));


% Create a short powerpoint

pptx = exportToPPTX('', 'Dimensions', [16, 9], ...
    'Title', ReportTitle, ...
    'Author', 'CPIR @ CCHMC');
pptx.addSlide();
pptx.addTextbox(sprintf('RBC Oscillation Analysis - Exploratory'));
pptx.addPicture(Binning_Montage_Fig, 'Position', [1.5, 5, 14, 4]);
pptx.addTable(tableData, 'Position', [.5, .5, 4, 4], ...
    'Vert', 'middle', 'Horiz', 'center', 'FontSize', 13, ...
    'ColumnWidth', [0.5, 0.5], 'BackgroundColor', 'w', ...
    'EdgeColor', 'k', 'LineWidth', 2);
pptx.addPicture(Histograms, 'Position', [5, .5, 11, 4]);
pptx.addTextbox(sprintf('Binned Oscillations'), 'FontSize', 14, ...
    'Position', [7, 4.75, 2, 1]);

cbar_Fig = figure('Name', 'Colorbar');
[cbar, map, alphachannel] = imread('img/colorbar_binnedosc_white.png');
h = imshow(cbar);
colormap(map);
set(cbar_Fig, 'Color', 'w');
pptx.addPicture(tightfig(cbar_Fig), 'Position', [.5, 5.8, 1.25, 3])
newFile = pptx.save(fullfile('Summaries/', ReportTitle));

%% Save important values to spreadsheet
% Summary Oscillation Values
disp('Saving Values of Interest to Summary Spreadsheet...')
[~, ~, ExcelFileRaw] = xlsread('Summaries/Keyhole_Summary_Auto.xlsx');
NextLine = size(ExcelFileRaw, 1) + 1;
NewLine = {Subject, datestr(date, 29), data_type, oscillation_stats.MeanRBCMeasure, oscillation_stats.Pct_Mean_Mean, ...
    oscillation_stats.Pct_Mean_STD, oscillation_stats.Pct_Mean_CV, oscillation_stats.RBC_Low_SNR, ...
    oscillation_stats.RBC_High_SNR, oscillation_stats.RBC_Tot_SNR, ...
    oscillation_stats.Bin1 + oscillation_stats.Bin2, ...
    oscillation_stats.Bin6 + oscillation_stats.Bin7 + oscillation_stats.Bin8, ...
    oscillation_stats.Bin1, oscillation_stats.Bin2, oscillation_stats.Bin3, oscillation_stats.Bin4, ...
    oscillation_stats.Bin5, oscillation_stats.Bin6, oscillation_stats.Bin7, oscillation_stats.Bin8, ...
    oscillation_stats.Bin9, rbc2bar_all.Tot, rbc2bar_all.High, rbc2bar_all.Low, Key_Radius, ...
    };
xlRange = ['A', num2str(NextLine)];
xlswrite('Summaries/Keyhole_Summary_Auto.xlsx', NewLine, 1, xlRange);

%% Save variables as mat

save(fullfile(subjFolder, strcat(Subject, '_osc.mat')), 'mask_reg', 'RBC_Osc');
close all;

end
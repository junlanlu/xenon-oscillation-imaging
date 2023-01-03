function [Mean_Pct_Diff, oscillation_stats, ...
    RBC_PctRBC_Key_Diff_Fig, Binning_Montage_Fig, Histogram_Fig] = ...
    analyze_display_detrend_keyhole_dissolvedPhase_elly(High_Image, ...
    Low_Image, Tot_Image, Gas_Image, Mask, rbc2bar_all, ScannerType, Subject, path)

%Function to analyze and Display images obtained using keyhole dissolved
%phase imaging.
%Pass Phase-shifted images so that RBC images are in the real, and Barrier
%are in the -Imaginary
%
% Updates:
% * changed code for color bin routine
% * pixels with low RBC SNR now grey instead of white

%% Definitions
NewCMap = summer;
NewCMap(1, :) = [0, 0, 0];

%% Specify Healthy Cohort Means and SDs for Oscillations, RBC, Barrier,
% and RBC/Barrier

%Sometimes the absolute signal intensity is super low... scale to get to
%something a little more friendly to work with
if mean(abs(Gas_Image(Mask == 1))) < 1e-4
    Gas_Image = Gas_Image * 1e8;
    Low_Image = Low_Image * 1e8;
    High_Image = High_Image * 1e8;
    Tot_Image = Tot_Image * 1e8;
end

if contains(ScannerType, 'GE')

    healthy_mean_Osc = 8.4997;
    healthy_std_Osc = 11.2616;

    healthy_mean_RBCBar = 0.4835;
    healthy_std_RBCBar = 0.2457;
    %Load in Healthy Cohort Distribution for Later QQ Plots
    load('GE_Healthy.mat');

    Mean_Disp_Range = [healthy_mean_Osc - 2 * healthy_std_Osc, ...
        healthy_mean_Osc + 4 * healthy_std_Osc];
else
    % - These values are for Detrending Raw Data using RBC Fit
    % rather than dissolved phase fit (Probably not how we want to do it)
    % healthy_mean_Mean = 9.4912;
    % healthy_std_Mean = 8.9182;
    healthy_mean_Osc = 9.8104;
    healthy_std_Osc = 9.0112;

    healthy_mean_RBCBar = 0.5463;
    healthy_std_RBCBar = 0.2507;

    Mean_Disp_Range = [healthy_mean_Osc - 2 * healthy_std_Osc, ...
        healthy_mean_Osc + 4 * healthy_std_Osc];
end

%% Separate into image components
RBC_High = real(High_Image);
if mean(RBC_High(Mask == 1)) < 0
    RBC_High = -RBC_High;
end
RBC_Low = real(Low_Image);
if mean(RBC_Low(Mask == 1)) < 0
    RBC_Low = -RBC_Low;
end
RBC_Tot = real(Tot_Image);
if mean(RBC_Tot(Mask == 1)) < 0
    RBC_Tot = -RBC_Tot;
end

Bar_Tot = imag(Tot_Image);
if mean(Bar_Tot(Mask == 1)) < 0
    Bar_Tot = -Bar_Tot;
end

%Get SNR
RBC_High_SNR = Tools.imageSNR(RBC_High, Mask, 8, 0.75*8^3);
RBC_Low_SNR = Tools.imageSNR(RBC_Low, Mask, 8, 0.75*8^3);
[RBC_SNR, ~, ~, RBC_Noise] = Tools.imageSNR(RBC_Tot, Mask, 8, 0.75*8^3);

%Create a Mask accounting for RBC Defects - use a cutoff of 1.5?
RBC_Mask = Mask;
RBC_Mask(RBC_Tot < 0.5*RBC_Noise) = 0;

%% Display High and Low Key Images (RBC) and Pct Difference Map
First_Slice = floor(length(RBC_High)*.25);
Last_Slice = floor(length(RBC_High)*.75);
%Let's display six slices
step = floor((Last_Slice - First_Slice)/5);
Slices = First_Slice:step:(First_Slice + 5 * step);

RBC_PctRBC_Key_Diff_Fig = figure('Name', 'Difference between High and Low Key ');
set(RBC_PctRBC_Key_Diff_Fig, 'color', 'white', 'Units', 'inches', 'Position', [0.25, 0.25, 12, 7])
[ha, ~] = tight_subplot(4, length(Slices), 0.01, 0.01, [0.05, 0.01]);
maxvox = max(max([abs(RBC_High(:)), abs(RBC_Low(:))]));

for slice = 1:length(Slices)
    axes(ha(slice));
    imagesc(((((squeeze(RBC_High(:, Slices(slice), :)))))))
    colormap(gray)
    caxis manual
    caxis([0, maxvox])
    axis square
    axis off
end
for slice = 1:length(Slices)
    axes(ha(slice+length(Slices)));
    imagesc(((((squeeze(RBC_Low(:, Slices(slice), :)))))))
    colormap(gray)
    caxis manual
    caxis([0, maxvox])
    axis square
    axis off
end

axes(ha(1));
ax = gca;
ax.FontSize = 12;
ylabel('RBC High Key')
ax.YLabel.Visible = 'on';

axes(ha(1+1*length(Slices)));
ax = gca;
ax.FontSize = 12;
ylabel('RBC Low Key')
ax.YLabel.Visible = 'on';

Mean_Pct_Diff = (RBC_High - RBC_Low) / mean(RBC_Tot(RBC_Mask == 1)) * 100;
Mean_Pct_Diff = Mean_Pct_Diff .* RBC_Mask;
Tmp_Mean_Pct_Diff = Mean_Pct_Diff;
Tmp_Mean_Pct_Diff(Mean_Pct_Diff == 0) = -500;

for slice = 1:length(Slices)
    axes(ha(slice+2*length(Slices)));
    imagesc(((((squeeze(Tmp_Mean_Pct_Diff(:, Slices(slice), :)))))))
    colormap(ha(slice+2*length(Slices)), NewCMap)
    caxis manual
    caxis(Mean_Disp_Range)
    axis square
    axis off
end

axes(ha(1+2*length(Slices)));
c = colorbar('south');
c.Color = 'w';
axes(ha(1+2*length(Slices)));
ax = gca;
ax.FontSize = 12;
ylabel('RBC Oscillation')
ax.YLabel.Visible = 'on';

%% Bin Oscillation Image
Tools.GX_defineColormaps;
MeanBinMap = Tools.eight_bin_image(Mean_Pct_Diff, Mask, healthy_thresh_Osc);
MeanBinMap(Mask == 1 & RBC_Mask == 0) = 9;

%Nine Bin Map for RBC Oscillations
EightBinMap_Osc = [0, 0, 0; ...
    1, 0, 0; ... %Red
    1, 0.7143, 0; ... %Orange
    0.4, 0.7, 0.4; ... %Green1
    0, 1, 0; ... %Green2
    184 / 255, 226 / 255, 145 / 255; ... %Green3
    243 / 255, 205 / 255, 213 / 255; ... %Light Pink
    225 / 255, 129 / 255, 162 / 255; ... %Med Pink
    197 / 255, 27 / 255, 125 / 255; ... %Dark Pink
    0.33, 0.33, 0.33; ... % Gray (For regions of RBC Defect)
    ];

for slice = 1:length(Slices)
    axes(ha(slice+3*length(Slices)));
    imagesc(((((squeeze(MeanBinMap(:, Slices(slice), :)))))))
    colormap(ha(slice+3*length(Slices)), EightBinMap_Osc)
    caxis manual
    caxis([0, 9])
    axis square
    axis off
end
colorbar(ha(3*length(Slices)+1), 'north', 'Color', [1, 1, 1])
axes(ha(1+3*length(Slices)));
ax = gca;
ax.FontSize = 22;
ylabel('Binned Oscillations')
ax.YLabel.Visible = 'on';

%% Image with Binned RBC, Bar, and Osc

Binning_Montage_Fig = figure('Name', 'RBC Binning Montage');

set(Binning_Montage_Fig, 'color', 'white', 'Units', 'inches', 'Position', ...
    [0.5, 0.5, 12, 7])

RBC_Mask = Mask;
RBC_Mask(RBC_Tot < 1.5*RBC_Noise) = 0;

Mean_Pct_Diff = (RBC_High - RBC_Low) / mean(RBC_Tot(RBC_Mask == 1)) * 100;
Mean_Pct_Diff = Mean_Pct_Diff .* RBC_Mask;

[ind_start, ~, ind_end] = decideStartInterval(permute(Mask, [1, 3, 2]));
montage(permute(MeanBinMap(:, int16(linspace(ind_start, ind_end, 16)), :), ...
    [1, 3, 2]), 'size', [2, 8], 'DisplayRange', [-0.5, 9.5]);

colormap(EightBinMap_Osc);

% Gas_Image(axial, coronal, sagittal)

%% Get Means and STDs for both measures of RBC oscillation pct.

outlier_thresh = [-35, 70];

Mean_Pct_Diff_Data = Mean_Pct_Diff;
Mean_Pct_Diff_Data(RBC_Mask == 0) = [];
Mean_Pct_Diff_Data(Mean_Pct_Diff_Data < outlier_thresh(1)) = [];
Mean_Pct_Diff_Data(Mean_Pct_Diff_Data > outlier_thresh(2)) = [];

Pct_Mean_Mean = mean(Mean_Pct_Diff_Data);
Pct_Mean_STD = std(Mean_Pct_Diff_Data);
Pct_Mean_CV = Pct_Mean_STD / Pct_Mean_Mean;

mean_High_data = RBC_High(Mask == 1);
mean_Low_data = RBC_Low(Mask == 1);
mean_Tot_data = RBC_Tot(Mask == 1);

mean_High_data(mean_High_data < 0) = [];
mean_Low_data(mean_Low_data < 0) = [];
mean_Tot_data(mean_Tot_data < 0) = [];

mean_High_data(isoutlier(mean_High_data)) = [];
mean_Low_data(isoutlier(mean_Low_data)) = [];
mean_Tot_data(isoutlier(mean_Tot_data)) = [];

mean_High = mean(mean_High_data);
mean_Low = mean(mean_Low_data);
mean_Tot = mean(mean_Tot_data);

MeanRBCMeasure = ((mean_High - mean_Low) / mean_Tot) * 100;

%% Get Histograms for RBC Differences
Mean_Pct_Diff = (RBC_High - RBC_Low) / mean(RBC_Tot(RBC_Mask == 1)) * 100;
Mean_Pct_Diff = Mean_Pct_Diff .* RBC_Mask;
raw_oscillations = Mean_Pct_Diff;
raw_oscillations(raw_oscillations == 0) = [];

if contains(Subject, 'GE')
    Mean_Edges = linspace(-40, 80, 100);
else
    Mean_Edges = linspace(-30, 60, 100);

end

Healthy_Osc_Dist = normpdf(Mean_Edges, healthy_thresh_Osc(3), 9.0);
Healthy_Osc_Dist = Healthy_Osc_Dist / sum(Healthy_Osc_Dist(:));

Histogram_Fig = figure('Name', 'Histograms');
set(Histogram_Fig, 'color', 'white', 'Units', 'inches', 'Position', ...
    [0.25, 0.25, 8, 3])

%RBC/Barrier
subplot(1, 2, 1)
RBCBarrierEdges = linspace(0, healthy_mean_RBCBar+10*healthy_std_RBCBar, 100);
HeatlhyRBCBarrierDist = normpdf(RBCBarrierEdges, healthy_mean_RBCBar, healthy_std_RBCBar);
NormFactor = sum(HeatlhyRBCBarrierDist(:));
HeatlhyRBCBarrierDist = HeatlhyRBCBarrierDist ./ NormFactor;
RBC_Bar_tmp = RBC_Tot ./ Bar_Tot;
RBC_Bar_tmp(Mask == 0) = [];
hold on
histogram(RBC_Bar_tmp, RBCBarrierEdges, 'Normalization', 'probability', 'FaceColor', [0.75, 0, 0.75], 'FaceAlpha', 1);
plot(RBCBarrierEdges, HeatlhyRBCBarrierDist, 'k--');
axis([0, RBCBarrierEdges(end-1) + RBCBarrierEdges(end-1) - RBCBarrierEdges(end-2), 0, inf])
legend('Subject', 'Healthy', 'Location', 'NorthEast');
title('RBC:Barrier Ratio Histogram')

%Oscillations Histogram
subplot(1, 2, 2)
hold on
histogram(raw_oscillations, Mean_Edges, 'Normalization', 'probability', 'FaceColor', [.1, 0.67, 0.67], 'FaceAlpha', 1);
plot(Mean_Edges, Healthy_Osc_Dist, 'k--');
ylabel('Bin Count')
xlabel('Percent Oscillation (of mean non-keyholed RBC)')
legend('Subject', 'Healthy', 'Location', 'NorthEast');
hold off
title('Cardio-Pulmonary Oscillations Histogram');
pbaspect([2, 1, 1]);

%Find the percentage of voxels in each particular bin
NumVox = length(find(Mask == 1));
Bin1 = length(find(MeanBinMap == 1)) / NumVox;
Bin2 = length(find(MeanBinMap == 2)) / NumVox;
Bin3 = length(find(MeanBinMap == 3)) / NumVox;
Bin4 = length(find(MeanBinMap == 4)) / NumVox;
Bin5 = length(find(MeanBinMap == 5)) / NumVox;
Bin6 = length(find(MeanBinMap == 6)) / NumVox;
Bin7 = length(find(MeanBinMap == 7)) / NumVox;
Bin8 = length(find(MeanBinMap == 8)) / NumVox;
Bin9 = length(find(MeanBinMap == 9)) / NumVox;
% Output a struct
oscillation_stats = struct('Subject', Subject, 'Date', datestr(date, 29), 'MeanRBCMeasure', MeanRBCMeasure, 'Pct_Mean_Mean', ...
    Pct_Mean_Mean, 'Pct_Mean_STD', Pct_Mean_STD, 'Pct_Mean_CV', Pct_Mean_CV, 'RBC_Low_SNR', RBC_Low_SNR, 'RBC_High_SNR', RBC_High_SNR, ... .
    'RBC_Tot_SNR', RBC_SNR, 'Bin1', Bin1, 'Bin2', Bin2, 'Bin3', Bin3, 'Bin4', Bin4, 'Bin5', Bin5, 'Bin6', Bin6', 'Bin7', Bin7, 'Bin8', Bin8, 'Bin9', Bin9);

%% Save workspace
% save(fullfile(path,'Detrended Gas Exchange Keyhole Workspace V3.mat'));

%% Save 3D RGB nifti images

Tools.GX_defineColormaps;
MeanBinMap3 = Tools.Bin2RGB(MeanBinMap, index2color_eightbin_map);
Tools.save3DRGB2nii(MeanBinMap3, fullfile(path, append('osc_Sub', Subject, '.nii')));
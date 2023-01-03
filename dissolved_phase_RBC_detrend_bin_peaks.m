function [High_Key, Low_Key, Dis2GasFID, RBC2Bar, Shift_Raw_FIDs_Fig, Smooth_Detrend_Fig, Key_Radius] = dissolved_phase_RBC_detrend_bin_peaks(Traj, GasFID, DisFID, TR, RBC2Barrier, params)
%Function to Bin Dissolved Phase Data based on RBC oscillations. This
%method is different than the original paper. It will apply a band pass
%filter and use a peak finding algorithm rather than looking at threshold
%Arguments: Traj - trajectories in 3xNPtsxNProj format
%GasFID and DisFID - FIDs for Gas and Dissolved Phase in NPtsxNProj Matrix
%TR - Repetition time (For making nicer figures)
%RBC2Barrier - RBC to Barrier Ratio (for separating RBC and Barrier)

%Returns:
%High_Key - KSpace with keyholing for High RBC Bin (Using Detrended Data)
%Low_Key - KSpace with keyholing for Low RBC Bin (Using Detrended Data)
%RBC2Bar - Structure holding High, Low, Total RBC/Barrier Ratios and
%Ratios of high to low total signal
%Shift_Raw_FIDs_Fig - Figure Handle for figure showing raw k0 points before
%and after phase shift
%Smooth_Detrend_Fig - Figure Handle for figure showing binning workflow

%% Get Number of Projections
NProj = size(DisFID, 2);

%% Perform Phase Shift

%Empirically, this is the way to properly phase shift data on each platform
switch params.ScannerType
    case 'GE'
        [Sh_DisFID, RBC2Bar.PhaseShift] = SinglePointDixonV2_FID(DisFID, -RBC2Barrier, GasFID);
    case 'Siemens'
        [Sh_DisFID, RBC2Bar.PhaseShift] = SinglePointDixonV2_FID_elly((DisFID), RBC2Barrier, GasFID);
    case 'Philips'
        [Sh_DisFID, RBC2Bar.PhaseShift] = SinglePointDixonV2_FID(DisFID, -RBC2Barrier, GasFID);
end
%Store Phase Shift Angle - I used this for some other stuff... probably not
%needed, but not hurting anything.
RBC2Bar.meanAngle = mean(angle(Sh_DisFID(1, :)));
RBC2Bar.InitAngle = mean(angle(DisFID(1, :)));
RBC2Bar.GasAngle = mean(angle(GasFID(1, :)));
RBC2Bar.PhaseChange = RBC2Bar.meanAngle - RBC2Bar.InitAngle + RBC2Bar.GasAngle;
Time_Axis = (0:(NProj - 1)) * TR;

%% Display k0 points

Shift_Raw_FIDs_Fig = figure('Name', 'Dissolved Phase Mag Real Imaginary k0 Before and After phase shift');

%Not sure if TR is ever specified in Matt's WIP code - if not, make sure it
%gets set and use that variable name here:
%Give RBC, Barrier, and Total new names to cut down on typing
Raw_Mag_k0 = abs(DisFID(1, :));
Raw_Real_k0 = real(DisFID(1, :));
Raw_Imag_k0 = imag(DisFID(1, :));
Raw_GasMag_k0 = abs(GasFID(1, :));

Sh_Mag_k0 = abs(Sh_DisFID(1, :));
Sh_Real_k0 = real(Sh_DisFID(1, :));
Sh_Imag_k0 = imag(Sh_DisFID(1, :));

subplot(3, 2, 1);
plot(Time_Axis, Raw_Mag_k0, '.-', 'Color', [.49, .13, .53])
title('Magnitude Dissolved Phase k0 (Before Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])
subplot(3, 2, 3);
plot(Time_Axis, Raw_Real_k0, 'r.-')
title('Real Dissolved Phase k0 (Before Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])

subplot(3, 2, 5);
plot(Time_Axis, Raw_Imag_k0, 'b.-')
title('Imaginary Dissolved Phase k0 (Before Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])

subplot(3, 2, 2);
plot(Time_Axis, Sh_Mag_k0, '.-', 'Color', [.49, .13, .53])
title('Magnitude Dissolved Phase k0 (After Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])

subplot(3, 2, 4);
plot(Time_Axis, Sh_Real_k0, 'r.-')
title('Real Dissolved Phase k0 (After Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])

subplot(3, 2, 6);
plot(Time_Axis, Sh_Imag_k0, 'b.-')
title('Imaginary Dissolved Phase k0 (After Phase Shift)')
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')
ax = gca;
ax.FontSize = 12;
xlim([0, 14.5])

set(Shift_Raw_FIDs_Fig, 'color', 'white', 'Units', 'inches', 'Position', [0.25, 1, 12, 7])

%% Smooth and Detrend Data
%Assume Oscillation rate of approx 1 Hz... Thus, there would be NProj*TR
%oscillations. That means the number of points per osc is
%NProj/NProj*TR = 1/TR. Then, we want about 1/5th that number, so...
Sm_Window = floor(1/TR/5);
Samp_Rate = 1 / TR;

%Make sure that we have positive RBC and Barrier Values
Mag_Raw_Dis = Sh_Mag_k0;
if mean(Sh_Real_k0) > 0
    RBC_Raw = Sh_Real_k0;
else
    RBC_Raw = -Sh_Real_k0;
end
if mean(Sh_Imag_k0) > 0
    Bar_Raw = Sh_Imag_k0;
else
    Bar_Raw = -Sh_Imag_k0;
end

%To make my life easier, normalize the RBC Data by dividing by k0
RBC2Gas = Sh_Real_k0 ./ Raw_GasMag_k0;

%Smooth RBC Data using smooth function
RBC_Smoothed_avg = smooth(Sh_Real_k0, Sm_Window);
% Smooth RBC Data using bandpass filter
RBC_Smoothed_BP = bandpass(RBC2Gas, [0.5, 2.5], Samp_Rate);
% RBC_Smoothed = RBC_Smoothed + mean(RBC2Gas);
%Display a detrending and binning workflow figure
Smooth_Detrend_Fig = figure('Name', "Smoothing and Detrending Workflow");
set(Smooth_Detrend_Fig, 'color', 'white', 'Units', 'inches', 'Position', [0.25, 1, 15, 7])

% plot the original decomposed RBC signal
subplot(2, 4, 1)
plot(Time_Axis, Sh_Real_k0, 'k.-');
title('1. RBC Raw k0 Points')
axis square;
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

% plot the RBC2Gas
subplot(2, 4, 2)
plot(Time_Axis, RBC2Gas, 'k.-', 'LineWidth', 2);
title('2. Detrended RBC k0 Points')
axis square
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

% plot the bandpass filtered RBC2Gas
subplot(2, 4, 3)
plot(Time_Axis, RBC_Smoothed_BP, 'k.-');
hold on
title('Smoothed using bandpass filter')
axis square
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

FFT_RBC = abs(fftshift(fft(RBC_Smoothed_BP)));
SampF = 1 / TR;
%Get x frequency space
Freq_Axis = SampF * ((-(NProj / 2) + 1):(NProj / 2)) / NProj;

subplot(2, 4, 4)
plot(Freq_Axis(round(((NProj / 2) + 1):end)), FFT_RBC(round(((NProj / 2) + 1):end)), 'k');
xlim([0, 8])
title({'5. Fourier Transform', 'of smoothed and detrended RBC'})
axis square
xlabel('Frequency (Hz)')
ylabel('Intensity (Arb. Units)')

%% Compare to smooth, detrend, stretch
Detrend_Fit = fit(Time_Axis', RBC_Smoothed_avg, 'exp2'); %,'StartPoint',[.2 -.5 .7 -.1],'Upper',[1 0 1 1],'Lower',[0 -10 0 -10]);
Detrend_Pts = Detrend_Fit.a * exp(Detrend_Fit.b*Time_Axis) + Detrend_Fit.c * exp(Detrend_Fit.d*Time_Axis);
Detrend_Pts = Detrend_Pts';
Detrend1_Smooth_RBC = (RBC_Smoothed_avg + 1) ./ (Detrend_Pts + 1) - 1;
polynomial_fitting = polyfit(Time_Axis', Detrend1_Smooth_RBC, 10);
sub_pts = polyval(polynomial_fitting, Time_Axis');;
Detrend2_Smooth_RBC = Detrend1_Smooth_RBC - sub_pts;
RBCStretch = Tools.stretchdata(Detrend2_Smooth_RBC);
Threshold = 0.6;
High_Bin_Index_stretch = find(RBCStretch > Threshold);
Low_Bin_Index_stretch = find(RBCStretch < -Threshold);


Orange_Color = [248 / 256, 136 / 256, 12 / 156];

subplot(2, 4, 5)
plot(Time_Axis, RBCStretch, 'k.-', 'MarkerSize', 8)
hold on
plot(Time_Axis(High_Bin_Index_stretch), RBCStretch(High_Bin_Index_stretch), 'g.')
plot(Time_Axis(Low_Bin_Index_stretch), RBCStretch(Low_Bin_Index_stretch), '.', 'Color', Orange_Color)
hold off
title('5. Bin peaks using stretching')
axis square
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

% Bin peaks using peak finding
Freq_Axis2 = Freq_Axis(((NProj / 2) + 1):end);

[MaxPeak, StrongFreqIndex] = max(FFT_RBC(round(((NProj / 2) + 1):end)));
RBCOscillationFreq = Freq_Axis2(StrongFreqIndex);
text(2, 0.5*MaxPeak, ['Oscillation Frequency: ', num2str(RBCOscillationFreq), ' Hz']);
Heart_Rate = RBCOscillationFreq * 60; % get the heart rate
% call function to bin smoothed data into high and low indices
[High_Bin_Index, Low_Bin_Index] = bin_highlow(RBC_Smoothed_BP, Time_Axis, Heart_Rate, TR, NProj);

subplot(2, 4, 6)
plot(Time_Axis, RBC_Smoothed_BP, 'k.-', 'MarkerSize', 8)
hold on
plot(Time_Axis(High_Bin_Index), RBC_Smoothed_BP(High_Bin_Index), 'g.')
plot(Time_Axis(Low_Bin_Index), RBC_Smoothed_BP(Low_Bin_Index), '.', 'Color', Orange_Color)
hold off
title('6. Bin peaks using peak finding')
axis square
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

subplot(2, 4, 7)
plot(Time_Axis, RBC_Smoothed_BP, 'k.-', 'MarkerSize', 8)
hold on
plot(Time_Axis(High_Bin_Index_stretch), RBC_Smoothed_BP(High_Bin_Index_stretch), 'g.')
plot(Time_Axis(Low_Bin_Index_stretch), RBC_Smoothed_BP(Low_Bin_Index_stretch), '.', 'Color', Orange_Color)
hold off
title('6. Bin peaks comparison peak finding vs threshold')
axis square
xlabel('Time (s)')
ylabel('Intensity (Arb. Units)')

%% Get statistics to put on plot
% Mean oscillation amplitude in the first 6 seconds

Less6ind = find(Time_Axis < 6);
Time6 = Time_Axis(Less6ind);
Data6 = RBC_Smoothed_BP(Less6ind);
[maxima, maxlocs] = findpeaks(Data6, Time6, 'MinPeakDistance', 0.5);
[minima, minlocs] = findpeaks(-Data6, Time6, 'MinPeakDistance', 0.5);
Avg_Amp = (mean(maxima) + mean(minima)) / mean(RBC2Gas) * 100;

Mag_Raw_Sh_Dis = abs(Sh_DisFID(1, :)); %Sh_Mag_k0;
if mean(real(Sh_DisFID(1, :))) > 0
    det_RBC_Raw = real(Sh_DisFID(1, :));
else
    det_RBC_Raw = -real(Sh_DisFID(1, :));
end
if mean(imag(Sh_DisFID(1, :))) > 0
    det_Bar_Raw = imag(Sh_DisFID(1, :));
else
    det_Bar_Raw = -imag(Sh_DisFID(1, :));
end

%Get New RBC to Barrier ratio for high and low keys
mean_Mag_High_k0 = mean(Mag_Raw_Sh_Dis(High_Bin_Index));
mean_Mag_Low_k0 = mean(Mag_Raw_Sh_Dis(Low_Bin_Index));
mean_Mag_Tot_k0 = mean(Mag_Raw_Sh_Dis);

mean_RBC_High_k0 = mean(det_RBC_Raw(High_Bin_Index));
mean_RBC_Low_k0 = mean(det_RBC_Raw(Low_Bin_Index));
mean_RBC_Tot_k0 = mean(det_RBC_Raw);

mean_Bar_High_k0 = mean(det_Bar_Raw(High_Bin_Index));
mean_Bar_Low_k0 = mean(det_Bar_Raw(Low_Bin_Index));
mean_Bar_Tot_k0 = mean(det_Bar_Raw);

High_RBC2Barrier = mean_RBC_High_k0 / mean_Bar_High_k0;
Low_RBC2Barrier = mean_RBC_Low_k0 / mean_Bar_Low_k0;
Tot_RBC2Barrier = mean_RBC_Tot_k0 / mean_Bar_Tot_k0;

RBC2Bar.High = High_RBC2Barrier;
RBC2Bar.Low = Low_RBC2Barrier;
RBC2Bar.Tot = Tot_RBC2Barrier;
RBC2Bar.High2Low = mean_Mag_High_k0 / mean_Mag_Low_k0;
% subplot(2,4,8);
% xlim([0 10]);
% ylim([0 10])
dim = [.75, .35, .1, .1];
str = ['Mean Osc. Amplitude = ', num2str(Avg_Amp), newline, ...
    'RBC/Barrier = ', num2str(RBC2Barrier), newline, ...
    'High Key RBC/Barrier = ', num2str(High_RBC2Barrier), newline, ...
    'Low Key RBC/Barrier = ', num2str(Low_RBC2Barrier), newline, ...
    'Total RBC/Barrier = ', num2str(Tot_RBC2Barrier), newline, ...
    'Mean Tot Dissolved = ', num2str(mean_Mag_Tot_k0), newline, ...
    'Mean High Tot Dissolved = ', num2str(mean_Mag_High_k0), newline, ...
    'Mean Low Tot Dissolved = ', num2str(mean_Mag_Low_k0)];
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'HorizontalAlignment', 'center');

%% Elly is stuck here, why do we detrend the DisFID? - Peter's paper says it is to remove T1
% This is the reason the scaling between recons is different
% is there a way to re-scale but keep the same signal intensity?
Dis2GasFID = DisFID ./ Raw_GasMag_k0;
Sh_DisFID = Sh_DisFID ./ Raw_GasMag_k0;

%% Keyhole Data 1
%In order to decide on Keyhole Radius, calculate sampling percentages
radpts = 1:ceil(size(Traj, 2)/2);
Num_High_Pro = length(High_Bin_Index);
Num_Low_Pro = length(Low_Bin_Index);

Full_Samp_NPro = 4 * pi * radpts.^2;

High_Samp_Pct = Num_High_Pro ./ Full_Samp_NPro * 100;
Low_Samp_Pct = Num_Low_Pro ./ Full_Samp_NPro * 100;

%All the data I've seen so far, we are doubly sampled along the radius, so
%we can go out twice as far as what is calculated:
if ~isfield(params, 'key_radius') %Set Keyhole Radius to the radius at which we are at least 50% sampled
    Key_Radius = max(find((Low_Samp_Pct + High_Samp_Pct)/2 > 50));
    %If there's more than 1 k0 point, we can use a few extra points in the key
    Traj_Radius = squeeze(sqrt(Traj(1, :, :).^2+Traj(2, :, :).^2+Traj(3, :, :).^2));
    Num_k0Pts = length(find(Traj_Radius(:, 1) == 0));
    Key_Radius = Key_Radius * 2 + Num_k0Pts - 1;
else
    Key_Radius = params.key_radius;
end
Keyhole = Dis2GasFID;

if ~isfield(params, 'kspace_filter')
    params.kspace_filter = 'UFC';
end
switch params.kspace_filter
    case 'UFC'
        
        Keyhole(1:Key_Radius, :) = 0;

        High_Key = Keyhole;
        High_Key = High_Key * mean(abs(Dis2GasFID(1, High_Bin_Index))) ./ abs(Dis2GasFID(1, :));
        High_Key(1:Key_Radius, High_Bin_Index) = Dis2GasFID(1:Key_Radius, High_Bin_Index);

        Low_Key = Keyhole;
        Low_Key = Low_Key * mean(abs(Dis2GasFID(1, Low_Bin_Index))) ./ abs(Dis2GasFID(1, :));
        Low_Key(1:Key_Radius, Low_Bin_Index) = Dis2GasFID(1:Key_Radius, Low_Bin_Index);
    case 'UND'
        Keyhole(1:end, :) = 0;

        High_Key = Keyhole;
        High_Key = High_Key * mean(abs(Dis2GasFID(1, High_Bin_Index))) ./ abs(Dis2GasFID(1, :));
        High_Key(1:end, High_Bin_Index) = Dis2GasFID(1:end, High_Bin_Index);

        Low_Key = Keyhole;
        Low_Key = Low_Key * mean(abs(Dis2GasFID(1, Low_Bin_Index))) ./ abs(Dis2GasFID(1, :));
        Low_Key(1:end, Low_Bin_Index) = Dis2GasFID(1:end, Low_Bin_Index);
end

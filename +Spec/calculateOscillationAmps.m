function [amp, detrend, fitted, area_fit, rbcFit, gof] = calculateOscillationAmps(dyn,BHs,b,varargin)

% calculateOscillationAmps returns the oscillation amplitude of all RBC spectral
% parameters (area, freq, fwhm, fwhmG, phase), along with the amplitude of
% the barrier oscillations (b), and the hr. The signal is first run though
% a high pass filter, and then fit to a sine function to determine the
% amplitude.
%
% INPUTS 
% dyn: structure containing dynamic spectroscopy fit information
% BHs: vector with time of BH start and BH end
% b: high pass filter
% 
% Optional: 
%
% getOscillationInfo(dyn,BHs,b,'oscType',method)
% available moethods are:
% 'sine' - sine wave identification of amplitude
% 'peaks' - peak finding identification of amplitude

p = inputParser;
addParameter(p,'oscType','sine',...
   @(x) any(validatestring(x,{'sine','peaks'})));
parse(p,varargin{:});
   
if isfield(dyn, 'fwhmL')
    dyn.fwhm = dyn.fwhmL;
    dyn.fwhmG(1,:) = [];
end 

% Remove first data point to avoid ringing when filtering 
dyn.t(1,:) = []; dyn.area(1,:) = []; dyn.freq(1,:) = []; dyn.fwhm(1,:) = [];
dyn.phase(1,:) = []; %dyn.snrs(1) = [];

if BHs(2) > dyn.t(end)
    BHs(2) = dyn.t(end);
end 

[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);


%% Oscillation Information

clear opts
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [0 0 -Inf];
opts.Upper = [Inf Inf Inf];
xData = dyn.t(BHstart:BHend)';

% [rbcFit, gof] = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,1),'exp2'); % biexponential decay
[rbcFit, gof] = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,1),'exp1');
rbcNorm = rbcFit(dyn.t(:,1));

area_detrend = filtfilt(b, 1, (dyn.area(:,1)-rbcNorm)./rbcNorm);
area_detrend = area_detrend(BHstart:BHend,1);

[area_fit, area_gof] = fit( xData, area_detrend, ft, opts );
areaVal = area_fit(dyn.t(BHstart:BHend,1));
amp_area = abs(area_fit.a1*2);
hr = area_fit.b1/(2*pi)*60;

freq_detrend = filtfilt(b,1,dyn.freq(:,1)-mean(dyn.freq(BHstart:BHend,1)));
freq_detrend = freq_detrend(BHstart:BHend,1);
    
fwhm_detrend = filtfilt(b,1,dyn.fwhm(:,1)-mean(dyn.fwhm(BHstart:BHend,1)));
fwhm_detrend = fwhm_detrend(BHstart:BHend,1);

dyn.phase(:,1) = rad2deg(unwrap(deg2rad(dyn.phase(:,1))));
phase_detrend = filtfilt(b,1,dyn.phase(:,1)-mean(dyn.phase(BHstart:BHend,1)));
phase_detrend = phase_detrend(BHstart:BHend,1);

bFit = fit(dyn.t(BHstart:BHend,1),dyn.area(BHstart:BHend,2),'exp1');
bNorm = bFit(dyn.t(:,1));
b_detrend = filtfilt(b,1,(dyn.area(:,2)-bNorm)./bNorm);
b_detrend = b_detrend(BHstart:BHend);
b_fit = fit( xData, b_detrend, ft, opts );
b_area = abs(b_fit.a1*2);

switch p.Results.oscType   
case 'sine' 
% use sine wave fitting to determine oscillation amplitudes

    % Set frequency to hr from amplitude
    opts.Lower = [0 area_fit.b1 -Inf];
    opts.Upper = [Inf area_fit.b1 Inf];

    [freq_fit, freq_gof] = fit( xData, freq_detrend, ft, opts );
    freqVal = freq_fit(dyn.t(BHstart:BHend,1));
    amp_freq = freq_fit.a1*2;

    fwhm_fit = fit( xData, fwhm_detrend, ft, opts );
    fwhmVal = fwhm_fit(dyn.t(BHstart:BHend,1));
    amp_fwhm = fwhm_fit.a1*2;

    phase_fit = fit( xData, phase_detrend, ft, opts );
    phaseVal = phase_fit(dyn.t(BHstart:BHend,1));
    amp_phase = phase_fit.a1*2;
    
    % Create Structure for Easy Saving of Data
    fnamesF = {'area','freq','fwhm','phase'};

    fitted.(fnamesF{1}) = areaVal;
    fitted.(fnamesF{2}) = freqVal;
    fitted.(fnamesF{3}) = fwhmVal;
    fitted.(fnamesF{4}) = phaseVal;
    
case 'peaks'
    % use peak identification to determine oscillation amplitudes
    
    [~, ay] = findpeaks(smooth(area_detrend,10).^2,'MinPeakDistance',hr/60/6*50);
    area_max = area_detrend(ay); ay = dyn.t(ay+BHstart-1,1);
    amp_area = median(area_max(area_max>0))-median(area_max(area_max<0));
    
    [~, fy] = findpeaks(smooth(freq_detrend,10).^2,'MinPeakDistance',hr/60/6*50);
    freq_max = freq_detrend(fy); fy = dyn.t(fy+BHstart-1,1);
    amp_freq = median(freq_max(freq_max>0))-median(freq_max(freq_max<0));
    
    [~, fwy] = findpeaks(smooth(fwhm_detrend,10).^2,'MinPeakDistance',hr/60/6*50);
    fwhm_max = fwhm_detrend(fwy); fwy = dyn.t(fwy+BHstart-1,1);
    amp_fwhm = median(fwhm_max(fwhm_max>0))-median(fwhm_max(fwhm_max<0));
    
    [~, py] = findpeaks(smooth(phase_detrend,10).^2,'MinPeakDistance',hr/60/6*50);
    phase_max = phase_detrend(py); py = dyn.t(py+BHstart-1,1);
    amp_phase = median(phase_max(phase_max>0))-median(phase_max(phase_max<0));
    
    % Calculate peak-to-peak based on sine RMSE
    % amp_area = std(area_detrend)*2*sqrt(2); 
    % amp_freq = std(freq_detrend)*2*sqrt(2);
    % amp_fwhm = std(fwhm_detrend)*2*sqrt(2);
    % amp_phase = std(phase_detrend)*2*sqrt(2);
   
    % Plot peak identification results
    figure(10), clf 
    subplot(4,1,1), hold on
    plot(dyn.t(BHstart:BHend,1),area_detrend*100,'Linewidth',3)
    plot(ay,area_max*100,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim(BHs), ylabel('Amplitude (%)')
    subplot(4,1,2), hold on
    plot(dyn.t(BHstart:BHend,1),freq_detrend,'Linewidth',3)
    plot(dyn.t(BHstart:BHend,1),smooth(freq_detrend,20),'Linewidth',2)
    plot(fy,freq_max,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim(BHs), ylabel('Freq (ppm)')
    subplot(4,1,3), hold on
    plot(dyn.t(BHstart:BHend,1),fwhm_detrend,'Linewidth',3)
    plot(dyn.t(BHstart:BHend,1),smooth(fwhm_detrend,20),'Linewidth',2)
    plot(fwy,fwhm_max,'.k','MarkerSize',30);
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim(BHs), ylabel('FWHM (ppm)')
    subplot(4,1,4), hold on
    plot(dyn.t(BHstart:BHend,1),phase_detrend,'Linewidth',3)
    plot(dyn.t(BHstart:BHend,1),smooth(phase_detrend,20),'Linewidth',2)
    plot(py,phase_max,'.k','MarkerSize',30); 
    plot(dyn.t(BHstart:BHend,1),zeros(1,length(BHstart:BHend)),'Linewidth',.25,'Color',[.5 .5 .5])
    xlim(BHs), xlabel('Time (s)'), ylabel('Phase ({\circ})')
    
    % Create Structure for Easy Saving of Data
    fnamesF = {'area','at','freq','ft','fwhm','fwt','phase','pt'};

    fitted.(fnamesF{1}) = area_max; fitted.(fnamesF{2}) = ay;
    fitted.(fnamesF{3}) = freq_max; fitted.(fnamesF{4}) = fy;
    fitted.(fnamesF{5}) = fwhm_max; fitted.(fnamesF{6}) = fwy;
    fitted.(fnamesF{7}) = phase_max; fitted.(fnamesF{8}) = py;
end 
    
% Create Structure for Easy Saving of Data
switch p.Results.oscType   
    case 'sine'     
        fnames = {'area','freq','fwhm','phase','hr','b','area_gof', 'freq_gof'}; 
        amp.(fnames{1}) = amp_area;
        amp.(fnames{2}) = amp_freq;
        amp.(fnames{3}) = amp_fwhm;
        amp.(fnames{4}) = amp_phase;
        amp.(fnames{5}) = hr;
        amp.(fnames{6}) = b_area;
        amp.(fnames{7}) = area_gof;
        amp.(fnames{8}) = freq_gof;
    case 'peaks'
        fnames = {'area','freq','fwhm','phase','hr','b','area_gof'};
        amp.(fnames{1}) = amp_area;
        amp.(fnames{2}) = amp_freq;
        amp.(fnames{3}) = amp_fwhm;
        amp.(fnames{4}) = amp_phase;
        amp.(fnames{5}) = hr;
        amp.(fnames{6}) = b_area;
        amp.(fnames{7}) = area_gof;
end 


% Create Structure for Easy Saving of Data
fnamesD = {'area','freq','fwhm','phase'};

detrend.(fnamesD{1}) = area_detrend;
detrend.(fnamesD{2}) = freq_detrend;
detrend.(fnamesD{3}) = fwhm_detrend;
detrend.(fnamesD{4}) = phase_detrend;

end 
 
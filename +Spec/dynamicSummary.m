function [amp, nmrFit_ppm] = dynamicSummary(raw_path,dyn_path,BHs,save_fig_flag,figname,varargin)
% raw_path: raw data file either .dat or .7
% dyn_path: dynamic spectroscopy matlab structure or structure location
% BHs: either a vector containing the time of inhale/exhale or the filepath
%      for a saved matlab vector with the BH information. Default
%      values are [2 10].
% save_fig_flag: logical input for saving figures. Figures are saved in the
%                same folder as dyn_path.
% figname: adds text after default figure name (e.g. '_2' -> dynV_2)

% Optional: 
%
% dynamicSummary(raw_path,dyn_path,BHs,save_fig_flag,figname,'oscType',method)
% available methods are:
% 'sine' - sine wave identification of amplitude
% 'peaks' - peak finding identification of amplitude

% Example inputs
% raw_path = 'D:\Elly\Documents\Duke\CIVM Research\Siemens Data\1 - Healthy Volunteers\18-10-01 SUBJECT 000-001C\2 - Raw\meas_MID00026_FID69955_Xe_fid_DynamicSpec_high_flip.dat';
% 
% dyn_name = 'dynV5_2';
% 
% dyn_path = [path(1:end-7),'4 - Spectra\',dyn_name,'.mat'];
% 
% BHs = [3 10];
% save_fig_flag = 0;
% 
% figname = '';
% fig_path = fileparts(dyn_path);
%%

p = inputParser;
addParameter(p,'suppressFigs',false, @islogical);
addParameter(p,'oscType','sine',...
   @(x) any(validatestring(x,{'sine','peaks'})));
parse(p,varargin{:});


%%
if isempty(raw_path)
    disp('Please select the raw dat file')
    raw_path = filepath();
    dyn_folder = fileparts(raw_path);
else 
    dyn_folder = fileparts(raw_path);
end 

if isempty(dyn_path)
    disp('Please select the dyn variable')
    dyn_path = filepath();
    load(dyn_path)
elseif ischar(dyn_path)
    load(dyn_path)
    [dyn_folder] = fileparts(dyn_path);
else
    dyn = dyn_path;
end

if isempty(BHs)
    BHs = [2 10];
elseif ischar(BHs)
    load(BHs);
elseif isvector(BHs)
end 
   
[BHstart, BHend] = findBHs(dyn.t(:,1), BHs);

if ~isfield(dyn,'fwhmG') % if the data is not Voigt, set all FWHM_G to 0
    dyn.fwhmG = zeros(size(dyn.t));
    dyn.fwhmL = dyn.fwhm;
    fitType = 'L';
else
    fitType = 'V';
end


%% Calculate Oscillations Amplitudes and Static Parameters

b = highpassfilter(length(dyn.area(:,1)));
[amp, detrend, fitted] = calculateOscillationAmps(dyn,BHs,b,'oscType',p.Results.oscType);
amp.snr = mean(dyn.snrs(BHstart:BHend));

% Display oscillation results
disp([10 'file: ',raw_path])
disp([10 'RBC Dynamic Oscillation Amplitudes:'])
disp('   Area      Freq      FWHM      Phase      HR      SNR');
disp([sprintf('%8.1f', amp.area*100), ' ' ...
    sprintf('%8.2f',amp.freq),  '  ' ...
    sprintf('%8.2f',amp.fwhm), '  ' ...
    sprintf('%8.1f',amp.phase), '   '...
    sprintf('%6.0f',amp.hr),'  ',...
    sprintf('%7.1f',amp.snr)]);

%% Calculate Oscillations Amplitudes and Static Parameters
[nmrFit, nmrFit_ppm, ~] = calculateStaticSpectroscopy(raw_path, BHs, fitType);

%% Plots

if p.Results.suppressFigs == false
    
    dynSave_path = [dyn_folder,'\dyn',fitType,figname,'.tif'];
    plotDynamics(dyn, BHs, save_fig_flag, dynSave_path);

    oscSave_path = [dyn_folder,'\dyn',fitType,'_oscs_new',figname,'.tif'];
    plotOscillations(dyn, BHs, detrend, fitted, save_fig_flag, oscSave_path)

    snrSave_path = [dyn_folder,'\dyn',fitType,'_error',figname,'.tif'];
    plotDynSNR(dyn, BHs, save_fig_flag, snrSave_path)
    
    staticSave_path = [dyn_folder,'\static',fitType,figname,'.tif'];
    plotStaticSpectroscopy(raw_path, nmrFit, save_fig_flag, staticSave_path)
end 

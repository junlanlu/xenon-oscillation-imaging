function [fids, dwell_time, npts, tr, xeFreqMHz] = readRawDyn(raw_path)

[~, ~, scanner] = fileparts(raw_path);

% Read in twix or P file and define associated variables
if strcmp(scanner,'.dat')
    % Twix file from Siemens
    twix = readtwix(raw_path);
    npts = size(twix.data,1);                          % Number of samples per FID                  
    dwell_time = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*10^-9;  % Receiver bandwidth (kHz); works for both calibration types                         
    tr = twix.hdr.Config.TR(1)*1E-6;                       % Time between each sample
    fids = twix.data;
    if isfield(twix.hdr.Config,'Frequency')
        % UVA Siemens File 
        xeFreqMHz = twix.hdr.Config.Frequency*10e-7; %34.093484
    elseif isfield(twix.hdr.Meas,'lFrequency')
        % Duke Siemens File
        xeFreqMHz = twix.hdr.Meas.lFrequency*10e-7; %34.091516
    end 


    % Magnetic Field Strength
    mag_fstregth = twix.hdr.Dicom.flMagneticFieldStrength;
    
    % This if, else block will read the exitation from different fields,
    % depending on our scan version
    
    if isfield(twix.hdr.Phoenix, 'sWiPMemBlock')
        if isfield(twix.hdr.Phoenix.sWiPMemBlock,'adFree')
            excitation = twix.hdr.Phoenix.sWiPMemBlock.adFree{1, 9};
        else
            excitation = 7091; % assumes 218ppm
        end
    elseif isfield(twix.hdr.Phoenix, 'sWipMemBlock')
        if isfield(twix.hdr.Phoenix.sWipMemBlock,'alFree')
            excitation = twix.hdr.Phoenix.sWipMemBlock.alFree{1, 5};
        end
    else
        excitation = 7091; % assumes 218ppm
    end
    
    gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
    
    % RF excitation will be in ppm, either 218ppm or 208 ppm
    rf_excitation = round(excitation/(gyro_ratio * mag_fstregth));

elseif strcmp(scanner,'.7')
    % Pfile from GE
    pfile = GE.Pfile.read(raw_path);
    MRI.DataProcessing.checkForOverranging(pfile); % Check for overranging
    pfile = MRI.DataProcessing.removeBaselineViews(pfile); % Remove baselines
    bw = 1000*pfile.rdb.rdb_hdr_user12;                    % Receiver bandwidth (kHz)
    dwell_time = 1/(2*bw);                                 % Time between each sample
    dwell_time = Math.nearestMultipleOf(dwell_time,0.000002);
    npts = pfile.rdb.rdb_hdr_frame_size;
    tr = pfile.image.tr*1E-6;
    fids = pfile.data;
    xeFreqMHz = 17.660445;
else 
    error('Unknown Raw File Type')
end 

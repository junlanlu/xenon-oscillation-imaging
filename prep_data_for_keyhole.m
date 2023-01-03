function [dissolved_data, gas_data, TR] = prep_data_for_keyhole(varargin)

if(nargin < 1 | ~exist(varargin{1}))
    disp('Select Dixon File');
    dixon_file = filepath();
else
    dixon_file = varargin{1};
end

% For most recent format (65)
output_image_sizeg = 64*[1 1 1];
output_image_sized = 64*[1 1 1];

[~, ~, scanner] = fileparts(dixon_file);

% Read in twix or P file and define associated variables
if strcmp(scanner,'.7')
    verbose = 0; 
    
    % Prepare to override values
    pfileOverride = GE.Pfile.Pfile();
    pfileOverride.rdb.rdb_hdr_user1  = 0.512; % pw_gxwa
    pfileOverride.rdb.rdb_hdr_user38 = 0.2;  % pw_gxwd/1000
    pfileOverride.rdb.rdb_hdr_user44 = 1.536; % pw_gxw/1000
    pfileOverride.rdb.rdb_hdr_user22 = 0.125; %toff
    pfileOverride.rdb.rdb_hdr_user32 = 0;  % Golden Means
    rmBline = 0;
    rmFirstGas = 1;
    rmFirstDis = 0;

    deltaf_gas = 0;
    deltaf_dissolved = 0; %hz

    % Gradient delays
    delays.x_delay = 0.000;
    delays.y_delay = 0.00;
    delays.z_delay = 0.000;

    %% Read Raw Pfile and process pfile
    pfile = GE.Pfile.read(dixon_file);

    % Convert from Pfile format
    pfile = convertLegacyPfile(pfile);

    % Override header values (optional)
    displayPfileHeaderInfo(pfile);
    pfile = overridePfile(pfile, pfileOverride);
    if(verbose)
        % Display key header info
        displayPfileHeaderInfo(pfile);
    end

    % Check for overranging
    MRI.DataProcessing.checkForOverranging(pfile);

    % Remove baselines
    if(rmBline)
        pfile = MRI.DataProcessing.removeBaselineViews(pfile);
    end

    %% Split pfile into 2: disolved and gas
    if(~exist('startDissolved') & ~exist('startGas'))
        if(pfile.rdb.rdb_hdr_ps_mps_freq == 176604450)
            % dissolved first
            startDissolved = 2;
            startGas = 1;
        else
            % gas first
            startDissolved = 1;
            startGas = 2;
        end
    end
    dissolved_pfile = pfile;
    dissolved_pfile.data = dissolved_pfile.data(:,startDissolved:2:end-1);
    if(rmFirstDis)
        dissolved_pfile.data = dissolved_pfile.data(:,2:end);
    end
    dissolved_pfile.rdb.rdb_hdr_user20 = size(dissolved_pfile.data,2);

    gas_pfile = pfile;
    gas_pfile.data = gas_pfile.data(:,startGas:2:end);
    if(rmFirstGas)
        gas_pfile.data = gas_pfile.data(:,2:end);
    end
    gas_pfile.rdb.rdb_hdr_user20 = size(gas_pfile.data,2);

    %% Calculate trajectories (will be same for both pfiles)
    % Calculate trajectory for a single radial ray
    radialDistanceg = MRI.Trajectories.Centric.Radial.calcRadialRay(gas_pfile, delays, output_image_sizeg);
    radialDistanced = MRI.Trajectories.Centric.Radial.calcRadialRay(dissolved_pfile, delays, output_image_sized);

    % Distribute rays onto 3d sphere
    trajd = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistanced, dissolved_pfile);
    trajg = MRI.Trajectories.Centric.Distribute.calculate3dTrajectories(radialDistanceg, gas_pfile);

    % Undo loopfactor
    [trajd, dissolved_pfile] = MRI.DataProcessing.undoloopfactor(trajd, dissolved_pfile);
    [trajg, gas_pfile] = MRI.DataProcessing.undoloopfactor(trajg, gas_pfile);

    % Demodulate signal
    bw = pfile.rdb.rdb_hdr_user12;                 % Receiver bandwidth (kHz)
    dwell_time = 1/(2*bw*1000);                                             % Time between each sample
    dwell_time = Math.nearestMultipleOf(dwell_time,0.000002); % Dwell time must be an integer multible of 2us
    nPts = pfile.rdb.rdb_hdr_frame_size;           % Number of sample points per frame/ray
    tMatd = repmat(dwell_time*(1:nPts)',[1 size(dissolved_pfile.data,2)]);
    tMatg = repmat(dwell_time*(1:nPts)',[1 size(gas_pfile.data,2)]);
    dissolved_pfile.data = dissolved_pfile.data.*exp(1i*2*pi*deltaf_dissolved*tMatd);
    gas_pfile.data = gas_pfile.data.*exp(1i*2*pi*deltaf_gas*tMatg);
    
    TR = 0.015;

    % [PAUSE here] this part of code is used for cincinnati data preparation
    dissolved_data = dissolved_pfile.data;
    gas_data = gas_pfile.data;
    traj = permute(trajd, [3 1 2]);
elseif strcmp(scanner,'.dat')
    % Twix file from Siemens
    twix_obj = mapVBVD(dixon_file);   
    raw_fids = squeeze(double(twix_obj.image.unsorted()));
    raw_fids(:,1) = [];
        
    if length(raw_fids) == 2002
        raw_fids(:,end-1:end) = [];
%         load('traj_siemens.mat');
    elseif length(raw_fids) == 2031
        raw_fids(:,end-31:end) = [];
%         load('traj_siemens2.mat');
    else 
        raw_fids(:,end-1:end) = [];
%         raw_fids(:,end-31:end) = [];
%         load('traj_siemens3.mat');
    end 
        
%     traj = [];
    
    dissolved_data = conj(raw_fids(:,3:2:end));
    gas_data = conj(raw_fids(:,2:2:end));
       
    for ii = 25:size(dissolved_data,1)
        for jj = 1:size(dissolved_data,2)
            if abs(dissolved_data(ii,jj)) > mean(abs(dissolved_data(1,:)))
                dissolved_data(ii,jj) = 0;
            end
        end
    end
        
    TR = 0.015;
else    
    error('Unknown Raw File Type')
end 
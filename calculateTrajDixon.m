function [data_dis, data_gas, traj_dis, traj_gas, data_type] = calculateTrajDixon(twix_obj, data, nSpec)

flag = 3; 
% determine number of points/ray, usually 64 or 128
npts = twix_obj.image.NCol;
if twix_obj.image.flagRemoveOS
    npts = npts/2;
end
nFrames = twix_obj.image.NAcq;
% calculate dwell time
if isfield(twix_obj.hdr.Meas,'alDwellTime')
    dwell_time = twix_obj.hdr.Meas.alDwellTime(1)*1E-3; % in usec maybe 2?
else
    dwell_time = cell2mat(twix_obj.hdr.Phoenix.sRXSPEC.alDwellTime(1)); % array
    dwell_time = dwell_time(1);
end

if twix_obj.image.flagRemoveOS
    dwell_time = dwell_time*2;
end

% determine oversampling
oversampling = 3; %maybe 2
% determine decay_time
decay_time = 60; % hard-coded

% determine ramp time and plat_time
if isfield(twix_obj.hdr.Meas,'RORampTime')
    ramp_time = twix_obj.hdr.Meas.RORampTime; % read ramp time from header

    % GD_userdef = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{4};
elseif isfield(twix_obj.hdr.Meas,'alRegridRampupTime')
    ramp_time = twix_obj.hdr.Meas.alRegridRampupTime(1);
    % GD_userdef = [];
else
    ramp_time = 100; % us
end 
% set minimum to 100us
if ramp_time == 0
    ramp_time = 100;
end
plat_time = 2500;
[yyyy, mm, dd] = Tools.getScanDate(twix_obj);
%% determine data points in data variable that are dixon
protocol_name = twix_obj.hdr.Config.ProtocolName;
[row, col] = size(data);
if contains(protocol_name, 'fast')
    disp('fast dixon');
    data_dixon = data(:,1:end-30);
    data_gas = data_dixon(:,1:2:end);
    data_dis_raw = data_dixon(:,2:2:end);
    nSkip_dis = 0;
    data_dis = data_dis_raw(:,nSkip_dis+1:end);
    nFrames = 2*size(data_gas, 2);
    del_x = -5; del_y = -5; del_z = -5;
    nSkip_start = 0;
    nSkip_end = 0;
    data_type = 'fast_dixon';
elseif contains(protocol_name, 'med')
    disp('medium dixon');
    data_dixon = data(:,1:end-30);
    data_gas = data_dixon(:,1:2:end);
    data_dis_raw = data_dixon(:,2:2:end);
    nSkip_dis = 0;
    data_dis = data_dis_raw(:,nSkip_dis+1:end);
    nFrames = 2*size(data_gas, 2);
    del_x = -5; del_y = -5; del_z = -5;
    nSkip_start = 0;
    nSkip_end = 0;
    data_type = 'medium_dixon';
elseif contains(protocol_name, '2008') || contains(protocol_name,'2007')
    disp('2007 or 2008 dixon sequence, intermediate transition');
    data_dixon = data(:,1:end)*exp(1j*pi/2);
    data_gas = data_dixon(:,1:2:end);
    data_dis_raw = data_dixon(:,2:2:end);
    nSkip_dis = 0;
    data_dis = data_dis_raw(:,nSkip_dis+1:end);
    nFrames = 2*size(data_gas, 2);
    del_x = 0; del_y = -4; del_z = -3;
    nSkip_start = 0;
    nSkip_end = 0;
    thre_snr = 1.5;
    data_type = 'normal_2007-2008';
elseif col == 2002 && datetime(str2double(yyyy), str2double(mm), str2double(dd)) > datetime(2018, 05, 01)
    disp('normal dixon: pre 2021 scanner upgrade highBW');
    nFrames = twix_obj.image.NAcq;
    nSkip_dis = 0;
    data_dixon = data(:,1:end-2);
    del_x = 0; del_y = -4; del_z = -3; % gradient delay
    data_gas = data_dixon(:,3:2:end);
    data_dis = data_dixon(:,4:2:end);
    nSkip_start = 1;
    nSkip_end = 1;
    data_type = 'normal_trio';
elseif col == 2002 && datetime(str2double(yyyy), str2double(mm), str2double(dd)) <= datetime(2018, 05, 01)
    disp('normal dixon: pre 2021 scanner upgrade lowBW');
    nFrames = twix_obj.image.NAcq;
    nSkip_dis = 0;
    data_dixon = data(:,1:end-2);
    del_x = 24; del_y = 22; del_z = 22; % gradient delay
    data_gas = data_dixon(:,3:2:end);
    data_dis = data_dixon(:,4:2:end);
    nSkip_start = 1;
    nSkip_end = 1;
    data_type = 'normal_trio';
elseif col == 2032
    disp('normal dixon: pre 2021 scanner upgrade w/ bonus spectra');
    nSkip_dis = 0;
    data_dixon = data(:,1:end-32);
    del_x = 0; del_y = -4; del_z = -3; % gradient delay
    data_gas = data_dixon(:,3:2:end);
    data_dis_raw = data_dixon(:,4:2:end);
    data_dis = data_dis_raw(:,nSkip_dis+1:end);
    nFrames = 2032;
    nSkip_start = 1;
    nSkip_end = 16;
    data_type = 'normal_trio_bonus';
elseif col == 2030
    disp('normal dixon: post 2021 scanner upgrade');
    data_dixon = data(:,1:end-30);
    nSkip_dis = 0;
    data_gas = data_dixon(:,1:2:end);
    data_dis_raw = data_dixon(:,2:2:end);
    data_dis = data_dis_raw(:,nSkip_dis+1:end);
    nFrames = 2000;
    del_x = -5; del_y = -5; del_z = -5;
    nSkip_start = 0;
    nSkip_end = 0;
    data_type = 'normal_prisma';
else
    disp('invalid sequence');
    ME = MException('myComponent:inputError','Sequence not supported');
    throw(ME);
end


%% set up trajectory


radialDistance_x = generate_radial_1D_traj(dwell_time, del_x, ramp_time, plat_time, decay_time, npts, oversampling);
radialDistance_y = generate_radial_1D_traj(dwell_time, del_y, ramp_time, plat_time, decay_time, npts, oversampling);
radialDistance_z = generate_radial_1D_traj(dwell_time, del_z, ramp_time, plat_time, decay_time, npts, oversampling);

if twix_obj.image.flagRemoveOS
    radialDistance_x = radialDistance_x*2;
    radialDistance_y = radialDistance_y*2;
    radialDistance_z = radialDistance_z*2;
end

m_lNumberOfFrames = 1;

%% generating trajectory
if(flag == 1 || flag == 3) %Spiral or haltonized spiral
    [x,y,z] = GX_f_gen_traj(floor(nFrames/2), m_lNumberOfFrames,flag); % 1 for spiral , others for Halton
         
else %Halton
    [x,y,z] = GX_f_gen_traj(round(nFrames/2), m_lNumberOfFrames,flag); % 1 for spiral , others for Halton
end

%% calculate trajectory
x = radialDistance_x(:,1)*x;
y = radialDistance_y(:,2)*y;
z = radialDistance_z(:,3)*z;

%dissolved traj, throw away the first one as the gas is contaminated
x_dis = x(:,1+nSkip_start+nSkip_dis:end-nSkip_end);
y_dis = y(:,1+nSkip_start+nSkip_dis:end-nSkip_end);
z_dis = z(:,1+nSkip_start+nSkip_dis:end-nSkip_end);
%data_dis = data(:,4:2:end-1); % dissolved

%gas traj
x_gas = x(:,1+nSkip_start:end-nSkip_end); 
y_gas = y(:,1+nSkip_start:end-nSkip_end);
z_gas = z(:,1+nSkip_start:end-nSkip_end);
%data_gas = data(:,3:2:end-2); % gas

%% Throw away the rays with spiky noise
nFrames_dis_raw = size(data_dis,2);

if ~ exist('thre_snr')
    thre_snr = 0.7; % default value of 0.7
end
good_indices_dis = removeNoiseRays( data_dis,x_dis,y_dis,z_dis, thre_snr); 
good_indices_gas = removeNoiseRays( data_gas,x_gas,y_gas,z_gas, thre_snr);

good_indices = intersect(good_indices_dis, good_indices_gas);
fprintf(2,['Dissolved-phase and gas data cleanse: threw away noisy rays ',num2str(nFrames_dis_raw-length(good_indices)),'\n']);
nFrames_gas = length(good_indices);
nFrames_dis = nFrames_gas;
x_gas = x_gas(:,good_indices); 
y_gas = y_gas(:,good_indices); 
z_gas = z_gas(:,good_indices);  
data_gas = data_gas(:, good_indices);

x_dis = x_dis(:,good_indices); 
y_dis = y_dis(:,good_indices); 
z_dis = z_dis(:,good_indices);  
data_dis = data_dis(:, good_indices);
assert(nFrames_gas == nFrames_dis);

%% Vectorize data and trajectory
traj_scale = npts/128;

data_gas= reshape(data_gas, [npts*nFrames_gas 1]);
traj_gas = reshape(x_gas, [npts*nFrames_gas 1]);
traj_gas(:,2) = reshape(y_gas, [npts*nFrames_gas 1]);
traj_gas(:,3) = reshape(z_gas, [npts*nFrames_gas 1]);
traj_gas = 0.5*traj_scale*traj_gas;

data_dis = reshape(data_dis, [npts*nFrames_dis 1]);
traj_dis = reshape(x_dis, [npts*nFrames_dis 1]);
traj_dis(:,2) = reshape(y_dis, [npts*nFrames_dis 1]);
traj_dis(:,3) = reshape(z_dis, [npts*nFrames_dis 1]);
traj_dis = 0.5*traj_scale*traj_dis;




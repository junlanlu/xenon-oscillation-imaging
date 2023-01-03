%% Siemens
% Mean oscillation
healthy_thresh_Osc = [-3.6436, .2917, 9.3565, 17.9514, 28.7061, 42.6798, 61.8248];
% old values
% healthy_thresh_Osc = [-6.3404    1.9544   10.2492   18.5440   26.8388   35.1337   43.4285];
healthy_lambda_Osc_boxcox = -1.7541;
healthy_constant_Osc_boxcox = 100;
healthy_mean_Osc_boxcox = 0.57;
healthy_std_Osc_boxcox = 1.8798e-05;

healthy_mean_Osc = 10.25; % percent
healthy_std_Osc = 8.30; % percent
healthy_high_Osc = 2.83; % percent
healthy_defect_Osc = 2.31; % percent
healthy_defectlow_Osc = 15.13; % percent
% Mean RBC
healthy_mean_RBCBar = 0.5463;
healthy_std_RBCBar = 0.2507;


% Define the dict for 8 oscillation bins
% index2color_eightbin = dict;
% index2color_eightbin(0) = [0, 0, 0];
% index2color_eightbin(1) = [1, 0, 0];
% index2color_eightbin(2) = [1, 0.7143, 0];
% index2color_eightbin(3) = [0.4, 0.7, 0.4];
% index2color_eightbin(4) = [0, 1, 0];
% index2color_eightbin(5) = [184.0/255.0, 226.0/255.0, 145.0/255.0];
% index2color_eightbin(6) = [243.0/255.0, 205.0/255.0, 213.0/255.0];
% index2color_eightbin(7) = [225.0/255.0, 129.0/255.0, 162.0/255.0];
% index2color_eightbin(8) = [197.0/255.0, 27.0/255.0, 125.0/255.0];
% index2color_eightbin(9) = [0.33, 0.33, 0.33];
keyset = 0:9;
valueset = {...
        [0, 0, 0],...
        [1, 0, 0],...
        [1, 0.7143, 0],...
        [0.4, 0.7, 0.4],...
        [0, 1, 0],...
        [184.0/255.0, 226.0/255.0, 145.0/255.0],...
        [243.0/255.0, 205.0/255.0, 213.0/255.0],...
        [225.0/255.0, 129.0/255.0, 162.0/255.0],...
        [197.0/255.0, 27.0/255.0, 125.0/255.0],...
        [0.33, 0.33, 0.33]...
    };

% index2color_eightbin_map = containers.Map(keyset, valueset);
index2color_eightbin_map = valueset;


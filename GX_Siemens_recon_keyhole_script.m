%% Get the config file
clear all; close all;
if ~exist('pathbundle', 'var') 
    disp('Locate the config file');
    config_folder = '+Filelists\option\';
    if exist(config_folder, 'dir')
        pathbundle.config_file = filepath(config_folder);
    else
        pathbundle.config_file = filepath();
    end
end  
%% load the params from pathbundle
[filepath, name, ext] = fileparts(pathbundle.config_file);
c = name;
fh = str2func(c);
params = fh();
%% Run the pipeline
GX_Siemens_recon_keyhole(params);

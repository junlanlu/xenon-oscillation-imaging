function params = example()
    params = {};
    params.ScannerType = 'Siemens';
    params.subject_name = '007-02-001A';
    params.dixon_file = 'D:\Patients\007-02-001A\meas_MID00459_FID03003_xe_radial_Dixon_cor_fast_2105__67.dat';
    params.mask_file = 'D:\Patients\007-02-001A\007-02-001A_gx.mat';
    params.spect_file = 'D:\Patients\007-02-001A\meas_MID00458_FID03002_fid_xe_calibration_2108__67.dat';
    [params.filepath, params.filename] = fileparts(params.dixon_file);
    params.N_Begin_Cutoff = 100;
    params.N_End_Cutoff = 0;
    params.key_radius = 9;
end
function [traj_dis, traj_gas] = generate_traj(trajMap)
    %% get trajectory parameters from trajMap
    npts = trajMap('npts');
    nFrames = trajMap('nFrames');
    traj_type = trajMap('traj_type');
    dwell_time = trajMap('dwell_time');
    oversampling = trajMap('oversampling');
    ramp_time = trajMap('ramp_time');
    plat_time = trajMap('plat_time');
    decay_time = trajMap('decay_time');
    del_x = trajMap('del_x');
    del_y = trajMap('del_y');
    del_z = trajMap('del_z');


    %% generating trajectory
    radialDistance_x = generate_radial_1D_traj(dwell_time, del_x, ramp_time, plat_time, decay_time, npts, oversampling);
    radialDistance_y = generate_radial_1D_traj(dwell_time, del_y, ramp_time, plat_time, decay_time, npts, oversampling);
    radialDistance_z = generate_radial_1D_traj(dwell_time, del_z, ramp_time, plat_time, decay_time, npts, oversampling);

    m_lNumberOfFrames = 1;
    [x,y,z] = GX_f_gen_traj(floor(nFrames/2), m_lNumberOfFrames, traj_type);
    
end
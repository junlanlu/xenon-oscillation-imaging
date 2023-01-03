function recon = reconDixonLowRes(ImageSize,data,traj)

kernel.sharpness = 0.14;
kernel.extent = 9*kernel.sharpness;
overgrid_factor = 3;
output_image_size = ImageSize*[1 1 1];
nDcfIter = 15;
deapodizeImage = false;%true();
cropOvergriddedImage = true();
verbose = false();

%  Choose kernel, proximity object, and then create system model
kernelObj = Recon.SysModel.Kernel.Gaussian(kernel.sharpness, kernel.extent, verbose);
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);

clear kernelObj;
systemObj = Recon.SysModel.MatrixSystemModel(traj, overgrid_factor, ...
    output_image_size, proxObj, verbose);

if ImageSize <= 64
    high_overgrid_factor = 4;
    systemObj_dcf = Recon.SysModel.MatrixSystemModel(traj, high_overgrid_factor, ...
    output_image_size, proxObj, verbose);
    % Choose density compensation function (DCF)
    dcfObj = Recon.DCF.Iterative(systemObj_dcf, nDcfIter, verbose);
else 
    dcfObj = Recon.DCF.Iterative(systemObj, nDcfIter, verbose);
end 

% Choose Reconstruction Model
reconObj = Recon.ReconModel.LSQGridded(systemObj, dcfObj, verbose);
clear modelObj;
clear dcfObj;
reconObj.crop = cropOvergriddedImage;
reconObj.deapodize = deapodizeImage;


% Reconstruct image using trajectories in pixel units
recon = reconObj.reconstruct(data, traj);

# xenon-oscillation-imaging

My version of oscillation imaging of 129Xe originally developed by Niedbalski et al.

## Setup

1. You need to install the [non-Cartesian reconstruction package](https://github.com/ScottHaileRobertson/Non-Cartesian-Reconstruction) developed by Scott Robertson
2. Install the [Duke CIVM toolkit](https://github.com/ScottHaileRobertson/Duke-CIVM-MRI-Tools) developed by Scott Robertson
3. Possible add the +Spec folder containing spectroscopy processing functions to the Matlab Path



## Usage

1. Create a config file for each subject, providing, at minimum, the file paths for the dynamic spectroscopy twix file, 1-point dixon dissolved phase imaging twix file, and .mat file containing the variable `mask_reg` defining the thoracic cavity volume (128x128x128 size)
2. Run the script  [GX_Siemens_recon_keyhole_script.m](GX_Siemens_recon_keyhole_script.m) and select the config file created

## Advanced Usage

For advanced users that don't want to follow the script wrapper, feel free to edit  [GX_Siemens_recon_keyhole.m](GX_Siemens_recon_keyhole.m) and turn this into a script. This file is the "main" file that performs the processing.

## Contributing

Pull requests and creating issues are welcome.

## Contact

Please contact junlan.lu@duke.edu or lujunlan98@gmail.com for questions
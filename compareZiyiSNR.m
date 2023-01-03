gasH = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\002106_highBW\nii_gas_002106_highBW_half.nii');
gas = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\002106_highBW\nii_g_002106_highBW_ori.nii');
mask = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\002106_highBW\mask_002106_highBW.nii');
gasH = flip(flip(flip(gasH,1),2),3);

maskH = imresize3(double(mask),[64 64 64],'nearest');
maskH(maskH <= 0.75) = 0; maskH(maskH > 0.75) = 1; 
maskH = logical(maskH);

[~, gas_SNR] = imageSNR(gas, mask, 8, 0.75*8^3)
[~, gasH_SNR] = imageSNR(gasH, maskH, 8, 0.75*8^3)

figure(1)
montage(gas(:,:,34:4:100)/max(gas(:)),'size',[2 8])
figure(2)
montage(mask(:,:,34:4:100),'size',[2 8])
figure(3)
montage(gasH(:,:,17:2:50)/max(gasH(:)),'size',[2 8])
figure(4)
montage(maskH(:,:,17:2:50),'size',[2 8])

%%

gasH = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\004021_highBW\nii_dis_004021_highBW_half.nii');
gas = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\004021_highBW\nii_dis_004021_highBW_ori.nii');
mask = niftiread('D:\Elly\Downloads\half_recon_scans\half_recon_scans\004021_highBW\mask_004021_highBW.nii');
gasH = flip(flip(flip(gasH,1),2),3);

maskH = imresize3(double(mask),[64 64 64],'nearest');
maskH(maskH <= 0.75) = 0; maskH(maskH > 0.75) = 1; 
maskH = logical(maskH);

[~, gas_SNR] = imageSNR(gas, mask, 8, 0.75*8^3)
[~, gasH_SNR] = imageSNR(gasH, maskH, 8, 0.75*8^3)

figure(1)
montage(gas(:,:,34:4:100)/max(gas(:)),'size',[2 8])
figure(2)
montage(mask(:,:,34:4:100),'size',[2 8])
figure(3)
montage(gasH(:,:,17:2:50)/max(gasH(:)),'size',[2 8])
figure(4)
montage(maskH(:,:,17:2:50),'size',[2 8])

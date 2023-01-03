function volume = save3DRGB2nii(volume, file_name)
% Function to convert a 3D binned volume to a nifti stack of RGB values
    color = uint8(volume*255); % need uint8 to save to RGB
    % some fancy and tricky re-arrange
%     color = permute(color, [2, 3, 0, 1]);
%     cline = reshape(color, 1, size(color));
%     color = reshape(cline, size(volume));
%   
    color = permute(color, [1, 3, 2, 4]);
    volume = permute(volume, [1, 3, 2, 4]);
    m = make_nii(volume);
    save_nii(m, file_name);
end
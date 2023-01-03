function RGBbin = Bin2RGB(binning, index2color_map)
    
    RGBbin = zeros([size(binning, 1), size(binning, 2), size(binning, 3), 3]);

    for i = 1:size(binning,1)
        for j = 1:size(binning, 2)
            for k = 1:size(binning,3)
                % RGBbin(i, j, k, :) = index2color_map(binning(i,j,k));
                RGBbin(i, j, k, :) = index2color_map{binning(i,j,k)+1};
             end
        end
    end

end
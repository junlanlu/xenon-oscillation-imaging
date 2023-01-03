function Binned_Image = six_bin_image(Input_Image,Mask,Healthy_Mean,Healthy_STD)

Binned_Image = zeros(size(Input_Image));

% < mean - 2*std
Binned_Image(Input_Image <= Healthy_Mean - 2 * Healthy_STD) = 1;  %1 - red
% mean - 2*std --> mean - 1*std
Binned_Image(Input_Image <= Healthy_Mean - 1 * Healthy_STD &...
    Input_Image > Healthy_Mean - 2 * Healthy_STD ) = 2;           %2 - orange
% mean - std --> mean
Binned_Image(Input_Image <= Healthy_Mean & ...
    Input_Image > Healthy_Mean - 1 * Healthy_STD) = 3;            %3 - green 1 
% mean --> mean + std
Binned_Image(Input_Image <= Healthy_Mean + 1 * Healthy_STD &...
    Input_Image > Healthy_Mean) = 4;                              %4 - green 2
% mean + std --> mean + 2*std
Binned_Image(Input_Image <= Healthy_Mean + 2 * Healthy_STD &...
    Input_Image > Healthy_Mean + 1 * Healthy_STD) = 5;            %5 - blue 1
% > mean + 2*std
Binned_Image(Input_Image > Healthy_Mean + 2 * Healthy_STD) = 6;   %6 - blue 2

Binned_Image(Mask == 0) = 0;                                      %0 - black

function Binned_Image = eight_bin_image(Input_Image,Mask,Healthy_Thresh)

Binned_Image = zeros(size(Input_Image));  

% < mean - 2*std
Binned_Image(Input_Image <= Healthy_Thresh(1)) = 1;  %1 - red
% mean - 2*std --> mean - 1*std
Binned_Image(Input_Image <= Healthy_Thresh(2) &...
    Input_Image > Healthy_Thresh(1) ) = 2;           %2 - orange
% mean - std --> mean
Binned_Image(Input_Image <= Healthy_Thresh(3) & ...
    Input_Image > Healthy_Thresh(2)) = 3;            %3 - green 1 
% mean --> mean + std
Binned_Image(Input_Image <= Healthy_Thresh(4) &...
    Input_Image > Healthy_Thresh(3)) = 4;             %4 - green 2;
% mean + std --> mean + 2*std
Binned_Image(Input_Image <= Healthy_Thresh(5) &...
    Input_Image >Healthy_Thresh(4)) = 5;            %5 - green 3
% mean + 2*std --> mean + 3*std
Binned_Image(Input_Image <= Healthy_Thresh(6) &...
    Input_Image > Healthy_Thresh(5)) = 6;            %6 - purple 1
% mean + 3*std --> mean + 4*std
Binned_Image(Input_Image <= Healthy_Thresh(7) &...
    Input_Image > Healthy_Thresh(6)) = 7;            %7 - purple 2
% > mean + 4*std
Binned_Image(Input_Image > Healthy_Thresh(7)) = 8;   %8 - purple 3

Binned_Image(Mask == 0) = 0;                                      %0 - black


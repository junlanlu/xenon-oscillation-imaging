

Mask = High_Res_Gas_Mask;
Low_Res_Mask = Low_Res_Gas_Mask;
RBC_Noise = std((RBC_Tot(Low_Res_Mask==0)));

ScannerType = 'Siemens';

if contains(ScannerType,'GE')
    
    healthy_mean_Osc = 8.4997;
    healthy_std_Osc = 11.2616;

else 

   %% healthy_mean_Mean = 9.4912; - These values are for Detrending Raw Data using RBC Fit rather than dissolved phase fit (Probably not how we want to do it)
   %% healthy_std_Mean = 8.9182;
    healthy_mean_Osc = 9.8104;
    healthy_std_Osc = 9.0112;
    
end

RBC_Mask = Mask;
RBC_Mask(RBC_Tot<1.5*RBC_Noise) = 0;

Mean_Pct_Diff = (RBC_High - RBC_Low)/mean(RBC_Tot(RBC_Mask==1))*100;
Mean_Pct_Diff = Mean_Pct_Diff.*RBC_Mask;

MeanBinMap = zeros(size(Mean_Pct_Diff));    
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc - 2 * healthy_std_Osc) = 1; %1 - red
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc - 1 * healthy_std_Osc &...
    Mean_Pct_Diff > healthy_mean_Osc - 2 * healthy_std_Osc ) = 2;        %2 - orange
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc & ...
    Mean_Pct_Diff > healthy_mean_Osc - 1 * healthy_std_Osc) = 3;         %3 - green 1 
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc + 1 * healthy_std_Osc &...
    Mean_Pct_Diff > healthy_mean_Osc) = 4;                               %4 - green 2;
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc + 2 * healthy_std_Osc &...
    Mean_Pct_Diff > healthy_mean_Osc + 1 * healthy_std_Osc) = 5;         %5 - green 3
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc + 3 * healthy_std_Osc &...
    Mean_Pct_Diff > healthy_mean_Osc + 2 * healthy_std_Osc) = 6;         %6 - purple 1
MeanBinMap(Mean_Pct_Diff <= healthy_mean_Osc + 4 * healthy_std_Osc &...
    Mean_Pct_Diff > healthy_mean_Osc + 3 * healthy_std_Osc) = 7;         %7 - purple 2
MeanBinMap(Mean_Pct_Diff > healthy_mean_Osc + 4 * healthy_std_Osc) = 8;  %8 - purple 3
MeanBinMap(Mask == 0) = 0;                                               %0 - black
MeanBinMap(Mask == 1 & RBC_Mask == 0) = 9;

%Nine Bin Map for RBC Oscillations
EightBinMap_Osc = [0 0 0
             1 0 0  %Red
             1 0.7143 0  %Orange
             0.4 0.7 0.4  %Green1
             0 1 0  %Green2
             184/255 226/255 145/255  %Green3
             243/255 205/255 213/255  %Light Pink
             225/255 129/255 162/255  %Med Pink
             197/255 27/255 125/255  %Dark Pink
             1 1 1 % White (For regions of RBC Defect)
             ];
         
MeanBinMap = flip(flip(flip(permute((MeanBinMap),[1, 3, 2]),1),2),3);
         
figure;
montage(MeanBinMap(:,:,34:3:100), 'size', [3 7], 'DisplayRange', [-0.5 8.5])
colormap(EightBinMap_Osc)
title('Binned Oscillations')


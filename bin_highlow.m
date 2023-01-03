function [High_Bin_Index, Low_Bin_Index] = bin_highlow(RBC_Smoothed, Time_Axis, Heart_Rate, TR, NProj)

%Function to take oscilatory data of varying amplitude and use peak finding
%and return the high and low bin indices
%[High_Bin_Index, Low_Bin_Index] = bin_highlow(data_smoothed, time_axis, Heart_Rate, TR, NProj);

%Find peaks again - this time for all of the data.
[maximaAll,maxlocsAll] = findpeaks(RBC_Smoothed,Time_Axis,'MinPeakDistance',0.5);
[minimaAll,minlocsAll] = findpeaks(-RBC_Smoothed,Time_Axis,'MinPeakDistance',0.5);

maxlocsAll = maxlocsAll(maximaAll>0);
minlocsAll = minlocsAll(minimaAll>0);

%Now, based on TR and HR, decide how many points to use around each peak
%Previously, I used 20% (using a sloppy binning method), so let's use 20%
%again here
HR_Bit = 0.2*(60/Heart_Rate); %time required (in s) for 20% of one heart beat
NPts = HR_Bit/TR;
%Force NPts to be even to make this easy
if mod(floor(NPts),2) ==0
    NPts = floor(NPts);
else
    NPts = ceil(NPts);
    if mod(floor(NPts),2) ~=0
        NPts = NPts + 1;
    end
end
%Now - we've found maxima and we know how many points we need. Just need to
%get those points now 
%Probably better ways to do this, but this should work
High_Bin_Index = [];
for i = 1:length(maxlocsAll)
    idx = find(Time_Axis==maxlocsAll(i));
    tmp_indx = (idx-NPts/2):(idx+NPts/2);
    High_Bin_Index = [High_Bin_Index tmp_indx];
end

Low_Bin_Index = [];
for i = 1:length(minlocsAll)
    idx = find(Time_Axis==minlocsAll(i));
    tmp_indx = (idx-NPts/2):(idx+NPts/2);
    Low_Bin_Index = [Low_Bin_Index tmp_indx];
end

%We might occasionally go beyond the indices allowed... fix that here
High_Bin_Index(High_Bin_Index<1) = [];
Low_Bin_Index(Low_Bin_Index<1) = [];
High_Bin_Index(High_Bin_Index>NProj) = [];
Low_Bin_Index(Low_Bin_Index>NProj) = [];
%Let's force these to have the same number of points - pretty easy to do in
%this case - pull points from the end since those will be low SNR anyway
if length(High_Bin_Index)>length(Low_Bin_Index)
    High_Bin_Index = High_Bin_Index(1:length(Low_Bin_Index));
elseif length(High_Bin_Index)<length(Low_Bin_Index)
    Low_Bin_Index = Low_Bin_Index(1:length(High_Bin_Index));
end

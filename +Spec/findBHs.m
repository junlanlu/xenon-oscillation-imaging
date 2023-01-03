function [BHstart, BHend] = findBHs(t, BHs)

BHstart = find(round(t,2) == BHs(1));
BHend = find(round(t,2) == BHs(2));
    
if isempty(BHend)
    BHend = length(t);
end 

if isempty(BHstart)
    BHstart = find(round(t,1) == BHs(1));
    BHstart = round(median(BHstart),1);
end
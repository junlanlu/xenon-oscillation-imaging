function [yyyy, mm, dd] = getScanDate(twix_obj)
% Read in twix or P file and get date of scan aquisition

scanDate = twix_obj.hdr.Phoenix.tReferenceImage0;   
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
yyyy = scanDate(1:4);
mm = scanDate(5:6);
dd = scanDate(7:8);

end 

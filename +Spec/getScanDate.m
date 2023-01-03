function [yyyy, mm, dd] = getScanDate(raw_path)

[~, ~, scanner] = fileparts(raw_path);

% Read in twix or P file and get date of scan aquisition
if strcmp(scanner,'.dat')
    twix = mapVBVD(raw_path);
    if isfield(twix.hdr.Meas, 'tFrameOfReference')
    % Duke Siemens File   
        scanDate = twix.hdr.Meas.tFrameOfReference;   
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{11};
        yyyy = scanDate(1:4);
        mm = scanDate(5:6);
        dd = scanDate(7:8);
    elseif isfield(twix.hdr.MeasYaps,'tReferenceImage1')
    % UVA Siemens File
        scanDate = twix.hdr.MeasYaps.tReferenceImage1;
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{10};
        yyyy = scanDate(1:4);
        mm = scanDate(5:6);
        dd = scanDate(7:8);
    end 
elseif strcmp(scanner,'.7')
    pfile = GE.Pfile.read(raw_path);
    if(pfile.exam.ex_datetime == 0)
        %Create exam and series dates in YYYY_MM_DD format
        yyyy = ['20',pfile.rdb.rdb_hdr_scan_date(8:9)'];
        mm = pfile.rdb.rdb_hdr_scan_date(1:2)';
        dd = pfile.rdb.rdb_hdr_scan_date(4:5)';
    else
        date_number = pfile.exam.ex_datetime/86400 + datenum(1970,1,1);
        scanDate = datestr(date_number,'yyyymmdd');
        yyyy = scanDate(1:4);
        mm = scanDate(5:6);
        dd = scanDate(7:8);
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make label images
%   data_dir
%       |-rtif  rxxxx.tif (realigned stack)
%       |-labelimg  rxxxx_labelimg.mat (labelimg), 
%       |-tcourse   rxxxx_tcourse.mat (raw timecourses from cell masks)
%
%       Kenichi Ohki  05.09.08


data_dir='K:\images\mouse090727';
analysis_dir='K:\images\mouse090727';

dirs=GetDirs(data_dir,analysis_dir);
[filelist, nFiles] = GetFileList(dirs,'labelimg');

for i=1:nFiles
    rfnames=GetFnames(filelist{i},data_dir,analysis_dir);
     if exist(rfnames.tcourse)==2
         continue;
     end
     load(rfnames.labelimg);
     array=readtiff(rfnames.rdata);
     timeCourses = stackGetTimeCourses(array, labelimg);
     save (rfnames.tcourse, 'timeCourses', '-v6');
 end
 clear all

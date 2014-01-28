%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make label images
%   data_dir
%       |-rtif  rxxxx.tif (realigned stack)
%       |-labelimg  rxxxx_labelimg.mat (labelimg), 
%       |-tcourse   rxxxx_tcourse.mat (raw timecourses from cell masks)
%
%       Kenichi Ohki  05.09.08


 data_dir='E:\users\kenichi\leicatemp\';
 %tif_dir=[data_dir,'tif\'];
 rtif_dir=[data_dir,'rtifpc\'];
 mkdir(data_dir, 'tcoursepc');
 label_dir=[data_dir,'labelimgpc\'];
 tcourse_dir=[data_dir,'tcoursepc\']
 [filelist, nFiles] = GetMatFileNames(label_dir,[],'labelimg');
 for i=1:nFiles
     temp=(filelist{i});
     baseFilename=getBaseFilename(temp,'_labelimg');
     fnames=filenameManager(baseFilename);
     labelimgfname=fullfile(label_dir, fnames.labelimg);
     rtiffname=fullfile(rtif_dir, fnames.data);
     tcoursefname=fullfile(tcourse_dir, fnames.tcourse);
     if exist(tcoursefname)==2
         continue;
     end
     load(labelimgfname);
     array=readtiff(rtiffname);
     timeCourses = stackGetTimeCourses(array, labelimg);
     save (tcoursefname, 'timeCourses', '-v6');
 end
 clear all

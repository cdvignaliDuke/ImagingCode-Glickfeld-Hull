 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make avg images, in case you forget to make avg at
% realignment
%   data_dir
%       |-tif   xxxx.mat (parameters)
%       |-rtif  rxxxx.tif (realigned stack)
%
%       Kenichi Ohki  04.28.08

data_dir='K:\images\mouse090709';

dirs=GetDirs(data_dir);

[filelist, nFiles] = GetTifFileNames(dirs.rtif);


for i=1:nFiles
    rfnames=GetFnames(filelist{i},data_dir);

    if exist(rfnames.avg)==2
        continue;
    end
    
    rstack=double(readtiff(rfnames.rdata));
    
    writetiff(mean(rstack,3),rfnames.avg);
end
 
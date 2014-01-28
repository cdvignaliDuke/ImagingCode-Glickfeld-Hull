 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make label images
%   data_dir
%       |-tif   xxxx.mat (parameters)
%       |-rtif  rxxxx.tif (realigned stack)
%       |-labelimg  rxxxx_labelimg.mat (labelimg), 
%                   rxxxx_shadelabel.tif (shaded cells to see results)
%
%       Kenichi Ohki  04.28.08

data_dir='K:\images\mouse090727';
analysis_dir='K:\images\mouse090727';

UseMP = 0; 
% if UseMP = 1, get Magnification from mat files, and set cell diameter automatically
%   for stacks obtained by MP scope

template=6;
r_th=0.3;
cell_size=2;
% You need to specify these numbers, if you don't use MP scope


dirs=GetDirs(data_dir, analysis_dir);

[filelist, nFiles] = GetFileList(dirs,'rdata');

for i=1:nFiles
    rfnames=GetFnames(filelist{i},data_dir,analysis_dir);

    if exist(rfnames.labelimg)==2
        continue;
    end
    img=double(readtiff(rfnames.avg));
    
    if UseMP
      fnames=GetFnames(filelist{i},data_dir,analysis_dir,'-r');
      load(fnames.mat);
      Mag=params.Magnification;
      Mag=str2num(Mag(2:length(Mag)));
      [labelimg,hp_img]=imFindCellsTM(img,Mag*4.5);
    else
        [labelimg,hp_img]=imFindCellsTM(img,template, r_th,cell_size, 1);
    end
     save(rfnames.labelimg, 'labelimg', 'hp_img');
     im=[imShade(hp_img, logical(labelimg)),repmat(hp_img-0.5,[1,1,3])];
     imwrite(im,rfnames.shadecells);
     close all;
 end
clear all
 
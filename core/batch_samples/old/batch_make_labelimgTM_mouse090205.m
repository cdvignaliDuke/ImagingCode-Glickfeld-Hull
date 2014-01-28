 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make label images
%   data_dir
%       |-tif   xxxx.mat (parameters)
%       |-rtif  rxxxx.tif (realigned stack)
%       |-labelimg  rxxxx_labelimg.mat (labelimg), 
%                   rxxxx_shadelabel.tif (shaded cells to see results)
% get Magnification from mat files, and set cell diameter automatically
%   for stacks obtained by MP scope
%
%       Kenichi Ohki  04.28.08

 data_dir='\\zmey\storlab\data\Kenichi\mouse090205\';
 tif_dir=[data_dir,'tif\'];
 rtif_dir=[data_dir,'rtif\'];
 label_dir=[data_dir,'labelimg\'];
 mkdir(data_dir, 'labelimg');
 [filelist, nFiles] = GetTifFileNames(rtif_dir,[],'avg');
 for i=1:nFiles
     temp=(filelist{i});
     baseFilename=getBaseFilename(temp,'_avg');
     fnames=filenameManager(baseFilename);
     fnames_before_align=filenameManager(baseFilename,'-r');
     avgfname=fullfile(rtif_dir,fnames.avg);
%      matfname=fullfile(tif_dir,fnames_before_align.mat);
     labelimgfname=fullfile(label_dir,fnames.labelimg);
     shadecellsfname=fullfile(label_dir,fnames.shadecells);
     if exist(labelimgfname)==2
         continue;
     end
     img=double(readtiff(avgfname));
%      load(matfname);
%      Mag=params.Magnification;
%      Mag=str2num(Mag(2:length(Mag)));
%     [labelimg,hp_img]=imFindCellsTM(img,Mag*4.5);
    template=8;
    r_th=0.3;
    cell_size=10;
     [labelimg,hp_img]=imFindCellsTM(img,template, r_th,cell_size);
     save(labelimgfname, 'labelimg', 'hp_img');
     im=[imShade(hp_img, logical(labelimg)),repmat(hp_img-0.5,[1,1,3])];
     imwrite(im,shadecellsfname);
     close all;
 end
clear all
 
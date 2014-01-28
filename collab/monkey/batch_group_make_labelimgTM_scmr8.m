 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batch file to make label images
%   data_dir
%       |-tif   xxxx.mat (parameters)
%       |-rtif  rxxxx_groupxx.tif (realigned stack)
%       |-labelimg  rxxxx_groupxx_labelimg.mat (labelimg), 
%                   rxxxx_groupxx_shadelabel.tif (shaded cells to see results)
% get Magnification from mat files, and set cell diameter automatically
%   for stacks obtained by MP scope
%
%       Kenichi Ohki  04.28.08


filelist='exp_list_scmr8_3.xls';
s=readXLSfilelist(filelist);

groups=[s.group];
Ngroups=max(groups);

data_dir='\\zmey\storlab\data\Soumya\scmr8\';
 tif_dir=[data_dir,'tif\'];
 rtif_dir=[data_dir,'rtif\'];
 label_dir=[data_dir,'labelimg\'];
% mkdir(data_dir, 'labelimg');

for i=1:Ngroups

    run_index=find(groups==i);
    if isempty(run_index)
        continue;
    end

    groupname=['group',num2str(i)];
    for j=1:length(run_index)
        rfname=s(run_index(j)).filename;
        avgfname=fullfile(rtif_dir, [rfname,'_',groupname,'_avg.tif']);
        if j==1
            img=double(readtiff(avgfname));
        else
            img=img+double(readtiff(avgfname));
        end
    end
        
     baseFilename=s(run_index(1)).filename;
     fnames_group=filenameManager([baseFilename,'_',groupname]);
     fnames_before_align=filenameManager(baseFilename,'-r');
     matfname=fullfile(tif_dir,fnames_before_align.mat);

     labelimgfname=fullfile(label_dir,fnames_group.labelimg);
     shadecellsfname=fullfile(label_dir,fnames_group.shadecells);
     
     
     if exist(labelimgfname)==2
         continue;
     end
     
     load(matfname);

     Mag=params.Magnification;
     Mag=str2num(Mag(2:length(Mag)));
     
     [labelimg,hp_img]=imFindCellsTM(img,Mag*4.5);
     
     for j=1:length(run_index)
         baseFilename=s(run_index(j)).filename;
         fnames_group=filenameManager([baseFilename,'_',groupname]);
         labelimgfname=fullfile(label_dir,fnames_group.labelimg);
         shadecellsfname=fullfile(label_dir,fnames_group.shadecells);
         save(labelimgfname, 'labelimg', 'hp_img');
         im=[imShade(hp_img, logical(labelimg)),repmat(hp_img-0.5,[1,1,3])];
         imwrite(im,shadecellsfname);
     end
         
     close all;
 end
clear all
 

function analyze_2stacks_stats_oristats (filelist, n, data_dir, out_dir)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze a pair of stacks 
%   stack1: any stack (no ori analysis)
%   stack2: ori stack
%
%   data_dir
%       |-tif   xxxx.mat (parameters)
%       |-rtif  rxxxx.tif (realigned stack)
%   out_dir
%                   rxxxx_labelimg.mat (labelimg), 
%                   rxxxx_shadelabel.tif (shaded cells to see results)
% get Magnification from mat files, and set cell diameter automatically
%   for stacks obtained by MP scope
%
%       Kenichi Ohki  01.13.09

%filelist='exp_list_scmr6_2.xls';
%data_dir='\\zmey\storlab\data\Soumya\scmr6\';
%out_dir='E:\users\kenichi\temp';

%n(1)=25;
%n(2)=26;

s=readXLSfilelist(filelist);

tif_dir=[data_dir,'tif\'];
rtif_dir=[data_dir,'rtif\'];

% filenames for stack1
fnames1=filenameManager(s(n(1)).filename);
fnames_before_align1=filenameManager(s(n(1)).filename,'-r');
avgfname1=fullfile(rtif_dir,fnames1.avg);
matfname1=fullfile(tif_dir,fnames_before_align1.mat);
labelimgfname1=fullfile(out_dir,fnames1.labelimg);
shadecellsfname1=fullfile(out_dir,fnames1.shadecells);
rtiffname1=fullfile(rtif_dir, fnames1.data);
tcoursefname1=fullfile(out_dir, fnames1.tcourse);
statsfname1=fullfile(out_dir, fnames1.stats);

% filenames for stack2
fnames2=filenameManager(s(n(2)).filename);
fnames_before_align2=filenameManager(s(n(2)).filename,'-r');
avgfname2=fullfile(rtif_dir,fnames2.avg);
matfname2=fullfile(tif_dir,fnames_before_align2.mat);
labelimgfname2=fullfile(out_dir,fnames2.labelimg);
shadecellsfname2=fullfile(out_dir,fnames2.shadecells);
rtiffname2=fullfile(rtif_dir, fnames2.data);
tcoursefname2=fullfile(out_dir, fnames2.tcourse);
Oristatsfname2=fullfile(out_dir, fnames2.Oristats);

% read avg images
img1=double(readtiff(avgfname1));
img2=double(readtiff(avgfname2));

% align two avg images
[outs, rimg2]=stackRegister(img2,img1,1);

% get average after aligning
img=(img1+rimg2)/2;

% read matfile to determine cell size
load(matfname1);
Mag=params.Magnification;
Mag=str2num(Mag(2:length(Mag)));
template_diam=Mag*4.5;

% outline cells by template matching
[labelimg,hp_img]=imFindCellsTM(img,template_diam,0.5,template_diam,1);

% highpass filtered images and label images for two stacks
hp_filter=fspecial('gaussian',ceil(template_diam*5),template_diam);
hp_img1=img1./imFilter2(hp_filter,img1);
hp_img2=img2./imFilter2(hp_filter,img2);
labelimg1=labelimg;
labelimg2=circshift(labelimg,[-outs(3),-outs(4)]);

% save 
labelimg=labelimg1;,hp_img=hp_img1;
save(labelimgfname1, 'labelimg', 'hp_img');
labelimg=labelimg2;,hp_img=hp_img2;
save(labelimgfname2, 'labelimg', 'hp_img');

% shade cells by label images
im1=[imShade(hp_img1, logical(labelimg1)),repmat(hp_img1-0.5,[1,1,3])];
im2=[imShade(hp_img2, logical(labelimg2)),repmat(hp_img2-0.5,[1,1,3])];
imwrite(im1,shadecellsfname1);
imwrite(im2,shadecellsfname2);

% get timecourses for stack1
array1=readtiff(rtiffname1);
timeCourses1 = stackGetTimeCourses(array1, labelimg1);

% get timecourses for stack2
array2=readtiff(rtiffname2);
timeCourses2 = stackGetTimeCourses(array2, labelimg1);

%save
timeCourses=timeCourses1;
save (tcoursefname1, 'timeCourses', '-v6');
timeCourses=timeCourses2;
save (tcoursefname2, 'timeCourses', '-v6');

% get stats for stack1
[ProcTimeCourses1,data_tables1]=tcProcess(timeCourses1,s(n(1)));
stats1 = tcStats(data_tables1, 0.05);

% get Oristats for stack2
[ProcTimeCourses2,data_tables2]=tcProcess(timeCourses2,s(n(2)));
Oristats2 = tcOriStats(data_tables2, 0.05);

%save
stats=stats1;
ProcTimeCourses=ProcTimeCourses1;
save (statsfname1 ,'stats','ProcTimeCourses');
Oristats=Oristats2;
ProcTimeCourses=ProcTimeCourses2;
save (Oristatsfname2 ,'stats','ProcTimeCourses');
save (fullfile(out_dir,[fnames1.stats,'_',fnames2.Oristats]), 'stats1', 'Oristats2','ProcTimeCourses1','ProcTimeCourses2','labelimg1', 'labelimg2');

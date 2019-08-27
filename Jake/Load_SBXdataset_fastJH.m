%% load 2P imaging data
<<<<<<< HEAD:ashley/uses scanbox/Load_SBXdataset_fast.asv
SubNum = '613';
date = '150508';
time = '1446';
ImgFolder = '003';
mouse = 'AW13';
fName = '003_000_000';
=======
SubNum = '925';
date = '150613';
time = '1833';
ImgFolder = '001';
mouse = 'img25';
fName = '925_000_000';
>>>>>>> 7e0db6ea1610cbdceb24c434872ee254e411c97b:Jake/Load_SBXdataset_fastJH.m

% load MWorks file
% CD = ['Z:\Data\2P_imaging\behavior'];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);

<<<<<<< HEAD:ashley/uses scanbox/Load_SBXdataset_fast.asv
% Set current directory to imaging data location
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' ImgFolder];
=======
% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['D:\Jake_temp\two photon' '\' date '_' mouse '\' mouse '\'];
>>>>>>> 7e0db6ea1610cbdceb24c434872ee254e411c97b:Jake/Load_SBXdataset_fastJH.m
cd(CD);
% CD = ['D:\Ashley_temp\' date '\' ImgFolder];
% cd(CD);
imgMatFile = [fName '.mat'];
load(imgMatFile);

% eyeMatFile = [fName '_eye.mat'];
% load(eyeMatFile);

%%

<<<<<<< HEAD:ashley/uses scanbox/Load_SBXdataset_fast.asv
%new datasets
nframes = 13577;
% nframes = info.config.frames;
=======
% %set the number of parallel streams (I think we have up to 8 cores)
% tic % tic and toc will tell you how long matlab took to load the data.
% nframes = 9000;
% ncores = 6;
% parpool = ncores;
% %set the size of the batch that each stream will read
% nbatch = nframes./ncores;
% %create a parallel loop to load the data
% parfor i=1:ncores
%     data(:,:,:,:,i) = sbxread(fName, ((i-1)*nbatch), nbatch);
% end
% %reshape data to make a 3D stack
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
% siz = size(data);
% data = reshape (data, siz(1), siz(2), siz(3)*ncores);
% toc

%%

%new datasets
nframes = info.config.frames;
>>>>>>> 7e0db6ea1610cbdceb24c434872ee254e411c97b:Jake/Load_SBXdataset_fastJH.m
tic
data = sbxread(fName,0,nframes);
toc
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
<<<<<<< HEAD:ashley/uses scanbox/Load_SBXdataset_fast.asv
=======
figure; imagesc(squeeze(mean(data(:,:,:),3))), truesize, axis off;
>>>>>>> 7e0db6ea1610cbdceb24c434872ee254e411c97b:Jake/Load_SBXdataset_fastJH.m

%% pre-yeti data (before 14/11/12)
% 
% %set the number of parallel streams (I think we have up to 8 cores)
% tic % tic and toc will tell you how long matlab took to load the data.
% nframes = info.config.frames;
% ncores = 6;
% parpool = ncores;
% %set the size of the batch that each stream will read
% nbatch = nframes./ncores;
% %create a parallel loop to load the data
% parfor i=1:ncores
%     data(:,:,:,:,i) = sbxread1(fName, ((i-1)*nbatch), nbatch);
% end
% %reshape data to make a 3D stack
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
% siz = size(data);
% data = reshape (data, siz(1), siz(2), siz(3)*ncores);
% toc
% 

% % pre-yeti data (before 14/11/12)
% %for splitting up data
% first half
% nframes = 27000;
% half = nframes./2;
% tic % tic and toc will tell you how long matlab took to load the data.
% ncores = 6;
% parpool = ncores;
% set the size of the batch that each stream will read
% nbatch = half./ncores;
% create a parallel loop to load the data
% parfor i=1:ncores
%     data1(:,:,:,:,i) = sbxread1(fName, ((i-1)*nbatch), nbatch);
% end
% reshape data to make a 3D stack
% pmt = 1; %1 = green 2 = red
% data1 = squeeze(data1(pmt,:,:,:,:));
% siz = size(data1);
% data1 = reshape (data1, siz(1), siz(2), siz(3)*ncores);
% toc
% second half
% tic
% create a parallel loop to load the data
% parfor i=1:ncores
%     data2(:,:,:,:,i) = sbxread1(fName, ((i-1)*nbatch)+half, nbatch);
% end
% reshape data to make a 3D stack
% pmt = 1; %1 = green 2 = red
% data2 = squeeze(data2(pmt,:,:,:,:));
% siz = size(data2);
% data2 = reshape (data2, siz(1), siz(2), siz(3)*ncores);
% toc

% %put together
% data = cat(3,data1,data2);
% clear data1
% clear data2

%%

%reshape image for with resonnant scan correction
% info.S = sparseint;
% for i = 1:nframes
%     dataSquish(:,:,i) = data(:,:,i)*info.S;
% end
% figure; imagesc(squeeze(mean(data(:,:,:),3))*info.S),truesize, axis off;

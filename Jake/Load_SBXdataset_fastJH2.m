%% load 2P imaging data
SubNum = '925';
date = '151117';
time = '1833';
ImgFolder = '001';
mouse = 'thy1';
fName = 'thy1_000_004';

% load MWorks file
% CD = ['Z:\Data\2P_imaging\behavior'];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['Z:\Data\2P_imaging\151117_thy1\thy1'];
cd(CD);
imgMatFile = [fName '.mat'];
load(imgMatFile);

%%

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
tic
data = sbxread(fName,0,nframes);
toc
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
figure; imagesc(squeeze(mean(data(:,:,:),3))), truesize, axis off;

%%
writetiff(data(:,:,1:1000), 'Z:\Data\2P_imaging\151117_thy1\thy1_tiff');

data_stable = mean(data(:,:,840:930),3); 
[out data_reg] = stackRegister(data(:,:,1:2000), data_stable);



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
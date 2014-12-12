%% load 2P imaging data
SubNum = '516';
date = '141017';
time = '1026';
ImgFolder = '004';
mouse = '516';
fName = '004_000_000';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['D:\Ashley_temp' '\' date '\' ImgFolder];
cd(CD);

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

% %new datasets
% nframes = 9000;
% tic
% data = sbxread1(fName,0,nframes);
% toc
% % pmt = 1; %1 = green 2 = red
% % data = squeeze(data(pmt,:,:,:,:));
% data = squeeze(data);


% pre-yeti data (before 14/11/12)
%for splitting up data
%first half
nframes = 27000;
half = nframes./2;
tic % tic and toc will tell you how long matlab took to load the data.
ncores = 6;
parpool = ncores;
%set the size of the batch that each stream will read
nbatch = half./ncores;
%create a parallel loop to load the data
parfor i=1:ncores
    data1(:,:,:,:,i) = sbxread1(fName, ((i-1)*nbatch), nbatch);
end
%reshape data to make a 3D stack
pmt = 1; %1 = green 2 = red
data1 = squeeze(data1(pmt,:,:,:,:));
siz = size(data1);
data1 = reshape (data1, siz(1), siz(2), siz(3)*ncores);
toc
%second half
tic
%create a parallel loop to load the data
parfor i=1:ncores
    data2(:,:,:,:,i) = sbxread1(fName, ((i-1)*nbatch)+half, nbatch);
end
%reshape data to make a 3D stack
pmt = 1; %1 = green 2 = red
data2 = squeeze(data2(pmt,:,:,:,:));
siz = size(data2);
data2 = reshape (data2, siz(1), siz(2), siz(3)*ncores);
toc

%put together
data = cat(3,data1,data2);
clear data1
clear data2

%%

%reshape image for with resonnant scan correction
info.S = sparseint;
for i = 1:nframes
    dataSquish(:,:,i) = data(:,:,i)*info.S;
end
figure; imagesc(squeeze(mean(data(:,:,:),3))*info.S),truesize, axis off;
%% load 2P imaging data
SubNum = '001';
date = '140815';
time = '0937';
ImgFolder = '001';
mouse = 'AW01';
fName = '001_000_000';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['D:\Ashley_temp' '\' date '\' ImgFolder];
cd(CD);

%set the number of parallel streams (I think we have up to 8 cores)
tic % tic and toc will tell you how long matlab took to load the data.
ncores = 3;
parpool = ncores;
%set the size of the batch that each stream will read
nbatch = 600;
%create a parallel loop to load the data
parfor i=1:ncores
    data(:,:,:,:,i) = sbxread(fName, (i-1)*nbatch, nbatch);
end
%reshape data to make a 3D stack
pmt = 1; %1 = green 2 = red
data = squeeze(data(pmt,:,:,:,:));
siz = size(data);
data = reshape (data, siz(1), siz(2), siz(3)*ncores);
toc
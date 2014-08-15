%% load 2P imaging data
SubNum = '001';
date = '140815';
time = '006';
ImgFolder = '006';
mouse = 'AW01';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['D:\Ashley_temp' '\' ImgFolder];
cd(CD);

%set the number of parallel streams (I think we have up to 8 cores)
tic % tic and toc will tell you how long matlab took to load the data.
ncores = 5;
parpool = ncores;
%set the size of the batch that each stream will read
nbatch = 600;
%create a parallel loop to load the data
fName = '006_000_000';
parfor i=1:ncores
    data(:,:,:,:,i) = sbxread(fName, (i-1)*nbatch, nbatch);
end
%reshape data to make a 3D stack
pmt = 1 %1 = green 2 = red
data = squeeze(data(pmt,:,:,:,:));
siz = size(data);
reshape (data, siz(1), siz(2), siz(3)*ncores);
toc
CD = 'Y:\home\andrew\Behavior\Data';
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to imaging data location
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
% CD = ['D:\Ashley_temp\' date '\' ImgFolder];
% cd(CD);
imgMatFile = [fName '.mat'];
load(imgMatFile);

% eyeMatFile = [fName '_eye.mat'];
% load(eyeMatFile);

%%

%new datasets
% nframes = 13577;
nframes = info.config.frames;
tic
data = sbxread(fName,0,nframes);
toc
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
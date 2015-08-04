SubNum = '614';
date = '150623';
time = '1235';
ImgFolder = '005';
mouse = 'AW14';
fName = '005_000_000';

% % load MWorks file
% CD = ['Z:\data\' mouse '\mworks\' date];
% cd(CD);
% mworks = ['data-' 'i' SubNum '-' date '-' time]; 
% load (mworks);

% Set current directory to crash folder
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% 
eyedata = sbxeyemotion([fName '_eye.mat']);
eyeArea = [eyedata.Area];
eyeCentroid = 

%%
load([fName '_eye.mat'])
eyemov = squeeze(data);
writetiff(eyemov,'eyemovie')
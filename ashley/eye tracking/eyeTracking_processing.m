SubNum = '614';
date = '150811';
time = '1115';
ImgFolder = '001';
mouse = 'AW14';
fName = '001_000_000';

% % load MWorks file
% CD = ['Z:\data\' mouse '\mworks\' date];
CD = ['Y:\home\andrew\Behavior\Data'];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to crash folder
directory = 'Y:\home\ashley';
CD = [directory '\data\' mouse '\eye tracking\' date '\' ImgFolder];
cd(CD);

%% 
eyedata = sbxeyemotion([fName '_eye.mat']);
eyeArea = [eyedata.Area];
eyeCentroid = reshape([eyedata.Centroid],2,length(eyeArea));

%%
load([fName '_eye.mat'])
eyemov = squeeze(data);
writetiff(eyemov,'eyemovie');
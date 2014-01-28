%% 
global DIRS;
DIRS.data = 'h:\users\vincent\data';
rootdir = 'mouse080311';
expt = setupexpt(rootdir,'mouse080311E');
% update_images(fullfile(expt.animal,expt.name))

%% read leica data
lPath = '//zstorlab/v1msq/data/histed/mouse080311/mouse080311_estim6';
array = readtiff(lPath, [], 'estim6_t', true);
img_red = double(readtiff(lPath, [], 'pre_anatomy.*ch00.tif'));
img_green = mean(array,3);

%% cell segmentation
[params]  = imFindCellsParamsSet
params.breakup_factor= 5;
params.min_area = 20;
params.min_center=10;
params.clear_border = 1;
params.ratio_th = 0.1
params.show_flag = 1;

% from mark
params.min_area = 80;
params.win_size = 31; %8;
params.noise_radius = 1;
params.junk_radius = 3; % 1;
params.ratio_th = 0.18; % 0.25;
params.min_center = params.min_area/2; %8;
params.con_ratio = 0.6; %5.0;
params.breakup_factor = 3;
params.clear_border = 2;
params.fast_breakup = 1;
params.presmooth_size = 3; 
    
[maskNeurons maskGlia, maskNeuropil]=scriptFindCells(img_green,img_red,params);
expt.masks.neurons = maskNeurons;
expt.masks.neuropil = maskNeuropil;
expt.masks.glia = maskGlia;

%% load raw data, load
%fn = fullfile(DIRS.analysis,rootdir,expt.name,expt.registeredfn);
array = readtiff(expt.dirs.rawgreenpn);
% array = readtiff(expt.dirs.rawredpn);
whos array
array(:,:,1:599)=[];
[w,h,nframes]=size(array);

expt = setupstim(expt,0,600,1,nframes);

array(:,:,expt.dur+1:end)=[];

%% calculate time courses 
raw = stackGetTimeCourses(array, expt.masks.neurons);
tcNeurons = scriptProcessTimeCourses(raw,expt.trialdur);
raw = stackGetTimeCourses(array, expt.masks.glia);
tcGlia = scriptProcessTimeCourses(raw,expt.trialdur);
raw = stackGetTimeCourses(array, expt.masks.neuropil);
tcNeuropil = scriptProcessTimeCourses(raw,expt.trialdur);
[tcNeurons.neuropilCorrect] = tcProjectAndSubtract(tcNeurons.df, tcNeuropil.df);
[tcNeurons.neuropilCorrectTrial] = tcTrialAverage(tcNeurons.neuropilCorrect, expt.trialdur);

%% figure time courses
figure;
subplot(3,1,1);
imagesc(tcNeurons.trialdf',[-10 20]);colorbar
title('Neurons Trial Average');
subplot(3,1,2);
delta = tcNeurons.trialdf-repmat(tcNeuropil.trialdf,1,tcNeurons.ncells);
imagesc(delta',[-20 20]);colorbar
title('Neurons Minus Neuropil');
subplot(3,1,3);
plot(mean(tcNeurons.trialdf,2),'g');
hold on
plot(tcNeuropil.trialdf,'k');
plot(mean(tcGlia.trialdf,2),'r');
plot(mean(tcNeurons.neuropilCorrectTrial,2),'c');h = colorbar;set(h,'visible','off')
ylabel('\DeltaF');
grid;legend('Neurons','Neuropil','Glia','Subtracted')



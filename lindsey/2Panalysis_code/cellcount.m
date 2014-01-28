%% set path name

pn = 'G:\users\lindsey\analysisLG\100112\BGcube_site1_GFP_pt5scale.tif';
anal_pn = 'G:\users\lindsey\analysisLG\100112\';

%% load stack in memory
stack_GFP = readtiff(pn);
stacks_GFP = single(stack_GFP);

%% adapt stacks
img_GFP_adapt = imScale(stackLocalContrastAdaptation(stacks_GFP, 30, 1));

%% find cell masks
bwimgcell = imCellEditInteractive_3DMA_090210(stacks_GFP,[],[],2);
mask_GFPneuron = bwlabeln(bwimgcell);
save(fullfile(anal_pn,'GFPneuronmasks_site1.mat'),'mask_GFPneuron');

%% reload masks
load 'g:\users\lindsey\analysisLG\100112\GFPneuronmasks_site1.mat';


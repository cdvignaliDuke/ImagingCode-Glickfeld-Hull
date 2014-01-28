%% set path name

pn_invivo = 'G:\users\lindsey\analysisLG\091215\in vivo\091215_invivo_RGB_flipped.tif';
pn_invitro = 'G:\users\lindsey\analysisLG\091215\in vitro\091215_invitro_red_um_pt5scale_8bit.tif';
anal_pn = 'G:\users\lindsey\analysisLG\091215\';

%% load stack in memory
stack_red_invivo = readtiff(pn_invivo);
stacks_red_invivo = single(stack_red_invivo);

stack_red_invitro = readtiff(pn_invitro);
stacks_red_invitro = single(stack_red_invitro);

%% adapt stacks
img_invivo_adapt = imScale(stackLocalContrastAdaptation(stacks_red_invivo, 30, 1));
img_invitro_adapt = imScale(stackLocalContrastAdaptation(stacks_red_invitro, 30, 1));
%% find cell masks
bwimgcell_invivo = imCellEditInteractive_3DMA_SY('G:\users\lindsey\analysisLG\091215\in vivo\091215_invivo_RGB_flipped.tif',[],[],2);
mask_neuron_invivo = bwlabeln(bwimgcell_invivo);
save(fullfile(anal_pn,'invivo_neuronmasks.mat'),'mask_neuron_invivo');

bwimgcell_invitro = imCellEditInteractive_3DMA_SY(stacks_red_invitro,[],[],2);
mask_neuron_invitro = bwlabeln(bwimgcell_invitro);
save(fullfile(anal_pn,'invitro_neuronmasks.mat'),'mask_neuron_invitro');

%% reload masks
load 'g:\users\lindsey\analysisLG\091215\invivo_neuronmasks.mat';
load 'g:\users\lindsey\analysisLG\091215\invitro_neuronmasks.mat';

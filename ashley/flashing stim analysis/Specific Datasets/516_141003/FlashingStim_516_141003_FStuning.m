%% flashing stim tuning matrix - matrix contains a mask_cell for cells that 
%respond to each orientation change, cells that respond to base grating,
%cells that respond to auditory change, cells that respond to base
%auditory; table contains details about each part of the matrix
%Note: only success trials are included
SubNum = '516';
date = '141003';
mouse = '516';
ImgFolder = '002+003+004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);

FStuning(1).experiment = 'Flashing Stim - Behavior';
FStuning(1).maskLabels = ['baseVisual' 'd64' 'd50' 'd40' 'd32' 'd25' 'd20' 'd16' 'd12' 'baseAuditory' 'dAuditory'];
FStuning(1).mask_cell = zeros(264,1250,(size(FStuning(1).maskLabels,2)));

FStuning(2).experiment = 'Direction Tuning';
FStuning(2).maskLabels = ['ori0' 'ori45' 'ori90'];
FStuning(2).mask_cell = zeros(264,1250,(size(FStuning(2).maskLabels,2)));

%% 
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(CD);
load('mask&TC3VisStimSuccess.mat');
FStuning(1).mask_cell(:,:,1) = mask_cell3VisStim;

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'\'Rsp2OriChange'];
cd(CD);
load('mask&TCTarget.mat');
FStuning(1).mask_cell(:,:,2:9) = mask_cellTarget;



clear all
close all

rc = behavConstsAV;
awFSAVdatasets_V1

ncells = 0;

for iexp = 1:size(expt,2)

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
down = 10;

fnbx = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,'final_mask.mat');

load(fnbx)

nc = length(unique(mask_cell(:)))-1;

ncells = ncells+nc;

end
% this script is to draw the distribution of raw fluorescence across cells
% for each single session. Doing this because we want to compare if the raw
% data is a lot different between the new 2P and the old one.
% use TCave. do average across time for each cell, and draw a distribution.

%% Section I: set paths and create analysis folders for each session
clear;
sessions = '190903_img1023'; 
days = '1023-190903_1';
image_analysis_base = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';
%image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'b','r','k','c'};
image_analysis_dest = [image_analysis_base, sessions, '\'];
%load data-- use 
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;

%%
%raw F: average across time for each cell
TC_cells = mean(TCave);
figure; hist(TC_cells,100);
xlabel('average fluorescence across time'); ylabel('number of cells');
title(['distribution of raw fluorescence' sessions]);
savefig([image_analysis_dest sessions '_hist_rawF']);

%bottom 10%, but what if the bottom 10% is always during running?
TCave_sort = sort(TCave,1,'ascend');
TCave_base = TCave_sort(1:ceil(length(TCave_sort)*0.1),:);
TCave_baseCells = mean(TCave_base);
figure; hist(TCave_baseCells,100);
xlabel('average 10% bottom fluorescence across time'); ylabel('number of cells');
title(['distribution of bottom fluorescence' sessions]);
savefig([image_analysis_dest sessions '_hist_bottomF']);

%signal-noise(mean/std)
TC_cells_std = std(TCave);
rawF_SNR = TC_cells./TC_cells_std;
figure; histogram(rawF_SNR,1:0.5:7);
xlabel('rawF SNR'); ylabel('number of cells');
title(['distribution of rawF SNR' sessions]);
savefig([image_analysis_dest sessions '_hist_rawF_SNR']);







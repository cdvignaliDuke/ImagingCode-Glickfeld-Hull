% load previous deconvolved session and delete cells manually. 
% some of the cells that have a FR above 0.4 still looks wierd. delete
% those cell by hand. For the follwoing sessions, just do it when doing
% deconvolution session by session

clear;
sessions = '190706_img1033_vermisVI'; 
days = '1033-190706-2131_1';
%image_analysis_base    = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on isilon in the movingDots analysis folder
color_code = {'b','r','k','c'};

%%
% the names (expecially the ones ends with '_cl') are confusing in this
% script, just double check what they are if you change anything.
image_analysis_dest = [image_analysis_base, sessions, '\'];
threshold = -4;
filename = dir([image_analysis_dest '*' num2str(threshold) '_TCave_cl.mat']);
TCoutput = load([image_analysis_dest filename.name]);
TCave = TCoutput.TCave_cl;% the TCave_cl here is after deleting cells that has a FR lower than 0.4
badPCs = TCoutput.badPCs;
FRthres = TCoutput.threshold;

filename2 = dir([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat']);
spkoutput = load([image_analysis_dest filename2.name]);
FRstay_cell = spkoutput.FRstay_cell_cl;
aveFR = spkoutput.aveFR_cl;
spk_logic = spkoutput.spk_logic_cl;
spk = spkoutput.spk_cl;
kernel = spkoutput.kernel_cl;
spk_peak = spkoutput.spk_peak_cl;
spk_inx = spkoutput.spk_inx_cl;

% load TCave to get initial TCaves, this tells you the total number of
% cells.
filename3 = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave_initial = load([image_analysis_dest 'getTC\' filename3.name]);
TCave_initial = TCave_initial.tc_avg;
ntotal = size(TCave_initial,2);
allcells_inx = 1:ntotal;
goodcells1 = setdiff(allcells_inx,badPCs); % cells that are left after delete the cells with FR<0.4
badPCs_Idelete = [9,10,34];
for i = 1: length(badPCs_Idelete)
    badPCs_Idelete_newinx(i) = find(badPCs_Idelete(i) == goodcells1); % this index tells you which ones to delete in the cleaned variables (variables end with '_cl' in "deconvolution_Jin_multisession") after deconvolution
end

TCave_cl = TCave;TCave_cl(:,badPCs_Idelete_newinx) = [];
FRstay_cell_cl = FRstay_cell;FRstay_cell_cl(badPCs_Idelete_newinx) = [];
aveFR_cl = mean(FRstay_cell_cl);
spk_logic_cl = spk_logic; spk_logic_cl(:,badPCs_Idelete_newinx) = [];
spk_cl = spk; spk_cl(:,badPCs_Idelete_newinx) = [];
kernel_cl = kernel; kernel_cl(:,badPCs_Idelete_newinx) = [];
spk_peak_cl = spk_peak; spk_peak_cl(badPCs_Idelete_newinx) = [];
spk_inx_cl = spk_inx; spk_inx_cl(badPCs_Idelete_newinx) = [];

hist_FR = figure;
histogram(FRstay_cell,'BinWidth',0.1);
title([sessions 'deconvolution-firing rate-stationary-threshold' num2str(threshold)]);
saveas(hist_FR,[image_analysis_dest sessions '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);

save([image_analysis_dest sessions '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'],...
    'TCave_cl','threshold','badPCs','badPCs_Idelete');
save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'FRstay_cell_cl', ...
    'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl','spk_inx_cl','FRthres','badPCs_Idelete');

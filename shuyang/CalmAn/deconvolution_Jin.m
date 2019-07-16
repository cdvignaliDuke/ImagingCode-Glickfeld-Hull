% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.

%% assign document paths and experimental sessions
clear;
sessions = '190429_img1027'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1027-190429_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%load data
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;

behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);

%% SECTION I deconvolution
% input: fluorescence trace: T*1, threshold for spike size (how big the spikes are supposed to be: normal set as -3(3 sigma, 3 times bigger than noise levels)).
% kernal: denoised trace, T*1 vector
% spikes: all of the values in a spike(thus, you will get more than 1 values if a spike is longer than a frame. e.g.: if a spike is 120ms long, you will get 4 values back to back which all belong to a single spike)
frames = 1:1:size(TCave,1);
threshold = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary

% write a loop, get the c and s for all cells and draw a GUI
kernel = zeros(size(TCave,1),size(TCave,2));
spk = zeros(size(TCave,1),size(TCave,2));
spk_peak = {};
spk_inx = {};
spk_logic = zeros(size(TCave,1),size(TCave,2));
num_spks_cell = zeros(1,size(TCave,2));
FRstay_cell = zeros(1,size(TCave,2));
for c = 1: size(TCave,2)
    [kernel(:,c), spk(:,c), options] = deconvolveCa(TCave(:,c), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold);
    % get only the peaks of each spike
    [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
    % spike logic
    spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
    num_spks_cell(c) = sum(spk_logic(frm_stay,c)==1);
    FRstay_cell(c)= num_spks_cell(c)/length(frm_stay)*30; % firing rate = # of spikes/duration(s)
end
aveFR = mean(FRstay_cell);
hist_FR = figure; 
histogram(FRstay_cell,'BinWidth',0.1);
title([sessions 'deconvolution-firing rate-stationary-threshold' num2str(threshold)]);
saveas(hist_FR,[image_analysis_dest sessions '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);

% there are a few cells has super low FR !!!!!

%save variables
save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'threshold',...
    'FRstay_cell', 'aveFR','num_spks_cell','options','spk_logic','spk','kernel','spk_peak','spk_inx');


%% plots
% write a GUI with subplots: TCave with red dots and kernel
[fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,sessions,threshold);
savefig([image_analysis_dest sessions '_GUI_TCave_deconvolution_threshold' num2str(threshold) '.fig']);

% spike logic 
% plot spike logic
figure;
imagesc(spk_logic);
xlabel('Cell #');
ylabel('Frame #');
title ([sessions ' deconvolution']);
savefig([image_analysis_dest sessions '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold) '.fig']);

% plot spike rates for each cell during stationary
figure; 
scatter(1:length(FRstay_cell), FRstay_cell, 'bo');
ylabel('Spike Rate (Hz)');
xlabel('Cell #');
title ([sessions ' deconvolution']);
savefig([image_analysis_dest sessions '_scatter_spikeRate_deconvolve_threshold' num2str(threshold) '.fig']);


%% find the ones that has super low FR (fake cells), delete those in TCave
badPCs1 = find(FRstay_cell <= 0.4);
% some wierd fake cells doesn't fire for a long time, but then fires a lot
% within a short time (the whole trace just looks noisy and wierd) --- even if the firing rate is bigger than 0.4, raw trace still looks wierd.
for n = 1: size(TCave,2)
% if it doesn't fire for more than half a minute, then throw away? 
% if the deconvovled trace doesn't go back to baseline for a few hundred
% frames, then throw away?
    spk_inx{n}


end
TCave_cl = TCave;
TCave_cl(:,badPCs) = [];
FRstay_cell_cl = FRstay_cell;
FRstay_cell_cl(badPCs) = [];
aveFR_cl = mean(FRstay_cell_cl);

save([image_analysis_dest sessions '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'], 'TCave_cl','threshold','badPCs');
save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'FRstay_cell_cl', 'aveFR_cl','badPCs','-append');

%%
% frames = 1:1:60000;
% % what pengcheng told me:
% [c1, s1, options1] = deconvolveCa(TCave(:,1), 'optimize_pars', true, ...
%     'optimize_b', true, 'method','foopsi', 'smin', -3.5);
% % get only the peaks of each spike 
% [peaks1,locs1] = findpeaks(s1);
% % spike logic
% s_logic1 = ismember(frames,locs1);
% s_logic1 = s_logic1';
% num_spks1 = sum(s_logic1==1);
% 
% % don't know why ar2 doesn't work well, spike numbers are not reasonable,
% % 13 times bigger than using ar1.
% [c2, s2, options2] = deconvolveCa(TCave(:,1), 'foopsi','ar2',...
%     'optimize_pars', true, 'optimize_b', true,'smin', -4);
% [peaks2,locs2] = findpeaks(s2);
% s_logic2 = ismember(frames,locs2);
% s_logic2 = s_logic2';
% num_spks2 = sum(s_logic2==1);
% 
% [c3, s3, options3] = deconvolveCa(TCave(:,1), 'optimize_pars', true, ...
%     'optimize_b', true, 'method','foopsi');
% s_logic3 = logical(s3);
% num_spks3 = sum(s_logic3==1);
% 
% %%
% plotred = TCave(:,1).*s_logic2;
% plotred(plotred==0) = NaN;
% figure; 
% subplot(2,1,1)
% plot(TCave(:,1)); hold on;
% plot(plotred,'ro','MarkerSize',8); ylabel('TCave');
% xlim([200,3500]);
% subplot(2,1,2)
% plot(c2);
% xlim([200,3500]);ylabel('c');xlabel('frame');


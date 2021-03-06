% spike identification using the standard derivations of the first
% derivative of raw traces. 
% function twoP_best_spkthreshold tells you what's the optimal std(i.e., 1 std ot 2 or 3), 
% then section II draws the GUI showing raw traces with identified events
% and first derivatives.
%% assign document paths and experimental sessions
clear;
% sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'};%,'200316_img1064_airpuff_2'};
% image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% % behavior analysis results
% days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1'};

sessions = '191115_img1042'; 
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';
image_analysis_dest = [image_analysis_base, sessions,'\'];
image_analysis_dest_deriv = [image_analysis_base, sessions, '\derivative\'];
if ~exist(image_analysis_dest_deriv)
    mkdir(image_analysis_dest_deriv);
end

% behavior analysis results 
days = '1042-191115_1';
behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I  SPIKE RATES
% get derivatives from raw fluorescence
TCave = load([image_analysis_dest sessions '_deconvolution_thresh-4_TCave_cl']); % use only good cells delected by deconvolution
TCave = TCave.TCave_cl;
deriv = diff(TCave,1,1); % derivatives shoud have positive and negative values
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
airpuffon = behav_output.airpuffon1;
% find the frames during stationary without airpuff (get rid of frames 300ms after airpuff onset)
% 1.stationary without airpuff
airpuff_period = [];
for a = 1:length(airpuffon)
    airpuff_period = cat(2,airpuff_period,airpuffon(a):airpuffon(a)+9);%300ms after airpuff onset is 10ish frames
end
stay_noairpuff = setdiff(frm_stay,airpuff_period);
%  decide the threshold for spikes,  get average number of spikes during stationary (bestFR)
[aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (deriv, stay_noairpuff,TCave);
savefig([image_analysis_dest_deriv sessions '_hist_stayFR_wdiffthresh']);
% 
% % plot spike logic
% figure;
% imagesc(spk_bi_cellmat);
% xlabel('Cell #');
% ylabel('Frame #');
% title (sessions);
% savefig([image_analysis_dest sessions '_scatter_spikeLogic']);
% 
% % plot spike rates for each cell during stationary
% figure; 
% scatter(1:length(FRstay_cells), FRstay_cells, 'bo');
% ylabel('Spike Rate (Hz)');
% xlabel('Cell #');
% title (sessions);
% savefig([image_analysis_dest sessions '_scatter_spikeRate']);

save([image_analysis_dest_deriv sessions,'_spikes.mat'],'aveFRsWstds','best_thres'...
    ,'bestFR','std_best','spk_inx','spk_bi_cellmat');

%% SECTION II plots to go back and look at the raw data

% % plot first derivatives with horizontal std lines --------------------------------------------------
% 
% listfrm = [1:300; 3301:3600;6801:7100;12051:12350;19001:19300;28501:28800; 16881:17180;  ...
%   21661:21960; 26001:26300 ];
% listfrm = listfrm';
% [fig_deriv] = GUI_w_hline(deriv, listfrm,std_deriv,std2,std2_5,std3);
% savefig([image_analysis_dest sessions '_GUIderiv_stdlines']);
% 
% % after function best_spkthreshold_2P, go back to raw fluorescence(TCave) and see if the time points of spikes make sense
% % plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
% listfrm = [1:300; 3301:3600;6801:7100;12051:12350;19001:19300;28501:28800;   16881:17180;  ...
%   21661:21960; 26001:26300 ];
% listfrm = listfrm';
% [fig] = GUI_rawTrace_wSpkevents(TCave,listfrm,spk_inx,sessions);  
% savefig([image_analysis_dest sessions '_GUI_TCave_w_spkEvent_stdbest']); 
% 
% % make GUI with suplots: raw fluorescence with red circles and first derivatives
% % with std lines
% 
% [fig] = GUI_rawTrace_nderivThresh(TCave,deriv,spk_bi_cellmat,sessions);
% savefig([image_analysis_dest sessions '_GUI_TCave_w_deriv']); 
% 
% % Plot raw fluorescence for the running triggered average windows -----
% [fig_runtrig] = GUI_w_vline(TCave,frms_runTrig_mat,spk_inx); %!!!!! you will get an error if you have less than 9 running windows
% savefig([image_analysis_dest sessions '_GUIraw_fluo_RUNtrig']);

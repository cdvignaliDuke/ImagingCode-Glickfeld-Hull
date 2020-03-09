%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191114_img1040'};
days = {'1040-191114_1'};
%sessions = {'190923_img1023','190903_img1025','190906_img1021','190908_img1024',...
 %   '190911_img1026','190911_img1029','190912_img1034','190923_img1023'};
%days = {'1023-190923_1','1025-190903_1','1021-190906_1','1024-190908_1',...
 %   '1026-190911_1','1029-190911_1','1034-190912_1','1023-190923_1'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder

%image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'c','r','y','g'};

%% plot PSTH for airpuff and determine time length want to use for calculating SNR
%%SNR = aveFR after airpuff delivery/aveFR before airpuff delivery, plot a bar graph for SNR during running vs stay
for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];    
    threshold = -4;
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    FR_airstimAve_stay = spk_deconv_output.FR_airstimAve_stay;
    FR_airstimAve_run = spk_deconv_output.FR_airstimAve_run;
    % bin every 3 frames -> 0.1s, calculate average FR during every 0.1s
    n = 3;% number of frames want to bin together
    vecEnd = n*floor(length(FR_airstimAve_run)/n);% vector length can't be divided by 3, floor to the nearest number that can be divided by 3.
    FRbin_airstimAve_stay1 = reshape(FR_airstimAve_stay(1:vecEnd),[n,30/n]);
    FRbin_airstimAve_stay = mean(FRbin_airstimAve_stay1);
    FRbin_airstimAve_run1 = reshape(FR_airstimAve_run(1:vecEnd),[n,30/n]);
    FRbin_airstimAve_run = mean(FRbin_airstimAve_run1);
    
    PSTH = figure;
    x = (0:0.1:1);
    subplot(1,2,1);
    histogram('BinEdges',x, 'BinCounts', FRbin_airstimAve_stay);
    vline(0.5,'r');
    title(['airpuff stim stay ' sessions{ii}]);  
    ylabel('firing rate');xlabel('time(s)');
    %ylim([0 2]);

    subplot(1,2,2)
    histogram('BinEdges',x, 'BinCounts', FRbin_airstimAve_run);
    xlabel('time(s)'); title('airpuff stim run');
    ylim([0 2]);vline(0.5,'r');
    saveas(PSTH, [image_analysis_dest '\' days{ii} '_airpuff_PSTH']);    
 
% calculate SNR: SNR = total # of spikes during 200ms before airpuff - total # of spikes after airpuff/0.5(std response^2+std baseline^2)
% didn't transfer total # of spikes into firing rate because the std will be a lot bigger (during running, some trials have only 1 spike during the whole 200ms, others have none. if transfer into FR, 1*5 = 5 
% for a single cell in multiple trials,it looks like: [0 0 0 0 0 5 0 0 0 0]. the variation is very big. 
% but if we do totla # of spikes, it is [0 0 0 0 1 0 0 0 0]), and the std is >0 and <1. 
% when averaging across trials, the spk probablity of each frame is >0 and <1. even if *5, still >0 and <1. 
% so if use aveFR, a small number is going to be divided by a big number,which will give us a small SNR.
% also I don't think it makes sense to infer what happens in 1 second by looking at a 200ms window. like you can calculate firing rate with data
% of multiple seconds, but not when you're looking at a window shorter than 1s

% choose 200ms as the time window (6 frames before and after airpuff)
    spk_airpuff_run_mat = spk_deconv_output.spk_airpuff_run_mat; % trials*frames*cells
    spk_airpuff_stay_mat = spk_deconv_output.spk_airpuff_stay_mat; % trials*frames*cells
    
    spk_airpuff_run_cells = spk_deconv_output.spk_airpuff_run_cells; % frames*cells
    spk_airpuff_stay_cells = spk_deconv_output.spk_airpuff_stay_cells;
    
    Pspk_bf_air_run = sum(spk_airpuff_run_cells(10:15,:)); 
    Pspk_aft_air_run = sum(spk_airpuff_run_cells(16:21,:));
    Pspk_befo_air_stay = sum(spk_airpuff_stay_cells(10:15,:));
    Pspk_aft_air_stay = sum(spk_airpuff_stay_cells(16:21,:));
   
%     % aveFR during200ms = sum of spike probablity for all frames during those 200ms/0.2; 
%     % this is the same as calculating the average firing rate for each frame first and then average all 6 frames during that 200ms
%     FR_run_befo_cells = sum(Pspk_bf_air_run,1)/0.2;% /0.2 = *30/6
    
% calculte standard deviation for each cell (trial by trial variation): I'm sure all of the calculation is right!
    % average across frames first: trial*cells --> std of trials
    ave_spk_run_befo_trials = squeeze(sum(spk_airpuff_run_mat(:,10:15,:),2));%after sum you get the total # of spks of 200ms for each trial each cell
    std_spk_run_befo_trials = std(ave_spk_run_befo_trials,0,1);
    
    ave_spk_run_aft_trials = squeeze(sum(spk_airpuff_run_mat(:,16:21,:),2));
    std_spk_run_aft_trials = std(ave_spk_run_aft_trials,0,1);
    
    ave_spk_stay_befo_trials = squeeze(sum(spk_airpuff_stay_mat(:,10:15,:),2));
    std_spk_stay_befo_trials = std(ave_spk_stay_befo_trials,0,1);
    
    ave_spk_stay_aft_trials = squeeze(sum(spk_airpuff_stay_mat(:,16:21,:),2));
    std_spk_stay_aft_trials = std(ave_spk_stay_aft_trials,0,1);
% SNR for each cell
    SNR_run_cells = (Pspk_aft_air_run - Pspk_bf_air_run)./sqrt((0.5*(std_spk_run_aft_trials.^2 + std_spk_run_befo_trials.^2)));
    SNR_stay_cells = (Pspk_aft_air_stay - Pspk_befo_air_stay)./sqrt((0.5*(std_spk_stay_aft_trials.^2 + std_spk_stay_befo_trials.^2)));
    % for SNR_run_cells = NaN, make it to 0
    SNR_run_cells(isnan(SNR_run_cells))=0;

% SNR across cells
    SNR_run = mean(SNR_run_cells);
    SNR_stay = mean(SNR_stay_cells);
    ste_SNR_run = std(SNR_run_cells)/sqrt(length(SNR_run_cells));
    ste_SNR_stay = std(SNR_stay_cells)/sqrt(length(SNR_stay_cells));
    
    SNR = figure; hold on;
    for i=1:length(SNR_run_cells)
        plot([1 3], [SNR_run_cells(i) SNR_stay_cells(i)], 'o-k');
    end
    e=errorbar([1 3], [SNR_run SNR_stay],[ste_SNR_run,ste_SNR_stay],...
        'or', 'Markerfacecolor','r'); 
    xticks([1 3]); xticklabels({'Run', 'Stay'}); xlim([0 4]);
    title([sessions{ii} 'airpuff SNR']); axis square;
    saveas(SNR, [image_analysis_dest '\' days{ii} '_airpuff_SNR']); 
    figure;hist(SNR_stay_cells);
    xlim([-1,1]); axis square;hold on;
    vline(0,'r');
    title('distribution of SNR of all cells during stataionary');
    

end

%% 
%

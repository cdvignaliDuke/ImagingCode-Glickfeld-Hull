% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.

%% paths
clear;
%sessions = {'190924_img1023'}; 
sessions = {'191206_img1038'}; 
%days = {'1023-190924_1'};
days = {'1038-191206_1'};
%image_analysis_base    = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on isilon in the movingDots analysis folder
color_code = {'b','r','k','c'};

%% deconvolution and plot
for i = 1:size(sessions,2)
    % paths and load data-------------------------------------------------------------------------------------------------------
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    % behavior analysis results
    %behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{i} '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    %load data
    filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
    TCave = load([image_analysis_dest 'getTC\' filename.name]);
    TCave = TCave.tc_avg;    
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frm_stay_cell = behav_output.frames_stay_cell;
    frm_stay = cell2mat(frm_stay_cell);
    airpuffon = behav_output.airpuffon1;
    % find the frames during stationary without airpuff (get rid of frames 300ms after airpuff onset)
    % 1.stationary without airpuff 
    airpuff_period = [];
    for a = 1:length(airpuffon)
        airpuff_period = cat(2,airpuff_period,airpuffon(a):airpuffon(a)+9);%300ms after airpuff onset is 10ish frames
    end
    stay_noairpuff = setdiff(frm_stay,airpuff_period);% find the frames in frame_stay but not in airpuff_period
   
    %Deconvolution------------------------------------------------------------------------------------------------------------------
    % input: fluorescence trace: T*1, threshold for spike size (how big the spikes are supposed to be: normal set as -3(3 sigma, 3 times bigger than noise levels)).
    % kernal: denoised trace, T*1 vector
    % spikes: all of the values in a spike(thus, you will get more than 1 values if a spike is longer than a frame. e.g.: if a spike is 120ms long, you will get 4 values back to back which all belong to a single spike)
    
    frames = 1:1:size(TCave,1);
    threshold = -4; %changing the threhold basically changes the identification of spikes when the peak amplitude is small (those small peaks), doesn't change anything with the bigger jittered ones. -3 or -3.5 gives a FR close to 1 during stationary
    
    kernel = zeros(size(TCave,1),size(TCave,2));
    spk = zeros(size(TCave,1),size(TCave,2));
    spk_peak = {};
    spk_inx = {};
    spk_logic = zeros(size(TCave,1),size(TCave,2));
    num_spks_cell = zeros(1,size(TCave,2));
    FRstay_cell = zeros(1,size(TCave,2));
    % for airpuff experiment, exlcude frames that airpuff happened and then calculate firing rate during stationary
    for c = 1: size(TCave,2)
        [kernel(:,c), spk(:,c), options] = deconvolveCa(TCave(:,c), 'optimize_pars', true, ...
            'optimize_b', true, 'method','foopsi', 'smin', threshold);
        % get only the peaks of each spike
        [spk_peak{c},spk_inx{c}] = findpeaks(spk(:,c));
        % spike logic
        spk_logic(:,c) = (ismember(frames,spk_inx{c}))';
        num_spks_cell(c) = sum(spk_logic(stay_noairpuff,c)==1);
        FRstay_cell(c)= num_spks_cell(c)/length(stay_noairpuff)*30; % firing rate = # of spikes/duration(s)
    end
    aveFR = mean(FRstay_cell);
    hist_FR = figure;
    histogram(FRstay_cell,'BinWidth',0.1);
    title([sessions{i} 'deconvolution-firing rate-stationary-threshold' num2str(threshold)]);
    saveas(hist_FR,[image_analysis_dest sessions{i} '_histFR_deconvolution_threshold' num2str(threshold) '.fig']);
    
    %save variables
    save([image_analysis_dest sessions{i} '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'threshold',...
        'FRstay_cell', 'aveFR','options','spk_logic','spk','kernel','spk_peak','spk_inx');
    
    % plots---------------------------------------------------------------------------------------------------------------------------
    % write a GUI with subplots: TCave with red dots and kernel
    [fig_deconvolve] = GUI_rawTrace_nDeconvolve(TCave,kernel,spk_logic,sessions{i},threshold);
    savefig([image_analysis_dest sessions{i} '_GUI_TCave_deconvolution_threshold' num2str(threshold) '.fig']);
    
    % spike logic
    % plot spike logic
    figure; 
    colormap gray;
    imagesc(imcomplement(spk_logic'));
    xlabel('Frame #');
    ylabel('Cell #');
    %xlim([25000 28992]);
    title ([sessions{i} ' deconvolution']);
    %print([image_analysis_dest sessions{i} '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold)],'-dpdf','-fillpage');
    % this figure is before deleting the bad cells!
    savefig([image_analysis_dest sessions{i} '_scatter_spikeLogic_deconvolve_threshold' num2str(threshold) '.fig']); 
    
    % plot spike rates for each cell during stationary
    figure;
    scatter(1:length(FRstay_cell), FRstay_cell,70, 'ko');
    ylabel('Spike Rate (Hz)');
    xlabel('Cell #'); axis square;
    title ([sessions{i} ' deconvolution']);
    savefig([image_analysis_dest sessions{i} '_scatter_spikeRate_deconvolve_threshold' num2str(threshold) '.fig']);
    
    
    %% delete the bad cells
    
    %---------------------------------------------------------------------------------------------------------------------------------
    % find the ones that has super low FR (fake cells), delete those in TCave, spk_logic,spk,kernel,spk_peak,spk_inx
    FRthres_low = 0.4;
    %FRthres_high = 1.2;
    badPCs1 = find(FRstay_cell<FRthres_low);
    %     %delete bad cells manually
    %     badPCs2 = [7,9,12,18,53,88,91,101];
    
    % use other criteria to sort out bad cells:
    %e.g. session:191115_img1041
    peak_std = zeros(1,size(spk_peak,2));
    peak_mean = zeros(1,size(spk_peak,2));
    for c = 1:size(spk_peak,2)                  % for each cell
        peak_std(c) = std(spk_peak{c});
        peak_mean(c) = mean(spk_peak{c});
    end
    peak_variation = peak_std./peak_mean; %if peak variation is too big, this cell might not be good.
    figure; plot(peak_variation);
    title('peak variation');
    savefig([image_analysis_dest sessions{i} '_CaeventPeak_variation' num2str(threshold) '.fig']);
    peak_vartoomuch = find(peak_variation>0.64);
    % cannot totally through out all of the cells in here, need to go back and
    % look at each cell that has a high variation and decide.
    % but this at least gives me a good pool to look at
    % cell # 18,9,7,4,1 can all be sorted out using this way
    % ----------------------------------------------------------------------------------------------
    % another feature is that the deconvolved signal doesn't go back to
    % baseline for a long time: cell#88,53,91,and 101
    nobase_maxperiod = zeros(1,size(spk_peak,2));
    for c = 1:size(spk_peak,2)
        nobase_maxperiod(c) = max(diff(find(kernel(:,c)<100))); %kernel values smaller than 100 is considered as basline. diff(find) gives you how many continous frames the denoised signal doesn't return to baseline
    end
    figure; plot(nobase_maxperiod);
    title('not return to baseline period');
    savefig([image_analysis_dest sessions{i} '_Caevent_nobase' num2str(threshold) '.fig']);
    %if the maximum not return to baseline period is bigger than n # of frames,
    %through out this cell
    nobase_thres = 1000;
    nobase_cell = find(nobase_maxperiod > nobase_thres);
    
    %----------------------------------------------------------------------------------------------------
    %a third feature is that the cells are firing a lot in a short period of
    %time
    %calculate the firing rate of a period of time and if the FR is high for a
    %long time then through the cell out first calculate # of spikes per second
    maxfire = zeros(1,size(spk_peak,2));
    binwidth = 900; % number of frames
    bins = floor(size(spk_logic,1)/binwidth);
    for c = 1:size(spk_peak,2)
        nfire = []; t = 1;
        while t < size(spk_logic,1)-binwidth
            nfire = [nfire sum(spk_logic(t:t+binwidth-1,c)==1)];
            t = t+1;
        end
        maxfire(c) = max(nfire);
    end
    figure; plot(maxfire); title('firing too much');
    savefig([image_analysis_dest sessions{i} '_Caevent_maxfire' num2str(threshold) '.fig']);
    figure; hist(maxfire); title('histogram # of events in 900 frames');
    savefig([image_analysis_dest sessions{i} '_Caevent_maxfire_hist' num2str(threshold) '.fig']);
    FR_toohigh = find(maxfire>65);
    
    badcells2 = union(peak_vartoomuch,nobase_cell);
    badPCs2 = union(badcells2,FR_toohigh);
    
    badPCs = union(badPCs1,badPCs2);
    TCave_cl = TCave;TCave_cl(:,badPCs) = [];
    FRstay_cell_cl = FRstay_cell;FRstay_cell_cl(badPCs) = [];
    aveFR_cl = mean(FRstay_cell_cl);
    spk_logic_cl = spk_logic; spk_logic_cl(:,badPCs) = [];
    spk_cl = spk; spk_cl(:,badPCs) = [];
    kernel_cl = kernel; kernel_cl(:,badPCs) = [];
    spk_peak_cl = spk_peak; spk_peak_cl(badPCs) = [];
    spk_inx_cl = spk_inx; spk_inx_cl(badPCs) = [];
    
    save([image_analysis_dest sessions{i} '_deconvolution_thresh', num2str(threshold), '_TCave_cl.mat'], 'TCave_cl','threshold','badPCs','badPCs1','badPCs2');
    save([image_analysis_dest sessions{i} '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],'FRstay_cell_cl', ...
        'aveFR_cl','badPCs','spk_logic_cl','spk_cl','kernel_cl','spk_peak_cl','spk_inx_cl','FRthres_low','-append');
    
    %% make a new mask figure
    % load mask3D final
    filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97.5_coor0.8_mask3D.mat']);
    mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
    mask3D = mask3D.mask3D;
    allcells = 1:size(TCave,2);
    goodcells = setdiff(allcells,badPCs);
    figure; imshow([image_analysis_dest 'getTC\' 'AVG_' sessions{i} '_000_rgstr_tiff_1_29999_50_ref18_' 'jpeg.jpg']); hold on;
    for m  = 1:length(goodcells)
        bound = cell2mat(bwboundaries(mask3D(:,:,goodcells(m))));
        randcolor = rand(1,4);
        plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
        text(mean(bound(:,2)),mean(bound(:,1)), ...
            num2str(goodcells(m)), 'color', 'y', 'FontSize', 8);
    end
    hold on;
    title([sessions ' masks only good cells']);
    savefig([image_analysis_dest sessions{i} '_mask_wdendrites_goodcells.fig']);
    
    
end

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


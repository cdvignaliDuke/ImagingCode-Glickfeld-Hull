% was initially thinking about sorting out cells by using signal-noise ratio after 
% identifying spike events using first derivatives. However,this method doesn't work well 
% because the SNR has a big range across sessions, and the neurons with a lower SNR isn't necessarily "bad cells". 
% cells that look noisy in raw traces can still have a high SNR. This can be because the spike identification 
% using derivative is not accurate in the first place. 

% after get spk_inx from "twoP_best_spikethreshold", 
% calculate SNR. SNR = Peak of Ca event/SD of baseline 
% here I use the max of Ca event / SD of tc_avg during stationary.
% delete the cells that has low SNRs

%% assign document paths and experimental sessions
clear;
sessions = '190503_img1029'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1029-190503_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION III  SPIKE RATES
% get derivatives from raw fluorescence
filename = dir([image_analysis_dest 'getTC\' '*' '_TCave.mat']);
TCave = load([image_analysis_dest 'getTC\' filename.name]);
TCave = TCave.tc_avg;
deriv = diff(TCave,1,1); % derivatives shoud have positive and negative values
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
%  decide the threshold for spikes,  get average number of spikes during stationary (bestFR)
[aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (deriv, frm_stay,TCave);
savefig([image_analysis_dest sessions '_hist_stayFR_wdiffthresh']);

%% plot spike logic
figure;
imagesc(spk_bi_cellmat);
xlabel('Cell #');
ylabel('Frame #');
title (sessions);
savefig([image_analysis_dest sessions '_scatter_spikeLogic']);

% plot spike rates for each cell during stationary
figure; 
scatter(1:length(FRstay_cells), FRstay_cells, 'bo');
ylabel('Spike Rate (Hz)');
xlabel('Cell #');
title (sessions);
savefig([image_analysis_dest sessions '_scatter_spikeRate']);

save([image_analysis_dest sessions,'_spikes.mat'],'aveFRsWstds','best_thres'...
    ,'bestFR','std_best','spk_inx','spk_bi_cellmat');

% after function best_spkthreshold_2P, go back to raw fluorescence(TCave) and see if the time points of spikes make sense
% plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
listfrm = [1:300; 3301:3600;6801:7100;12051:12350;19001:19300;28501:28800;   16881:17180;  ...
  21661:21960; 26001:26300 ];
listfrm = listfrm';
[fig] = GUI(TCave,listfrm,spk_inx,sessions);  
savefig([image_analysis_dest sessions '_GUI_TCave_w_spkEvent_stdbest']); 

%% Calculate SNR
% baseline standard deviation(BSD): 
TCave_stay = TCave(frm_stay,:);
BSD = std(TCave_stay, 0,1); %same number as # of cells

%get the peak values
pkvalues = cell(1,size(TCave,2));
SNR = zeros(1,size(TCave,2));
for c = 1:size(TCave,2)
    pkvalues{c} = TCave((spk_inx{c}),c);
    for v = 1:length(pkvalues{c})
        if pkvalues{c}(v)> TCave(spk_inx{c}(v)+1,c) % if the peak value is bigger than the next value in TCave
           continue
        else
           pkvalues{c}(v)= TCave(spk_inx{c}(v)+1,c);
        end
    end
    SNR(c) = mean(pkvalues{c}/BSD(c)); 
end
top1 = ceil(length(SNR)*0.8);           top2 = ceil(length(SNR)*0.9);
SNR_sort = sort(SNR,'descend');         
x1 = round(SNR_sort(top1),1);           x2 = round(SNR_sort(top2),1);
histSNR = figure;
histogram(SNR,'BinWidth',0.1); hold on;
vline(x1,'r');vline(x2,'r');
title(['histogram of signal-noise ratio' sessions]);
saveas(histSNR,[image_analysis_dest sessions '_SNR_histogram']); 


%%
neuron = 1:size(TCave,2);
threshold = 7.8;
badcells = neuron(SNR<threshold);
% delete bad cells in TCave, TCave_stay, mask image,
% spk_bi_cellmat,spk_inx, SNR, BSD,deriv,FRstay_cells,pkvalues

%figure; plot(TCave(:,c)); hold on;
%plot(spk_inx{c},pkvalues{c},'ro','MarkerSize', 6);hold on;
%xlim([2000 3000]);



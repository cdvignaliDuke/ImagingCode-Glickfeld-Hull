% the purpose of the script is to see if cells are more/less synchronized
% when behavior state changes. See paper "On-beam synchrony in the cerebellum as the mechanism for the timing and oordination of movement, 2009 PNAS"
% synchrony is basically calculated by doing pair-wise (of cells)
% correlation within each behavioral state, and if the correlation of these
% pairs increase/decrese, synchrony increases/decreases
%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042','200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
% behavior analysis results
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1','1064-200316_2'};

% behavior analysis results 
color_code = {'r','g','m','y','b'};

%% for each session, 1s running and 1s in stay

for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    dfOvF = dfOvF_output.dfOvF_btm_cl;
    
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{i} '\'];
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frms_airstim_stay = behav_output.frms_airstim_stay; % trials*frames, in each trial: 1-15:500ms before airpuff, 16:airpuff on, 17-31:after airpuff
    frms_airstim_stay = frms_airstim_stay'; %frames*trials
    frms_airstim_stay(end,:) = []; %delete the 31st frame in each trial so that time before airpuff and after airpuff is equal
    
    tlen = 15; % how long the time bin is
    airpuff = zeros(tlen,size(frms_airstim_stay,2));
    stay = zeros(tlen,size(frms_airstim_stay,2));
    % generate frame matrices for last n=tlen frames of running and n = tlen frames starting from 1s after running ends
    for n = 1: size(frms_airstim_stay,2)                                   
        airpuff(:,n) = frms_airstim_stay(tlen+1:end,n);     % frames*trials   
        stay(:,n) = frms_airstim_stay(1:tlen,n);
    end
    
    % generate dfOvF logic matrices    
    dfOvF_airpuff = zeros(size(airpuff,1),size(airpuff,2),size(dfOvF,2)); %frames*trials*cells
    dfOvF_stay = zeros(size(stay,1),size(stay,2),size(dfOvF,2)); 
    for c = 1:size(dfOvF,2)
        dfOvF_c = dfOvF(:,c);
        dfOvF_airpuff(:,:,c) = dfOvF_c(airpuff);
        dfOvF_stay(:,:,c) = dfOvF_c(stay);
    end
    
    %=================================================================================================================
    pairs = nchoosek((1:size(dfOvF,2)),2); %pairs of cells, pair*2
    %raw correlation: do correlation for each pair of cells in each trial. conditions:1s before running ends, 1-2s after running ends, dfOvF --> spk
    r_dfOvF_airpuff      = zeros(size(pairs,1),size(dfOvF_airpuff,2)); %pairs*trials
    pval_dfOvF_airpuff   = zeros(size(pairs,1),size(dfOvF_airpuff,2));
    r_dfOvF_stay    = zeros(size(pairs,1),size(dfOvF_airpuff,2));
    pval_dfOvF_stay = zeros(size(pairs,1),size(dfOvF_airpuff,2));
    
    for p = 1:size(pairs,1) %for each pair of cells
        for t = 1:size(dfOvF_airpuff,2) %for each trial
            [r_dfOvF_airpuff(p,t),pval_dfOvF_airpuff(p,t)] = corr(dfOvF_airpuff(:,t,pairs(p,1)),dfOvF_airpuff(:,t,pairs(p,2)));
            [r_dfOvF_stay(p,t),pval_dfOvF_stay(p,t)] = corr(dfOvF_stay(:,t,pairs(p,1)),dfOvF_stay(:,t,pairs(p,2)));
        end
    end
    
    r_dfOvF_airpuff_pairs = mean(r_dfOvF_airpuff,2);
    r_dfOvF_stay_pairs    = mean(r_dfOvF_stay,2);
    rawCorr = [r_dfOvF_stay_pairs,r_dfOvF_airpuff_pairs];
    
    rshiftdfOvF_airpuff  = zeros(size(pairs,1),size(dfOvF_airpuff,2)*(size(dfOvF_airpuff,2)-1));
    pshift_dfOvF_airpuff = zeros(size(pairs,1),size(dfOvF_airpuff,2)*(size(dfOvF_airpuff,2)-1));%pairs*trials combinations, combinations = ntrials*ntrials-1 (C ntrial,2)
    rshiftdfOvF_stay = zeros(size(pairs,1),size(dfOvF_airpuff,2)*(size(dfOvF_airpuff,2)-1));
    pshiftdfOvF_stay = zeros(size(pairs,1),size(dfOvF_airpuff,2)*(size(dfOvF_airpuff,2)-1));
    
    %shift predictor: this can take fucking long
    for p = 1:size(pairs,1)
        a = 0;
        p
        for t1 = 1:size(dfOvF_airpuff,2)
            for t2 = 1:size(dfOvF_airpuff,2)
                if t1~=t2 %for cellA trial1, corr with cellB trial2-N
                    a = a+1;
                    a
                    [rshiftdfOvF_airpuff(p,a),pshift_dfOvF_airpuff(p,a)] = corr(dfOvF_airpuff(:,t1,pairs(p,1)),dfOvF_airpuff(:,t2,pairs(p,2)));
                    [rshiftdfOvF_stay(p,a),pshiftdfOvF_stay(p,a)] = corr(dfOvF_stay(:,t1,pairs(p,1)),dfOvF_stay(:,t2,pairs(p,2)));
                end
            end
        end
    end
   
    rshiftdfOvF_airpuff_pairs    = mean(rshiftdfOvF_airpuff,2);
    rshiftdfOvF_stay_pairs       = mean(rshiftdfOvF_stay,2);
    
    eCorr_dfOvF_airpuff_pairs     = r_dfOvF_airpuff_pairs - rshiftdfOvF_airpuff_pairs;
    eCorr_dfOvF_stay_pairs    = r_dfOvF_stay_pairs - rshiftdfOvF_stay_pairs;
    eCorr = [eCorr_dfOvF_stay_pairs,eCorr_dfOvF_airpuff_pairs];
    
    x = [1,2]; x_plot = repmat(x,size(eCorr,1),1);
    eCorr_fig = figure;
    plot(x_plot',eCorr','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]);hold on;
    fast_errbar(x,eCorr,1,'color',[0.7922 0 0.1255]);
    xlabel ('behavioral state');xlim([0.5 2.5]); 
    axis square;
    ylim([-1 1]);
    set(gca,'XTick',x,'XTicklabel',{'stay','airpuff'});
    ylabel('correlation between pairs of cells');
    title(['correlation ' sessions{i}]);
    savefig([image_analysis_dest sessions{i}, '_synch_corr_tlen' num2str(tlen)]);
    
    save([image_analysis_dest sessions{i} '_dfOvF.mat'] ,'eCorr','tlen','-append');
end

%{
%% across sessions
eCorr_all = [];
for i = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    dfOvF_output = load([image_analysis_dest sessions{i} '_dfOvF.mat']);
    eCorr = dfOvF_output.eCorr;
    eCorr_all = cat(1,eCorr_all,eCorr);
end

x = [1,2]; x_plot = repmat(x,size(eCorr_all,1),1);
eCorr_fig = figure;
plot(x_plot',eCorr_all','.','LineStyle','none','MarkerSize',10,'color',[0.5294 0.5294 0.5294]);hold on;
fast_errbar(x,eCorr_all,1,'color',[0.7922 0 0.1255]);
xlabel ('behavioral state');xlim([0.5 2.5]);
axis square;
ylim([-1 1]);
set(gca,'XTick',x,'XTicklabel',{'run','stationary'});
ylabel('correlation between pairs of cells');
title('correlation across sessions');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_synch_corr.fig']);
%}


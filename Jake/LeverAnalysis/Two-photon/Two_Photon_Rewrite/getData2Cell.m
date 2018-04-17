function output = getData2Cell(dest_sub, output, id, day)

if exist(dest_sub)
    load([dest_sub '_cell_TCs.mat']);
    load([dest_sub '_F_risetime.mat']);
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub '_pvals.mat']);
    load([dest_sub '_cue_movies.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub '_cell_resp.mat']);
    
    if day == 1 || day == 5
        load([dest_sub 'same_cell.mat']);
        output.ncells{id} = length(same_cell_idx);
        allresp_cells_sc = 0*allresp_cells;
        allresp_cells_sc(same_cell_idx) = 1;
        allresp_cells = (allresp_cells & logical(allresp_cells_sc));
        
        OR_Cue_h2_sc = 0*OR_Cue_h2;
        OR_Cue_h2(same_cell_idx) = 1;
        OR_Cue_h2 = (OR_Cue_h2 & logical(OR_Cue_h2_sc));
        output.cue_cells{id} = cue_cells(same_cell_idx);
        output.rew_cells{id} = rew_cells(same_cell_idx);
        
        output.cue_cells_perc{id} = sum(output.cue_cells{id}) / length(same_cell_idx);
        output.rew_cells_perc{id} = sum(output.rew_cells{id}) / length(same_cell_idx);
        
        output.NR_F_onset{id} = (NR_riseIdx(same_cell_idx))' * double(ifi);
        output.NR_F_2_bout{id} = (NR_riseIdx(same_cell_idx) - mean(NR_lick_info.bout_onset))' * double(ifi);
        output.OR_resp_sc{id} = mean(OR_Cue_resp2(:,same_cell_idx)  - OR_Cue_base2(:,same_cell_idx),1);
        if day == 1
            output.NR_plot_ones{id} = ones(size(output.NR_F_onset{id}));
            
        else
            output.NR_plot_ones{id} = 2*ones(size(output.NR_F_onset{id}));
        end
    else
        output.ncells{id} = length(allresp_cells);
        output.cue_cells{id} = cue_cells;
        output.rew_cells{id} = rew_cells;
    end
    
    load([dest_sub 'parse_behavior.mat']);
    
    output.RS_cells{id} = allresp_cells;
    output.OR_resp_cells{id} = OR_Cue_h2;
    
   
    RS_cells{id} =  allresp_cells;
    output.tot_resp{id} = sum(RS_cells{id});
    
    output.cue_cells_perc{id} = sum(output.cue_cells{id}) / length(allresp_cells);
    output.rew_cells_perc{id} = sum(output.rew_cells{id}) / length(allresp_cells);
    
    % normal reward
    output.NR_TC{id} = avg_NR;
    output.NR_TC_RS{id} = avg_NR(RS_cells{id},:);
    
    
%     if day ~= 1
%         output.NR_TC_RewPos{id} = avg_NR(NR_Rew_resp_pos_cells,:);
%         output.NR_TC_RewNeg{id} = avg_NR(NR_Rew_resp_neg_cells,:);
%         output.NR_TC_ORcellNeg{id} = avg_NR(OR_Rew_resp_neg_cells,:);
%         
%         output.NR_TC_nolick_RewPos{id} = avg_NR_nolick(NR_Rew_resp_pos_cells,:);
%         output.NR_TC_nolick_RewNeg{id} = avg_NR_nolick(NR_Rew_resp_neg_cells,:);
%         output.NR_TC_nolick_ORcellNeg{id} = avg_NR_nolick(OR_Rew_resp_neg_cells,:);
%     end
    output.NR_TC_mean{id} = mean(avg_NR,1); %avg per animal
    output.NR_TC_sem{id} = std(avg_NR,1)./sqrt(size(avg_NR,1));
    output.NR_TC_RS_mean{id} = mean(avg_NR(RS_cells{id},:),1);
    output.NR_TC_RS_sem{id} = std(avg_NR(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
    
%     output.NR_resp_all{id} = mean(NR_Cue_resp - NR_Cue_base,1);
%     output.NR_resp_RS{id} = mean(NR_resp(:,NR_resp_cells)  - NR_base(:,NR_resp_cells),1);
   
    
    
%     output.NR_resp_all{id} = mean((NR_resp-NR_base),1);
%     output.NR_resp_mean{id} = mean(output.NR_resp_all{id},2);
%     output.NR_resp_sem{id} = std(output.NR_resp_all{id},[],2)./sqrt(output.ncells{id});
%     output.NR_resp_RS{id} = mean((NR_resp(:,RS_cells{id})-NR_base(:,RS_cells{id})),1);
%     output.NR_resp_RS_mean{id} = mean(output.NR_resp_RS{id},2);
%     output.NR_resp_RS_sem{id} = std(output.NR_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    
    %unexpected reward
    output.UR_TC{id} = avg_UR;
    output.UR_TC_RS{id} = avg_UR(RS_cells{id},:);
    
    output.UR_TC_mean{id} = mean(avg_UR,1); %avg per animal
    output.UR_TC_sem{id} = std(avg_UR,1)./sqrt(size(avg_UR,1));
    output.UR_TC_RS_mean{id} = mean(avg_UR(RS_cells{id},:),1);
    output.UR_TC_RS_sem{id} = std(avg_UR(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
    
%     output.UR_resp_all{id} = mean(UR_resp - UR_base,1);
    
%     output.UR_F_onset{id} = (UR_riseIdx - pre_cue_frames - 1)' * double(ifi);
    output.UR_F_2_bout{id} = (UR_riseIdx - mean(UR_lick_info.bout_onset))' * double(ifi);
%     output.UR_resp_all{id} = mean((UR_resp-UR_base),1);
%     output.UR_resp_mean{id} = mean(output.UR_resp_all{id},2);
%     output.UR_resp_sem{id} = std(output.UR_resp_all{id},[],2)./sqrt(output.ncells{id});
%     output.UR_resp_RS{id} = mean((UR_resp(:,RS_cells{id})-UR_base(:,RS_cells{id})),1);
%     output.UR_resp_RS_mean{id} = mean(output.UR_resp_RS{id},2);
%     output.UR_resp_RS_sem{id} = std(output.UR_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    
    %omitted reward
    output.OR_TC{id} = avg_OR;
    output.OR_TC_RS{id} = avg_OR(RS_cells{id},:);
    
    omitRewardIndx = logical(trial_outcome.omitRewardIndx);
    prevR_TC = []; prevNR_TC = []; trialLen = [];
    all_cue_resp = NR_Cue_resp - NR_Cue_base;
    all_cue_resp = [all_cue_resp; OR_Cue_resp - OR_Cue_base];
    
    all_movie_nolick = cat(1, NR_movie_nolick, OR_movie_nolick);
    omitTrialLeng = trial_outcome.trialLen(omitRewardIndx);
    trial_outcome.trialLen(omitRewardIndx) = [];
    trial_outcome.trialLen = [trial_outcome.trialLen omitTrialLeng];
    rj = [];
    for j = 1:size(all_movie_nolick,1)
        if  j > 1 && omitRewardIndx(j - 1 ) == 0
            prevR_TC = [prevR_TC; squeeze(all_movie_nolick(j, :, :))];
            trialLen = [trialLen trial_outcome.trialLen(j)];
        else
            prevNR_TC = [prevNR_TC; squeeze(all_movie_nolick(j, :, :))];
            rj = [rj;j];
        end
    end
    
    hasReward = ones(1, size(trial_outcome.omitRewardIndx,2));
    hasReward(logical(trial_outcome.omitRewardIndx)) = 0;
    omitRewardIndx = find(trial_outcome.omitRewardIndx);
    OR_prevR_TC = []; OR_prevNR_TC = [];
    for j = 1:size(OR_movie_nolick,1)
        if omitRewardIndx(j) - 1 > 0 && hasReward(omitRewardIndx(j) - 1) == 1
            OR_prevR_TC = [OR_prevR_TC; squeeze(OR_movie_nolick(j, :, :))];
           
        else
            OR_prevNR_TC = [OR_prevNR_TC; squeeze(OR_movie_nolick(j, :, :))];
           
        end
    end
    
    output.trialLen_prevR{id} = double(trialLen).*double(min((ifi)))/1000;
    all_cue_resp(rj,:) = [];
    output.all_resp_RS{id} = nanmean(all_cue_resp(:, RS_cells{id}),2); 
    output.all_TC_prevR{id} = prevR_TC;
    output.all_TC_prevNR{id} = prevNR_TC;
    output.OR_TC_prevR{id} = OR_prevR_TC;
    output.OR_TC_prevNR{id} = OR_prevNR_TC;
    
%     if day ~= 1
%         output.OR_TC_RewPos{id} = avg_OR(OR_Rew_resp_pos_cells,:);
%         output.OR_TC_RewNeg{id} = avg_OR(OR_Rew_resp_neg_cells,:);
%         
%         output.OR_TC_nolick_RewPos{id} = avg_OR_nolick(OR_Rew_resp_pos_cells,:);
%         output.OR_TC_nolick_RewNeg{id} = avg_OR_nolick(OR_Rew_resp_neg_cells,:);
%     end
    output.OR_TC_mean{id} = mean(avg_OR,1); %avg per animal
    output.OR_TC_sem{id} = std(avg_OR,1)./sqrt(size(avg_OR,1));
    output.OR_TC_RS_mean{id} = mean(avg_OR(RS_cells{id},:),1);
    output.OR_TC_RS_sem{id} = std(avg_OR(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
    
    output.OR_resp_all{id} = mean(OR_Cue_resp2 - OR_Cue_base2,1);
    
    
    
    
    
    output.OR_F_onset{id} = (OR_riseIdx)' * double(ifi);
    output.OR_F_2_bout{id} = (OR_riseIdx - mean(OR_lick_info.bout_onset))' * double(ifi);
    if day == 1
        output.OR_plot_ones{id} = ones(size(output.OR_F_onset{id}));
    else
        output.OR_plot_ones{id} = 2*ones(size(output.OR_F_onset{id}));
    end
    
    % frame info
    output.ifi{id} = ifi;
    output.pre_frames{id} = pre_cue_frames;
    output.post_frames{id} = post_cue_frames;
%     output.OR_resp_all{id} = mean((OR_resp-OR_base),1);
%     output.OR_resp_mean{id} = mean(output.OR_resp_all{id},2);
%     output.OR_resp_sem{id} = std(output.OR_resp_all{id},[],2)./sqrt(output.ncells{id});
%     output.OR_resp_RS{id} = mean((OR_resp(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
%     output.OR_resp_RS_mean{id} = mean(output.OR_resp_RS{id},2);
%     output.OR_resp_RS_sem{id} = std(output.OR_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
end
end
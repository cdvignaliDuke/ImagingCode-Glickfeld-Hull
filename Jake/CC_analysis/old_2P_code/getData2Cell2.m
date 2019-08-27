%function for loading variables for a given session and storing them in a cell array. 
%Rewritten to work without across days data on 9/7/18

function output = getData2Cell(dest_sub, output, id, day)
%check to make sure the directory exists
if exist(dest_sub)
    %load variables
    load([dest_sub '_cell_TCs.mat']);
    %load([dest_sub '_F_risetime.mat']);
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub '_pvals.mat']);
    load([dest_sub '_cue_movies.mat']);
    load([dest_sub '_cue_movies_lick.mat']);
    load([dest_sub '_cell_resp.mat']);
    load([dest_sub 'parse_behavior.mat']);
else
    return
end

%number of cells
output.ncells{id} = length(allresp_cells);
output.ncells_resp{id} = sum(allresp_cells);
RS_cells{id} =  allresp_cells;

%ALL responsive cells
output.RS_cells{id} = allresp_cells;
output.RS_cells_neg{id} = allresp_cells_neg;
output.RS_cells_pos{id} = allresp_cells_pos;
%fraction of all responsive cells
output.RS_cells{id} = allresp_cells/length(allresp_cells);
output.RS_cells_neg{id} = allresp_cells_neg/length(allresp_cells);
output.RS_cells_pos{id} = allresp_cells_pos/length(allresp_cells);

%CUE responsive cells
output.cue_cells{id} = cue_cells;
output.cue_cells_pos{id} = cue_cells_pos;
output.cue_cells_neg{id} = cue_cells_neg;
%fraction of cue responsive cells
output.cue_cells_perc{id} = cue_cells_perc;
output.cue_cells_perc_pos{id} = cue_cells_perc_pos;
output.cue_cells_perc_neg{id} = cue_cells_perc_neg;

%REWARD responsive cells
output.rew_cells{id} = rew_cells;
output.rew_cells_pos{id} = rew_cells_pos;
output.rew_cells_neg{id} = rew_cells_neg;
%fraction of reward responsive cells
output.rew_cells_perc{id} = rew_cells_perc;
output.rew_cells_perc_pos{id} = rew_cells_perc_pos;
output.rew_cells_perc_neg{id} = rew_cells_perc_neg;


%% normal reward

%NR cue responsive
output.NR_Cue_h_pos{id} = NR_Cue_h_pos;
output.NR_Cue_h_neg{id} = NR_Cue_h_neg;

%NR reward responsive
output.NR_Rew_h_pos{id} = NR_Rew_h_pos;
output.NR_Rew_h_neg{id} = NR_Rew_h_neg;

%% unexpected reward

%UR cue responsive
output.UR_Cue_h_pos{id} = UR_Cue_h_pos;
output.UR_Cue_h_neg{id} = UR_Cue_h_neg;

%UR reward responsive
output.UR_Rew_h_pos{id} = UR_Rew_h_pos;
output.UR_Rew_h_neg{id} = UR_Rew_h_neg;

%% omitted reward

%OR cue responsive
output.OR_Cue_h_pos{id} = OR_Cue_h_pos;
output.OR_Cue_h_neg{id} = OR_Cue_h_neg;
output.OR_Cue_h2_pos{id} = OR_Cue_h2_pos;
output.OR_Cue_h2_neg{id} = OR_Cue_h2_neg;

%OR reward responsive
output.OR_Rew_h_pos{id} = OR_Rew_h_pos;
output.OR_Rew_h_neg{id} = OR_Rew_h_neg;

%% frame info
output.ifi{id} = ifi;
output.pre_frames{id} = pre_cue_frames;
output.post_frames{id} = post_cue_frames;

%% old code Ziye wrote which I dont understand yet. 

% 
% omitRewardIndx = logical(trial_outcome.omitRewardIndx);
% prevR_TC = []; prevNR_TC = []; trialLen = [];
% all_cue_resp = NR_Cue_resp - NR_Cue_base;
% all_cue_resp = [all_cue_resp; OR_Cue_resp - OR_Cue_base];
% 
% all_movie_nolick = cat(1, NR_movie_nolick, OR_movie_nolick);
% omitTrialLeng = trial_outcome.trialLen(omitRewardIndx);
% trial_outcome.trialLen(omitRewardIndx) = [];
% trial_outcome.trialLen = [trial_outcome.trialLen omitTrialLeng];
% rj = [];
% for j = 1:size(all_movie_nolick,1)
%     if  j > 1 && omitRewardIndx(j - 1 ) == 0
%         prevR_TC = [prevR_TC; squeeze(all_movie_nolick(j, :, :))];
%         trialLen = [trialLen trial_outcome.trialLen(j)];
%     else
%         prevNR_TC = [prevNR_TC; squeeze(all_movie_nolick(j, :, :))];
%         rj = [rj;j];
%     end
% end
% 
% hasReward = ones(1, size(trial_outcome.omitRewardIndx,2));
% hasReward(logical(trial_outcome.omitRewardIndx)) = 0;
% omitRewardIndx = find(trial_outcome.omitRewardIndx);
% OR_prevR_TC = []; OR_prevNR_TC = [];
% for j = 1:size(OR_movie_nolick,1)
%     if omitRewardIndx(j) - 1 > 0 && hasReward(omitRewardIndx(j) - 1) == 1
%         OR_prevR_TC = [OR_prevR_TC; squeeze(OR_movie_nolick(j, :, :))];
%         
%     else
%         OR_prevNR_TC = [OR_prevNR_TC; squeeze(OR_movie_nolick(j, :, :))];
%         
%     end
% end
% 
% output.trialLen_prevR{id} = double(trialLen).*double(min((ifi)))/1000;
% all_cue_resp(rj,:) = [];
% output.all_resp_RS{id} = nanmean(all_cue_resp(:, RS_cells{id}),2);
% output.all_TC_prevR{id} = prevR_TC;
% output.all_TC_prevNR{id} = prevNR_TC;
% output.OR_TC_prevR{id} = OR_prevR_TC;
% output.OR_TC_prevNR{id} = OR_prevNR_TC;
% 
% %     if day ~= 1
% %         output.OR_TC_RewPos{id} = avg_OR(OR_Rew_resp_pos_cells,:);
% %         output.OR_TC_RewNeg{id} = avg_OR(OR_Rew_resp_neg_cells,:);
% %
% %         output.OR_TC_nolick_RewPos{id} = avg_OR_nolick(OR_Rew_resp_pos_cells,:);
% %         output.OR_TC_nolick_RewNeg{id} = avg_OR_nolick(OR_Rew_resp_neg_cells,:);
% %     end
% output.OR_TC_mean{id} = mean(avg_OR,1); %avg per animal
% output.OR_TC_sem{id} = std(avg_OR,1)./sqrt(size(avg_OR,1));
% output.OR_TC_RS_mean{id} = mean(avg_OR(RS_cells{id},:),1);
% output.OR_TC_RS_sem{id} = std(avg_OR(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
% 
% output.OR_resp_all{id} = mean(OR_Cue_resp2 - OR_Cue_base2,1);
% 
% 
% 
% 
% 
% output.OR_F_onset{id} = (OR_riseIdx)' * double(ifi);
% output.OR_F_2_bout{id} = (OR_riseIdx - mean(OR_lick_info.bout_onset))' * double(ifi);
% if day == 1
%     output.OR_plot_ones{id} = ones(size(output.OR_F_onset{id}));
% else
%     output.OR_plot_ones{id} = 2*ones(size(output.OR_F_onset{id}));
% end
% 
% 
% %     output.OR_resp_all{id} = mean((OR_resp-OR_base),1);
% %     output.OR_resp_mean{id} = mean(output.OR_resp_all{id},2);
% %     output.OR_resp_sem{id} = std(output.OR_resp_all{id},[],2)./sqrt(output.ncells{id});
% %     output.OR_resp_RS{id} = mean((OR_resp(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
% %     output.OR_resp_RS_mean{id} = mean(output.OR_resp_RS{id},2);
% %     output.OR_resp_RS_sem{id} = std(output.OR_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));


end



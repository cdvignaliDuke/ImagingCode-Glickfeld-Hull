%handles omission reward condition for lever task as well as licking
%averages. Correlation between licking and df/f magnitude. 

% second part of analysis for licking
% scatter for lick rate and dF/F
clear
file_info
out_base = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis\');
mouseID = mouseID(21:end);
date = date(21:end);
success_TC_notlick = []; success_TC_lick = []; licks_TC_lick = []; licks_TC_notlick = [];
sem_licks_TC_lick = []; sem_licks_TC_notlick = []; lickb_TC_lick = []; sem_lickb_TC_lick = [];
fail_TC_notlick = []; fail_TC_lick = [];
fail_lick_corr_sig = []; fail_lick_corr_notsig = [];fail_lick_corr_p = [];
success_lick_corr_sig = []; success_lick_corr_notsig = []; success_lick_corr_p = [];
success_lick_freq_all = []; successlong_lick_freq_all =[];
succ_hold_dur_all = []; fail_hold_dur_all = []; success_lick_all = []; fail_lick_all = [];
early_fail_lick_freq_all = []; late_fail_lick_freq_all =[]; early_fail_TC_RS = []; late_fail_TC_RS = [];
success_TC = []; success_TC_long = [];
faillong_lick_freq_all = []; fail_TC_long = [];
omitR_lick_freq_all = []; omitR_TC_long = []; omitRlong_lick_freq_all = []; omitR_TC = [];
omitR_prevR_lick_freq_all = []; omitR_prevNR_lick_freq_all = [];
itiRlong_lick_freq_all = []; itiR_TC_long = [];
omitR_TC_prevR_all = []; omitR_TC_prevNR_all = [];
for id =1:10%[2,4,6,11,9,10, 12, 13,15, 16] lick rate%[12 14 16]%[11 13 15]%1:size(mouseID,2)
    id
    
    for rID  = 1:2
        dest_sub  = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        %         dest_sub = ['Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure\', date{id}, '_', runID{rID}, '_', mouseID{id}, '\'];
        if exist(dest_sub)
            load([dest_sub '_cell_TCs.mat']);
            load([dest_sub '_cell_resp.mat']);
            %             load([dest_sub '_cell_resp_amp.mat']);
            load([dest_sub '_cell_categories.mat']);
            load([dest_sub '_release_movies.mat']);
            load([dest_sub 'parse_behavior.mat']);
            load([dest_sub '_spike_variance.mat']);
            load([dest_sub '_dff_lick.mat']);
            load([dest_sub '_lick_stats.mat']);
            load([dest_sub '_lick_movies.mat']);
            %             load([dest_sub '_omitR_previousTrial.mat'])
            
            ncells{id} = size(press_resp,2);
            
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
            
            lick_resp_cells{id} = lick_cells;
            lickandsuccess_cells{id} = lick_success_cells;
            licknotsuccess_cells{id} = lick_notsuccess_cells;
            successnotlick_cells{id} = success_notlick_cells;
            notlick_cells_all{id} = notlick_rel_cells;
            nrelease_cells{id} = length(release_resp_cells);
            
            success_lick_corr_all{id} = success_lick_corr;
            fail_lick_corr_all{id} = fail_lick_corr;
            itiR_lick_corr_all{id} = itiR_lick_corr;
            omitR_lick_corr_all{id} = omitR_lick_corr;
            
            tot_resp{id} = length(RS_cells{id});
            tot_success_cell{id} = length(success_resp_cells);
            tot_early_cell{id} = length(fail_resp_cells);
            tot_press_cell{id} = length(press_resp_cells);
            success_cell_only{id} = length(success_only_cells);
            success_fail_cell{id} = length(intersect(success_resp_cells, fail_resp_cells));
            notlick_cell_tot{id} = length(notlick_cells);
            
            success_lick_corr_sig = [success_lick_corr_sig; success_lick_corr(success_lick_corr(:,2) < 0.05, 1)];
            success_lick_corr_notsig = [success_lick_corr_notsig; success_lick_corr(success_lick_corr(:,2) >= 0.05, 1)];
            success_lick_corr_p = [success_lick_corr_p; success_lick_corr(:,2)];
            
            fail_lick_corr_sig = [fail_lick_corr_sig; fail_lick_corr(fail_lick_corr(:,2) < 0.05, 1)];
            fail_lick_corr_notsig = [fail_lick_corr_notsig; fail_lick_corr(fail_lick_corr(:,2) >= 0.05, 1)];
            fail_lick_corr_p = [fail_lick_corr_p; fail_lick_corr(:,2)];
            %             TC_length{id} = size(avg_release,2);
            %             TC_ifi{id} = ifi;
            %             pre_frames{id} = pre_release_frames;
            %             post_frames{id} = post_release_frames;
            %             pre_framesL{id} = pre_release_frames2;
            %             post_framesL{id} = post_release_frames2;
            %
            
            %             success_TC{id} = avg_success;
            %             success_TC_RS{id} = avg_success(RS_cells{id},:);
            %             success_TC_RS_L{id} = avg_success_long(RS_cells{id},:);
            %             success_TC_RL{id} = avg_success(RL_cells{id},:);
            %             success_TC_trans{id} = avg_success(trans_cell_ind,:);
            %             success_TC_sus{id} = avg_success(sustain_cell_ind,:);
            success_TC_notlick = [success_TC_notlick; avg_success(notlick_cells_all{id},:)];
            success_TC_lick = [success_TC_lick; avg_success(lick_resp_cells{id}, :)];
            [success_lick_freq, ~] = binLicks(success_lick, double(ifi));
            success_lick_freq_all = [success_lick_freq_all; success_lick_freq];
            success_TC = [success_TC; avg_success(RS_cells{id},:)];
            
            [successlong_lick_freq, sem_successlong_lick] = binLicks(success_lick_long, double(ifi));
            successlong_lick_freq_all = [successlong_lick_freq_all; successlong_lick_freq];
            success_TC_long = [success_TC_long; avg_success_long(RS_cells{id},:)];
            
            success_lick_all = [success_lick_all mean(success_lick(:,pre_release_frames:pre_release_frames+500/double(ifi)),2)'*1000/double(ifi)];
            fail_lick_all = [fail_lick_all mean(fail_lick(:,pre_release_frames:pre_release_frames+500/double(ifi)),2)'*1000/double(ifi)];
            succ_hold_dur_all = [succ_hold_dur_all trial_outcome.succ_hold_dur/1000];
            fail_hold_dur_all = [fail_hold_dur_all trial_outcome.fail_hold_dur/1000];
            
            if ~isempty(notlick_cells_all{id})
                success_resp_RS_notlick{id} = mean((success_resp(:,notlick_cells_all{id})-success_base(:,notlick_cells_all{id})),1);
                fail_resp_RS_notlick{id} = mean((fail_resp(:,notlick_cells_all{id})-fail_base(:,notlick_cells_all{id})),1);
            else
                success_resp_RS_notlick{id} = nan;
                fail_resp_RS_notlick{id} = nan;
            end
            %             success_TC_mean{id} = mean(avg_success,1);
            %             success_TC_sem{id} = std(avg_success,1)./sqrt(size(avg_success,1));
            %             success_TC_RS_mean{id} = mean(avg_success(RS_cells{id},:),1);
            %             success_TC_RS_sem{id} = std(avg_success(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            %             success_TC_RL_mean{id} = mean(avg_success(release_resp_cells,:),1);
            %             success_TC_RL_sem{id} = std(avg_success(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            %             success_TC_lick{id} = avg_success(lick_success_cells,:);
            %             success_TC_nolick{id} = avg_success(success_notlick_cells,:);
            %             nosuccess_TC_lick{id} = avg_success(lick_notsuccess_cells,:);
            % %             lick_TC{id} = avg_lick(lick_resp_cells,:);
            %
            %             success_resp_all{id} = mean((success_resp-success_base),1);
            %             success_resp_all_2{id} = mean((success_resp-success_base),2)';
            %             success_resp_mean{id} = mean(success_resp_all{id},2);
            %             success_resp_sem{id} = std(success_resp_all{id},[],2)./sqrt(ncells{id});
            success_OR_resp_RS_lateWin{id} = mean((success_resp2(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            success_OR_resp_RS_earlyWin{id} = mean((success_resp1(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            %             success_resp_RS_mean{id} = mean(success_resp_RS{id},2);
            %             success_resp_RS_sem{id} = std(success_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            %             success_resp_RL{id} = mean((success_resp(:,release_resp_cells)-success_base(:,release_resp_cells)),1);
            %             success_resp_RL_mean{id} = mean(success_resp_RL{id},2);
            %             success_resp_RL_sem{id} = std(success_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            %             success_resp_RS_succonly{id} = mean((success_resp(:,success_only_cells)-success_base(:,success_only_cells)),1);
            %             success_resp_RS_succfail{id} = mean((success_resp(:,fail_and_success_cells)-success_base(:,fail_and_success_cells)),1);
            %             success_resp_SC{id} = mean((success_resp(:,success_resp_cells)-success_base(:,success_resp_cells)),1);
            %             success_resp_SC_mean{id} = mean(success_resp_SC{id},2);
            %             success_resp_SC_sem{id} = std(success_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            %             success_resp_RL_trans{id} = mean((success_resp(:,trans_cell_ind)-success_base(:,trans_cell_ind)),1);
            %             success_resp_RL_sus{id} = mean((success_resp(:,sustain_cell_ind)-success_base(:,sustain_cell_ind)),1);
            %             success_resp_RS_2{id} = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),2)';
            %
            %             %             success_resp_peak_mean{id} = mean(mean(success_peak));
            %
            %             tooFast_TC{id} = avg_tooFast;
            %             tooFast_TC_RS{id} = avg_tooFast(RS_cells{id},:);
            %             tooFast_TC_RS_L{id} = avg_tooFast_long(RS_cells{id},:);
            %             tooFast_TC_RL{id} = avg_tooFast(RL_cells{id},:);
            %             tooFast_resp_all{id} = mean((tooFast_resp-tooFast_base),1);
            %             tooFast_resp_RS{id} = mean((tooFast_resp(:,RS_cells{id})-tooFast_base(:,RS_cells{id})),1);
            %             tooFast_resp_RL{id} = mean((tooFast_resp(:,release_resp_cells)-tooFast_base(:,release_resp_cells)),1);
            %
            %             fail_TC{id} = avg_fail;
            %             fail_TC_RS{id} = avg_fail(RS_cells{id},:);
            %             fail_TC_RS_L{id} = avg_fail_long(RS_cells{id},:);
            %             fail_TC_RL{id} = avg_fail(RL_cells{id},:);
            %             fail_TC_trans{id} = avg_fail(trans_cell_ind,:);
            %             fail_TC_sus{id} = avg_fail(sustain_cell_ind,:);
            %             fail_TC_mean{id} = mean(avg_fail,1);
            fail_TC_notlick = [fail_TC_notlick; avg_fail(notlick_cells_all{id},:)];
            fail_TC_lick = [fail_TC_lick; avg_fail(lick_resp_cells{id}, :)];
            [fail_lick_freq, ~] = binLicks(fail_lick_long, double(ifi));
            faillong_lick_freq_all = [faillong_lick_freq_all; fail_lick_freq];
            fail_TC_long = [fail_TC_long; avg_fail_long(RS_cells{id},:)];
            %             fail_TC_sem{id} = std(avg_fail,1)./sqrt(size(avg_fail,1));
            %             fail_TC_RS_mean{id} = mean(avg_fail(RS_cells{id},:),1);
            %             fail_TC_RS_sem{id} = std(avg_fail(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            %             fail_TC_RL_mean{id} = mean(avg_fail(release_resp_cells,:),1);
            %             fail_TC_RL_sem{id} = std(avg_fail(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            %             fail_resp_all{id} = mean((fail_resp-fail_base),1);
            %             fail_resp_all_2{id} = mean((fail_resp-fail_base),2)';
            %             fail_resp_mean{id} = mean(fail_resp_all{id},2);
            %             fail_resp_sem{id} = std(fail_resp_all{id},[],2)./sqrt(ncells{id});
            %             fail_resp_RS{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),1);
            %             fail_resp_RS_mean{id} = mean(fail_resp_RS{id},2);
            %             fail_resp_RS_sem{id} = std(fail_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            %             fail_resp_RL{id} = mean((fail_resp(:,release_resp_cells)-fail_base(:,release_resp_cells)),1);
            %             fail_resp_RL_mean{id} = mean(fail_resp_RL{id},2);
            %             fail_resp_RL_sem{id} = std(fail_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            %             fail_resp_SC{id} = mean((fail_resp(:,success_resp_cells)-fail_base(:,success_resp_cells)),1);
            %             fail_resp_SC_mean{id} = mean(fail_resp_SC{id},2);
            %             fail_resp_SC_sem{id} = std(fail_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            %             fail_resp_RL_trans{id} = mean((fail_resp(:,trans_cell_ind)-fail_base(:,trans_cell_ind)),1);
            %             fail_resp_RL_sus{id} = mean((fail_resp(:,sustain_cell_ind)-fail_base(:,sustain_cell_ind)),1);
            %             fail_resp_RS_2{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),2)';
            %
            %               [omitR_lick_freq, ~] = binLicks(omitR_lick, double(ifi));
            %               omitR_lick_freq_all = [omitR_lick_freq_all; omitR_lick_freq];
            if ~isempty(avg_OR)
                omitR_TC_RS = [omitR_TC; avg_OR(RS_cells{id},:)];
                
                [omitR_lick_freq, ~] = binLicks(omitR_lick_long, double(ifi));
                omitRlong_lick_freq_all = [omitRlong_lick_freq_all; omitR_lick_freq];
                omitR_TC_long = [omitR_TC_long; avg_OR_long(RS_cells{id},:)];
                omitR_resp_RS_lateWin{id} = mean((OR_resp2(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
                omitR_resp_RS_earlyWin{id} = mean((OR_resp1(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
            end
            %               avg_omitR_prevR_TC = squeeze(mean(omitR_prevR_TC,3));
            %               avg_omitR_prevNR_TC = squeeze(mean(omitR_prevNR_TC, 3));
            %               omitR_TC_prevR_all = [omitR_TC_prevR_all; omitR_prevR_TC(RS_cells{id},:)];
            %               omitR_TC_prevNR_all = [omitR_TC_prevNR_all; omitR_prevNR_TC(RS_cells{id},:)];
            %
            %               [omitR_prevR_lick_freq, ~] = binLicks(omitR_prevR_lick, double(ifi));
            %               omitR_prevR_lick_freq_all = [omitR_prevR_lick_freq_all; omitR_prevR_lick_freq];
            %
            %               [omitR_prevNR_lick_freq, ~] = binLicks(omitR_prevNR_lick, double(ifi));
            %               omitR_prevNR_lick_freq_all = [omitR_prevNR_lick_freq_all; omitR_prevNR_lick_freq];
            %               [itiR_lick_freq, ~] = binLicks(itiR_lick_long, double(ifi));
            %               itiRlong_lick_freq_all = [itiRlong_lick_freq_all; itiR_lick_freq];
            %               itiR_TC_long = [itiR_TC_long; avg_IR_long(RS_cells{id},:)];
            
            early_fail_TC_RS = [early_fail_TC_RS; avg_early_fail(RS_cells{id},:)];
            late_fail_TC_RS = [late_fail_TC_RS; avg_late_fail(RS_cells{id},:)];
            
            [early_fail_lick_freq, ~] = binLicks(early_fail_lick, double(ifi));
            early_fail_lick_freq_all = [early_fail_lick_freq_all; early_fail_lick_freq];
            
            [late_fail_lick_freq, ~] = binLicks(late_fail_lick, double(ifi));
            late_fail_lick_freq_all = [late_fail_lick_freq_all; late_fail_lick_freq];
            
            %               nlickF = 3;
            %               tlick = ((-pre_release_frames:nlickF:post_release_frames + 1.5).*double(100/3));
            %                   figure;
            %                   bar(tlick(1:end-1), late_fail_lick_freq, 'k');
            %                   hold on;
            %                   bar(tlick(1:end-1), early_fail_lick_freq, 'b');
            
            
            licks_TC_lick = [licks_TC_lick; avg_single_lick(lick_resp_cells{id},:)];
            sem_licks_TC_lick = [sem_licks_TC_lick; sem_single_lick(lick_resp_cells{id},:)];
            lickb_TC_lick = [lickb_TC_lick; avg_lick(lick_resp_cells{id},:)];
            sem_lickb_TC_lick = [sem_lickb_TC_lick; sem_lick(lick_resp_cells{id},:)];
            licks_TC_notlick = [licks_TC_notlick; avg_single_lick(notlick_cells_all{id},:)];
            sem_licks_TC_notlick = [sem_licks_TC_notlick; sem_single_lick(notlick_cells_all{id},:)];
            
            licks_movie{id} = single_lick_movie(:,lick_resp_cells{id},:,:);
            licks_movie_notlick{id} = single_lick_movie(:,notlick_cells_all{id},:,:);
            
            
            %             lickb_resp_all{id} = mean((lickb_resp-lickb_base),1);
            %             lickb_resp_mean{id} = mean(lickb_resp_all{id},2);
            %             lickb_resp_sem{id} = std(lickb_resp_all{id},[],2)./sqrt(ncells{id});
            %             lickb_resp_RS{id} = mean((lickb_resp(:,RS_cells{id})-lickb_base(:,RS_cells{id})),1);
            %             lickb_resp_RS_mean{id} = mean(lickb_resp_RS{id},2);
            %             lickb_resp_RS_sem{id} = std(lickb_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            lickb_resp_lick{id} = mean((lickb_resp(:,lick_resp_cells{id})-lickb_base(:,lick_resp_cells{id})),1);
            %             lickb_resp_LK_mean{id} = mean(lickb_resp_LK{id},2);
            %             lickb_resp_LK_sem{id} = std(lickb_resp_LK{id},[],2)./sqrt(length(lick_resp_cells{id}));
            %
            %             licks_resp_all{id} = mean((licks_resp-licks_base),1);
            %             licks_resp_mean{id} = mean(licks_resp_all{id},2);
            %             licks_resp_sem{id} = std(licks_resp_all{id},[],2)./sqrt(ncells{id});
            %             licks_resp_RS{id} = mean((licks_resp(:,RS_cells{id})-licks_base(:,RS_cells{id})),1);
            %             licks_resp_RS_mean{id} = mean(licks_resp_RS{id},2);
            %             licks_resp_RS_sem{id} = std(licks_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            licks_resp_lick{id} = mean((licks_resp(:,lick_resp_cells{id})-licks_base(:,lick_resp_cells{id})),1);
            %             licks_resp_LK_mean{id} = mean(licks_resp_LK{id},2);
            %             licks_resp_LK_sem{id} = std(licks_resp_LK{id},[],2)./sqrt(length(lick_resp_cells{id}));
            %
            %             lick_rate_success{id} = lick_data.lickRateS;
            %             lick_rate_fail{id} = lick_data.lickRateF(lick_data.lickRateF<20);
            %             temp_fail_resp_RS_2 = fail_resp_RS_2{id};
            %             fail_resp_RS_2{id} = temp_fail_resp_RS_2(lick_data.lickRateF<20);
            %
            %             s_resp_RS_2_wLick{id} = success_resp_RS_2{id}(lick_data.lickRateS > 0);
            %             s_resp_RS_2_wNoLick{id} = success_resp_RS_2{id}(lick_data.lickRateS == 0);
            %
            %             f_resp_RS_2_wLick{id} = fail_resp_RS_2{id}(lick_data.lickRateF > 0 & lick_data.lickRateF<20);
            %             f_resp_RS_2_wNoLick{id} = fail_resp_RS_2{id}(lick_data.lickRateF == 0);
            %
            % %             lick_Rate_success_2hz = lick_data.lickRateS(lick_data.lickRateS==2);
            % %             success_resp_rate_2hz = success_resp_RS_2{id};
            % %             success_resp_rate_2hz = success_resp_rate_2hz(lick_data.lickRateS==2);
            %
            %             mean_lickRate_succ{id} = mean(lick_data.lickRateS)*ones(size(RS_cells{id}));
            %             mean_lickRate_fail{id} = mean(lick_data.lickRateF)*ones(size(RS_cells{id}));
            %
            %             succ_temp = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),2);
            %             succ_temp(lick_data.lickRateS == 0) = [];
            %             success_resp_withLick{id} = succ_temp;
            %
            %             fail_temp = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),2);
            %             fail_temp(lick_data.lickRateF == 0) = [];
            %             fail_resp_withLick{id} = fail_temp;
            %
            %             early_percent{id} = sum(lick_data.lickRateF > 0)/length(lick_data.lickRateF)*100;
            %             success_resp_mean{id} = mean(success_resp_RS_2{id});
            %             fail_resp_mean{id} = mean(fail_resp_RS_2{id});
            %
            %             % to match lick rate for plotting resp amplitude
            %
            % %             success_resp_LK{id} = mean((success_resp(:,lick_resp_cells{id})-success_base(:,lick_resp_cells{id})),1);
            % %             success_resp_LK_mean{id} = mean(success_resp_LK{id},2);
            % %             success_resp_LK_sem{id} = std(success_resp_LK{id},[],2)./sqrt(length(lick_resp_cells{id}));
            % %
            % %             fail_resp_LK{id} = mean((fail_resp(:,lick_resp_cells{id})-fail_base(:,lick_resp_cells{id})),1);
            % %             fail_resp_LK_mean{id} = mean(fail_resp_LK{id},2);
            % %             fail_resp_LK_sem{id} = std(fail_resp_LK{id},[],2)./sqrt(length(lick_resp_cells{id}));
            %
            %
            %
            %             succ_hold_dur{id} = trial_outcome.succ_hold_dur;
            %             fail_hold_dur{id} = trial_outcome.fail_hold_dur;
            %
            %
            %             peak2release_1{id} = peak2release_std1;
            %             peak2release_2{id} = peak2release_std2;
            %             peak2cue_1{id} = peak2cue_std1;
            %             peak2cue_2{id} = peak2cue_std2;
            
        end
    end
end

total_cells = sum(cell2mat(ncells));
total_resp = sum(cell2mat(tot_resp));
total_release = sum(cell2mat(nrelease_cells));
total_lick_cell = length(cell2mat(lick_resp_cells));
total_notlick_cell = length(cell2mat(notlick_cells_all));

tt = (-pre_release_frames:post_release_frames)*double(ifi);
tt2 = (-pre_release_frames2:post_release_frames2)*double(ifi);

nlickF = 3;
tlick = ((-pre_release_frames:nlickF:post_release_frames + 1.5).*double(100/3));
tlick2 = ((-pre_release_frames2:nlickF:post_release_frames2 + 1.5).*double(100/3));

col_mat_s = [0.5 0.5 0.5];

% scatter plot of dff for success and omit
fig = figure;
[h_os, p_os, stats_os] = scatter_plot(mouseID(1:15), success_OR_resp_RS_lateWin, omitR_resp_RS_lateWin, col_mat_s);
hold on
xlim([-0.05 0.4]); ylim([-0.05 0.4])
x = -0.05:0.1:0.4; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlabel('success dF/F')
ylabel('omission dF/F')
title(['p value = ', num2str(p_os)]);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_successOmitR_lateWindow_amp_scatter.fig']);
print([out_base 'Summary_successOmitR_lateWindow_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_successOmitR_lateWindow_amp_scatter.pdf'], '-dpdf');

% scatter plot of dff for success and omit using early window
fig = figure;
[h_os2, p_os2, stats_os2] = scatter_plot(mouseID(1:15), success_OR_resp_RS_earlyWin, omitR_resp_RS_earlyWin, col_mat_s);
hold on
x = -0.05:0.1:0.5; y=x;
plot(x,y,'-k')
xlim([-0.05 0.4]); ylim([-0.05 0.4])
vline(0,'--k')
hline(0,'--k')
xlabel('success dF/F')
ylabel('omission dF/F')
title(['p value = ', num2str(p_os2), ', early Window']);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_successOmitR_earlyWindow_amp_scatter_100165.fig']);
print([out_base 'Summary_successOmitR_earlyWindow_amp_scatter_100165.eps'], '-depsc');
print([out_base 'Summary_successOmitR_earlyWindow_amp_scatter_100165.pdf'], '-dpdf');

% fig = figure;
% hold on;
% shadedErrorBar(tt, nanmean(early_fail_TC_RS,1)-0.01, nanstd(early_fail_TC_RS,[],1)./sqrt(size(early_fail_TC_RS,1)), 'm');
% shadedErrorBar(tt, nanmean(late_fail_TC_RS,1)+0.01, nanstd(late_fail_TC_RS,[],1)./sqrt(size(late_fail_TC_RS,1)), 'r');
% title('Early early TC- magenta and late early TC- red');
% xlabel('Time (s)')
% ylabel('dF/F')
% xlim([-0.5 1.5])

fig = figure;
hold on;
shadedErrorBar(tt/1000, nanmean(omitR_TC_RS,1)+0.002, nanstd(omitR_TC_RS,[],1)./sqrt(size(omitR_TC_RS,1)), 'm');
shadedErrorBar(tt/1000, nanmean(late_fail_TC_RS,1)+0.018, nanstd(late_fail_TC_RS,[],1)./sqrt(size(late_fail_TC_RS,1)), 'r');
title('omitR TC- magenta and late early TC- red');
xlabel('Time (s)')
ylabel('dF/F')
xlim([-0.5 1])
saveas(fig, [out_base 'Summary_avgLateEarlyOmit_TCs.fig']);
print([out_base 'Summary_avgLateEarlyOmit_TCs.eps'], '-depsc');
print([out_base 'Summary_avgLateEarlyOmit_TCs.pdf'], '-dpdf');

% plot lick rate for early earlies and late earlies
fig = figure;
hold on
bar(tlick(1:end-1), nanmean(early_fail_lick_freq_all), 'b');
errorbar(tlick(1:end-1), nanmean(early_fail_lick_freq_all), nanstd(early_fail_lick_freq_all)/sqrt(size(early_fail_lick_freq_all,1)), '.r');
bar(tlick(1:end-1), nanmean(late_fail_lick_freq_all), 'k');
errorbar(tlick(1:end-1), nanmean(late_fail_lick_freq_all), nanstd(late_fail_lick_freq_all)/sqrt(size(late_fail_lick_freq_all,1)), '.r');
xlabel('time(ms)'); ylabel('lick rate(Hz)')
title('lick rate for early earlies- blue and late earlies- black');
saveas(fig, [out_base 'Summary_lickRate_earlylateEarly.fig']);
print([out_base 'Summary_lickRate_earlylateEarly.eps'], '-depsc');
print([out_base 'Summary_lickRate_earlylateEarly.pdf'], '-dpdf');

% plot binned lick rate vs hold time for early and correct
sResp = success_lick_all;
sHold = succ_hold_dur_all;
sResp(isnan(sHold)) = [];
sHold(isnan(sHold)) = [];
[LinearCoeffS1, fitS1] = polyfit(sHold, sResp, 1);
CorrSfit1 = polyval(LinearCoeffS1, sHold);
[BS,idx] = histc(sHold,0:0.25:round(max(sHold)));
sHold_bins = accumarray(idx(:),sHold,[],@mean);
sResp_std = accumarray(idx(:),sResp,[],@sem);
sResp_bins = accumarray(idx(:),sResp,[],@mean);
[LinearCoeffS2, fitS2] = polyfit(sHold_bins(BS>=10), sResp_bins(BS>=10), 1);
CorrSfit2 = polyval(LinearCoeffS2, sHold_bins(BS>=10));

fResp = fail_lick_all;
fResp(isnan(fResp)) = [];
fHold = fail_hold_dur_all;
fHold(isnan(fHold)) = [];
[LinearCoeffF1, fitF1] = polyfit(fHold, fResp, 1);
CorrFfit1 = polyval(LinearCoeffF1, fHold);
[BF,idx] = histc(fHold,0:0.25:round(max(fHold)));
fHold_bins = accumarray(idx(:),fHold,[],@mean);
fResp_std = accumarray(idx(:),fResp,[],@sem);
fResp_bins = accumarray(idx(:),fResp,[],@mean);
[LinearCoeffF2, fitF2] = polyfit(fHold_bins(BF>=10), fResp_bins(BF>=10), 1);
CorrFfit2 = polyval(LinearCoeffF2, fHold_bins(BF>=10));

fig = figure;
% scatter_plot(mouseID, succ_hold_dur, success_resp_RS_2, col_mat_s);
hold on
errorbar(sHold_bins(BS>=10), sResp_bins(BS>=10), sResp_std(BS>=10), '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
% plot(sHold_bins(BS>=10), CorrSfit2, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);

errorbar(fHold_bins(BF>=10), fResp_bins(BF>=10), fResp_std(BF>=10), '.', 'MarkerSize', 10, 'color', [0.9,0,0])
% plot(fHold_bins(BF>=10), CorrFfit2, 'color', [0.9, 0,0], 'Linewidth', 1.5);

%plot(x,y,'-k')
title('licking rate black- correct, red- early');
xlabel('hold time(s)');
ylabel('mean lick rate (Hz)');
saveas(fig, [out_base 'Summary_lickRate_holdtime_trial.fig']);
print([out_base 'Summary_lickRate_holdtime_trial.eps'], '-depsc');
print([out_base 'Summary_lickRate_holdtime_trial.pdf'], '-dpdf');


% plot time course of single lick and lick bout
fig = figure;
hold on
shadedErrorBar(tt, mean(licks_TC_lick,1), std(licks_TC_lick,[],1)./sqrt(size(licks_TC_lick,1)), 'b');
shadedErrorBar(tt, mean(lickb_TC_lick,1)+0.002, std(lickb_TC_lick,[],1)./sqrt(size(lickb_TC_lick,1)), 'g');
xlabel('Time (ms)')
ylabel('dff')
% xlim([-600 1000]); ylim([0 10])
title('Average time course for single lick- blue and lick bout- green for lick resp cells');
saveas(fig, [out_base 'Summary_TC_licks.fig']);
print([out_base 'Summary_TC_licks.eps'], '-depsc');
print([out_base 'Summary_TC_licks.pdf'], '-dpdf');

% scatter plot of dff for single lick and lick bout
fig = figure;
[h_licksb, p_licksb, stats_licksb] = scatter_plot(mouseID, licks_resp_lick, lickb_resp_lick, col_mat_s);
hold on
x = -0.2:0.1:0.4; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlim([-0.2 0.4]); ylim([-0.2 0.4])
xlabel('single lick dF/F')
ylabel('lick bout dF/F')
title(['p value = ', num2str(p_licksb)]);
% title(['lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_lick_amp_scatter.fig']);
print([out_base 'Summary_lick_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_lick_amp_scatter.pdf'], '-dpdf');



fig = figure;
[h_sfRS, p_sfRS, stats_sfRS] = scatter_plot(mouseID, success_resp_RS_notlick, fail_resp_RS_notlick, col_mat_s);
hold on
x = 0:0.1:0.2; y=x;
plot(x,y,'-k')
vline(0,'--k')
hline(0,'--k')
xlabel('Correct dF/F')
ylabel('Early dF/F')
title(['Not lick but release resp cells- n = ' num2str(sum(cell2mat(notlick_cell_tot))) ' out of total resp cells n = ' num2str(total_resp)])
saveas(fig, [out_base 'Summary_SuccessFail_LickNoLick_scatter.fig']);
print([out_base 'Summary_SuccessFail_LickNoLick_scatter.eps'], '-depsc');
print([out_base 'Summary_SuccessFail_LickNoLick_scatter.pdf'], '-dpdf');


fig = figure;
hold on;
shadedErrorBar(tt, mean(omitR_TC_prevNR_all,1)+0.05, std(omitR_TC_prevNR_all,[],1)./sqrt(size(omitR_TC_prevNR_all,1)), 'r');
shadedErrorBar(tt, mean(omitR_TC_prevR_all,1), std(omitR_TC_prevR_all,[],1)./sqrt(size(omitR_TC_prevR_all,1)), 'b');
bar(tlick(1:end-1), mean(omitR_prevR_lick_freq_all)/500, 'b');
bar(tlick(1:end-1), mean(omitR_prevNR_lick_freq_all)/500, 'r');
xlabel('Time (ms)');xlim([-500 500])
ylabel('dF/F')
title('omit before Reward- blue, omit before no reward- red');
saveas(fig, [out_base 'Summary_TC_omitR_prevR.fig']);
print([out_base 'Summary_TC_omitR_prevR.eps'], '-depsc');
print([out_base 'Summary_TC_omitR_prevR.pdf'], '-dpdf');

%figure for previous reward/no reward
fStart = round(250/double(ifi)) + pre_release_frames + 1;
fEnd = round(350/double(ifi)) + pre_release_frames + 1;
omitR_dff_prevNR = omitR_TC_prevNR_all(:, fStart:fEnd)+0.02;
omitR_dff_prevR = omitR_TC_prevR_all(:, fStart:fEnd);

fig = figure;
scatter(mean(omitR_dff_prevNR,2), mean(omitR_dff_prevR,2),4, 'MarkerEdgeColor', [0.5 0.5 0.5])
x = -0.3:0.7;
y = x;
hold on; plot(x,y, '-k')
h = errorbarxy(mean(mean(omitR_dff_prevNR)), mean(mean(omitR_dff_prevR)),  std(mean(omitR_dff_prevNR,2))./sqrt(size(omitR_dff_prevNR,1)), std(mean(omitR_dff_prevR,2))./sqrt(size(omitR_dff_prevR,1)),{'o', [0.9 0 0], [0.9 0 0], [0.9 0 0]});
set(h.hMain,'LineWidth', 1);
[h,p,ci,stats] = ttest(mean(omitR_dff_prevNR,2),  mean(omitR_dff_prevR,2));
xlabel('omit after No reward dF/F')
ylabel('omit after reward dF/F')
title(['nCells = ', num2str(size(omitR_dff_prevNR,1)), ', p = ', num2str(p)]);
saveas(fig, [out_base 'Summary_omitR_prevR_scatter.fig']);
print([out_base 'Summary_omitR_prevR_scatter.eps'], '-depsc');
print([out_base 'Summary_omitR_prevR_scatter.pdf'], '-dpdf');

fig = figure; p = panel(); p.pack(2,2);
p(1,1).select()
hold on;xlim([-600 1000]);
shadedErrorBar(tt2, mean(success_TC_long,1), std(success_TC_long,[],1)./sqrt(size(success_TC_long,1)), 'k');
shadedErrorBar(tt2, mean(omitR_TC_long,1), std(omitR_TC_long,[],1)./sqrt(size(omitR_TC_long,1)), 'r');
ylabel('dF/F')
title(['Responsive cells- n = ' num2str(total_resp)])
p(2,1).select()
hold on
bar(tlick2(1:end-1), mean(successlong_lick_freq_all), 'k'); errorbar(tlick2(1:end-1), mean(successlong_lick_freq_all), std(successlong_lick_freq_all,[],1)/sqrt(size(successlong_lick_freq_all,1)), '.k');
bar(tlick2(1:end-1), mean(omitRlong_lick_freq_all), 'r'); errorbar(tlick2(1:end-1), mean(omitRlong_lick_freq_all), std(omitRlong_lick_freq_all,[],1)/sqrt(size(omitRlong_lick_freq_all,1)), '.r');
xlabel('Time (ms)')
ylabel('Lick Rate (Hz)')
xlim([-600 1000]); ylim([0 10])

p(1,2).select()
hold on;xlim([-600 1000]);
shadedErrorBar(tt2, mean(fail_TC_long,1), std(fail_TC_long,[],1)./sqrt(size(fail_TC_long,1)), 'b');
shadedErrorBar(tt2, mean(omitR_TC_long,1), std(omitR_TC_long,[],1)./sqrt(size(omitR_TC_long,1)), 'r');
ylabel('dF/F')
p(2,2).select()
hold on
bar(tlick2(1:end-1), mean(faillong_lick_freq_all), 'b'); errorbar(tlick2(1:end-1), mean(faillong_lick_freq_all), std(faillong_lick_freq_all,[],1)/sqrt(size(faillong_lick_freq_all,1)), '.b');
bar(tlick2(1:end-1), mean(omitRlong_lick_freq_all), 'r'); errorbar(tlick2(1:end-1), mean(omitRlong_lick_freq_all), std(omitRlong_lick_freq_all,[],1)/sqrt(size(omitRlong_lick_freq_all,1)), '.r');
xlabel('Time (ms)')
ylabel('Lick Rate (Hz)')
xlim([-600 1000]); ylim([0 10])

saveas(fig, [out_base 'Summary_TC_omitR.fig']);
print([out_base 'Summary_TC_omitR.eps'], '-depsc');
print([out_base 'Summary_TC_omitR.pdf'], '-dpdf');


%%% for IR
fig = figure; p = panel(); p.pack(2,1);
p(1,1).select()
hold on;
shadedErrorBar(tt2, mean(success_TC_long,1), std(success_TC_long,[],1)./sqrt(size(success_TC_long,1)), 'k');
shadedErrorBar(tt2, mean(itiR_TC_long,1)+0.008, std(itiR_TC_long,[],1)./sqrt(size(itiR_TC_long,1)), 'g');
ylabel('dF/F')
xlim([-600 1000]);
title(['Responsive cells- n = ' num2str(total_resp)])
p(2,1).select()
hold on
bar(tlick2(1:end-1), mean(successlong_lick_freq_all), 'k'); %errorbar(tlick2(1:end-1), successlong_lick_freq_all/500, sem_successlong_lick/500, '.k');
bar(tlick2(1:end-1), mean(itiRlong_lick_freq_all), 'g'); %errorbar(tlick2(1:end-1), omitRlong_lick_freq/500, sem_omitRlong_lick/500, '.r');
xlabel('Time (ms)')
ylabel('Lick Rate (Hz)')
xlim([-600 1000]); ylim([0 10])


saveas(fig, [out_base 'Summary_TC_itiR.fig']);
print([out_base 'Summary_TC_itiR.eps'], '-depsc');
print([out_base 'Summary_TC_itiR.pdf'], '-dpdf');


fig = figure;
subplot(1,2,1)
hold on;
shadedErrorBar(tt, mean(licks_TC_notlick,1)- 0.009, std(licks_TC_notlick,[],1)./sqrt(size(licks_TC_notlick,1)), 'b');
shadedErrorBar(tt, mean(fail_TC_notlick,1), std(fail_TC_notlick,[],1)./sqrt(size(fail_TC_notlick,1)), 'r');
shadedErrorBar(tt, mean(success_TC_notlick,1)+0.01, std(success_TC_notlick,[],1)./sqrt(size(success_TC_notlick,1)), 'k');
title(['Not lick but release resp cells- n = ' num2str(sum(cell2mat(notlick_cell_tot))) ' out of total resp cells n = ' num2str(total_resp)])
xlabel('Time (ms)')
ylabel('dF/F')
ylim([-0.01 0.1]);
axis square
% saveas(fig, [out_base 'Summary_SuccessFail_notLick.fig']);
% print([out_base 'Summary_SuccessFail_notLick.eps'], '-depsc');
% print([out_base 'Summary_SuccessFail_notLick.pdf'], '-dpdf');

subplot(1,2,2)
hold on;
shadedErrorBar(tt, mean(licks_TC_lick,1), std(licks_TC_lick,[],1)./sqrt(size(licks_TC_lick,1)), 'b');
shadedErrorBar(tt, mean(fail_TC_lick,1)+0.004, std(fail_TC_lick,[],1)./sqrt(size(fail_TC_lick,1)), 'r');
shadedErrorBar(tt, mean(success_TC_lick,1)+ 0.011, std(success_TC_lick,[],1)./sqrt(size(success_TC_lick,1)), 'k');
title(['Lick resp cells- n = ' num2str(total_lick_cell) ' out of total resp cells n = ' num2str(total_resp)])
xlabel('Time (ms)')
ylabel('dF/F')
ylim([-0.01 0.1]);
axis square
saveas(fig, [out_base 'Summary_SuccessFail_LickNoLick.fig']);
print([out_base 'Summary_SuccessFail_LickNoLick.eps'], '-depsc');
print([out_base 'Summary_SuccessFail_LickNoLick.pdf'], '-dpdf');

fig = figure;
hold on
for ic = 1:size(mouseID,2)
    scatter(20*ones(size(success_lick_corr_sig,1),1), success_lick_corr_sig, 20, 'MarkerFaceColor', [0, 0, 0],'MarkerEdgeColor', [0 0 0]);
    scatter(20*ones(size(success_lick_corr_notsig,1),1), success_lick_corr_notsig, 20, 'MarkerEdgeColor', [0.5 0.5 0.5]);
    scatter(40*ones(size(fail_lick_corr_sig,1),1), fail_lick_corr_sig, 20, 'MarkerFaceColor', [0, 0, 0], 'MarkerEdgeColor', [0 0 0]);
    scatter(40*ones(size(fail_lick_corr_notsig,1),1), fail_lick_corr_notsig, 20, 'MarkerEdgeColor', [0.5 0.5 0.5]);
end
col = [0.9 0 0];
errorbarxy(20, mean(success_lick_corr(:,1)),  0, std(success_lick_corr(:,1))/sqrt(size(success_lick_corr,1)),{'o', col, col, col});
hold on
errorbarxy(40, mean(fail_lick_corr(:,1)),  0, std(fail_lick_corr(:,1))/sqrt(size(fail_lick_corr,1)),{'o', col, col, col});
xlim([0 60])
set(gca, 'XTick', [20, 40], 'XTickLabel', {'correct', 'early'});
ylabel('correlation coeffients')
saveas(fig, [out_base '_Summary_Correlation_lick_dFF_scatter.fig']);
print([out_base '_Summary_Correlation_lick_dFF_scatter.eps'], '-depsc');
print([out_base '_Summary_Correlation_lick_dFF_scatter.pdf'], '-dpdf');

fig = figure;
subplot(2,2,1);
hold on
for ic = 1:size(mouseID,2)
    scatter(nanmean(success_lick_corr_all{ic},2), nanmean(fail_lick_corr_all{ic},2), 4, 'MarkerEdgeColor', [0.5 0.5 0.5])
end
x = -0.2:0.1:0.4;
y = x; plot(x,y, '-k');
xlabel('correlation between success and lick rate');
ylabel('correlation between fail and lick rate');
title('mean correlation for each cell');
saveas(fig, [out_base 'Summary_LickRate_dFF_correlation.fig']);
print([out_base 'Summary_LickRate_dFF_correlation.eps'], '-depsc');
print([out_base 'Summary_LickRate_dFF_correlation.pdf'], '-dpdf');

subplot(2,2,2);
hold on
for ic = 1:size(mouseID,2)
    scatter(nanmean(success_lick_corr_all{ic},2), nanmean(omitR_lick_corr_all{ic},2), 4, 'MarkerEdgeColor', [0.5 0.5 0.5])
end


col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0
    1 0.5 0.5
    0 0.8 0.2];

fig = figure;
scatter_plot(mouseID, success_lick_corr_all, fail_lick_corr_all, [0.5 0.5 0.5])

fig = figure;
subplot(1,2,1)
scatter_plot(mouseID, success_resp_RS_2, lick_rate_success, [0.5 0.5 0.5]);
hold on
scatter_plot(mouseID, fail_resp_RS_2, lick_rate_fail, [0.9 0 0]);
xlabel('dF/F'); ylabel('lick rate(HZ)');
title(['Release resp cells - ', num2str(total_resp)]);

subplot(1,2,2)
col_mat = repmat([0.5 0.5 0.5], id,1);
scatter_plot(mouseID, success_resp_RS_2, lick_rate_success, col_mat);
hold on
col_mat = repmat([0.9 0 0], id,1);
scatter_plot(mouseID, fail_resp_RS_2, lick_rate_fail, col_mat);
xlabel('dF/F'); ylabel('lick rate(HZ)');
title('Release resp cells');

saveas(fig, [out_base 'Summary_LickRate_dFF_scatter_trial.fig']);
print([out_base 'Summary_LickRate_dFF_scatter_trial.eps'], '-depsc');
print([out_base 'Summary_LickRate_dFF_scatter_trial.pdf'], '-dpdf');


fig = figure;
s_resp_RS_2 = cell2mat(success_resp_RS_2);
f_resp_RS_2 = cell2mat(fail_resp_RS_2);
s_lickRate = cell2mat(lick_rate_success);
f_lickRate = cell2mat(lick_rate_fail);
col_mat2 = [0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2];

for i = 1:4
    if sum(f_lickRate == 2*i) < 10
        s_resp_RS_2_tot{i} = [];
        f_resp_RS_2_tot{i} = [];
    else
        s_resp_RS_2_tot{i} = s_resp_RS_2(s_lickRate == 2*i);
        f_resp_RS_2_tot{i} = f_resp_RS_2(f_lickRate == 2*i);
    end
end

hold on
scatter_plot({'2','4','6','8'}, s_resp_RS_2_tot, f_resp_RS_2_tot, col_mat2);
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'o', 'MarkerFaceColor', col_mat2(1,:), 'MarkerEdgeColor', col_mat2(1,:));
h(2) = plot(NaN,NaN,'o', 'MarkerFaceColor', col_mat2(2,:), 'MarkerEdgeColor', col_mat2(2,:));
h(3) = plot(NaN,NaN,'o', 'MarkerFaceColor', col_mat2(3,:), 'MarkerEdgeColor', col_mat2(3,:));
h(4) = plot(NaN,NaN,'o', 'MarkerFaceColor', col_mat2(4,:), 'MarkerEdgeColor', col_mat2(4,:));
legend(h, '2Hz','4Hz','6Hz', '8Hz');

x = 0:0.01:0.08;
y = x;
plot(x,y,'k-');

xlabel('Correct dF/F'); ylabel('Early dF/F');
title('Response Amplitude for Correct and Early with Matching Lick Rate');
saveas(fig, [out_base 'Summary_resp_amp_lickRate.fig']);
print([out_base 'Summary_resp_amp_lickRate.eps'], '-depsc');
print([out_base 'Summary_resp_amp_lickRate.pdf'], '-dpdf');

fig = figure;
subplot(1,2,1)
hold on
col_mat = repmat([0.5 0.5 0.5], id,1);
scatter_plot(mouseID, s_resp_RS_2_wLick, s_resp_RS_2_wNoLick, col_mat);
x = -0.05:0.01:0.15; y = x;
plot(x,y,'k-');
xlabel('Correct dF/F with lick'); ylabel('Correct dF/F without lick');

subplot(1,2,2)
hold on
col_mat = repmat([0.9 0 0], id,1);
scatter_plot(mouseID, f_resp_RS_2_wLick, f_resp_RS_2_wNoLick, col_mat);
x = -0.05:0.01:0.15; y = x;
plot(x,y,'k-');
xlabel('Early dF/F with lick'); ylabel('Early dF/F without lick');
supertitle('Response Amplitude for Correct and Early with/without Lick');
saveas(fig, [out_base 'Summary_resp_amp_lickNolick.fig']);
print([out_base 'Summary_resp_amp_lickNolick.eps'], '-depsc');
print([out_base 'Summary_resp_amp_lickNolick.pdf'], '-dpdf');

fig =figure;
hold on
for i = 1:size(mouseID,2)
    scatter(success_resp_mean{i}, early_percent{i}, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
    scatter(fail_resp_mean{i}, early_percent{i},'MarkerFaceColor', [0.9 0 0], 'MarkerEdgeColor', [0.9 0 0]);
end
% x = 0:1:45; y = x;
% plot(x,y,'k-');
xlabel('dF/F'); ylabel('%early trial with licks')
saveas(fig, [out_base 'Summary_resp_amp_earlyLick.fig']);
print([out_base 'Summary_resp_amp_earlyLick.eps'], '-depsc');
print([out_base 'Summary_resp_amp_earlyLick.pdf'], '-dpdf');


fig = figure;
subplot(1,2,1)
scatter_plot(mouseID, success_resp_RS, mean_lickRate_succ, [0.5 0.5 0.5]);
hold on
scatter_plot(mouseID, fail_resp_RS, mean_lickRate_fail, [0.9 0 0]);
xlabel('dF/F'); ylabel('lick rate(HZ)');
title(['Release resp cells - ', num2str(total_resp)]);

subplot(1,2,2)
col_mat = repmat([0.5 0.5 0.5], id,1);
scatter_plot(mouseID, success_resp_RS, mean_lickRate_succ, col_mat);
hold on
col_mat = repmat([0.9 0 0], id,1);
scatter_plot(mouseID, fail_resp_RS, mean_lickRate_fail, col_mat);
xlabel('dF/F'); ylabel('lick rate(HZ)');
title('Release resp cells');

saveas(fig, [out_base 'Summary_LickRate_dFF_scatter_cell.fig']);
print([out_base 'Summary_LickRate_dFF_scatter_cell.eps'], '-depsc');
print([out_base 'Summary_LickRate_dFF_scatter_cell.pdf'], '-dpdf');


%

% plot histogram of lick-triggered avg dFF
fig = figure;
subplot(1,2,1)
histogram(cell2mat(lickb_resp_all));
xlabel('lick bout triggered Peak dF/F'); ylabel('cell count');
title('lick bout');

subplot(1,2,2)
histogram(cell2mat(licks_resp_all));
xlabel('single lick triggered Peak dF/F'); ylabel('cell count');
title('single licks');

saveas(fig, [out_base 'Summary_lick_resp_hist.fig']);
print([out_base 'Summary_lick_resp_hist.eps'], '-depsc');
print([out_base 'Summary_lick_resp_hist.pdf'], '-dpdf');



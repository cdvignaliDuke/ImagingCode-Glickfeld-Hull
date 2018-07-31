%% second part of analysis for licking
% scatter for lick rate and dF/F
clear
file_info
out_base = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake');
mouseID = mouseID(21:end);
date = date(21:end);
success_TC_notlick = []; success_TC_lick = []; 
licks_TC_lick = []; licks_TC_notlick = []; licks_TC = [];
lickb_TC_lick = []; lickb_TC_notlick = []; lickb_TC = [];
lickb_resp_lick = []; licks_resp_lick =[]; lickb_resp_notlick = []; licks_resp_notlick =[];
fail_TC_notlick = []; fail_TC_lick = [];
early_fail_TC_notlick = []; early_fail_TC_lick = [];
late_fail_TC_notlick = []; late_fail_TC_lick = [];
itiR_TC_notlick = []; itiR_TC_lick = []; itiR_TC = [];
licks_TC_notlick_IR = []; licks_TC_lick_IR = []; licks_TC_IR = [];
lickb_TC_notlick_IR = []; lickb_TC_lick_IR = []; lickb_TC_IR = [];
fail_lick_corr_sig = []; fail_lick_corr_notsig = [];fail_lick_corr_p = [];
success_lick_corr_sig = []; success_lick_corr_notsig = []; success_lick_corr_p = [];
success_lick_freq_all = []; successlong_lick_freq_all =[]; success_lick_freq_IR = [];
succ_hold_dur_all = []; fail_hold_dur_all = []; success_lick_all = []; fail_lick_all = [];
early_fail_lick_freq_all = []; late_fail_lick_freq_all =[]; early_fail_TC_RS = []; late_fail_TC_RS = [];
success_TC = []; success_TC_long = [];
success_TC_IR_lick = []; success_TC_IR_notlick = []; success_TC_IR = [];
faillong_lick_freq_all = []; fail_TC_long = [];
omitR_lick_freq_all = []; omitR_TC_long = []; omitRlong_lick_freq_all = []; omitR_TC = [];
omitR_prevR_lick_freq_all = []; omitR_prevNR_lick_freq_all = [];
itiRlong_lick_freq_all = []; itiR_TC_long = [];
omitR_TC_prevR_all = []; omitR_TC_prevNR_all = [];
resp_success_notlick = []; resp_success_lick = []; 
resp_fail_notlick = []; resp_fail_lick = [];
resp_early_fail_notlick = []; resp_early_fail_lick = [];
resp_late_fail_notlick = []; resp_late_fail_lick = [];
resp_licks_notlick = []; resp_licks_lick = []; resp_licks_all = [];
resp_lickb_notlick = []; resp_lickb_lick = []; resp_lickb_all = [];
resp_success_IR = []; resp_IR = [];
for id =1:size(mouseID,2) %1:10%[2,4,6,11,9,10, 12, 13,15, 16] lick rate%[12 14 16]%[11 13 15]%1:size(mouseID,2)
    id
    
    for rID  = 1:2
        dest_sub  = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        %         dest_sub = ['Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure\', date{id}, '_', runID{rID}, '_', mouseID{id}, '\'];
        if exist(dest_sub)
            load([dest_sub '_cell_TCs_LG.mat']);
            load([dest_sub '_cell_resp_LG.mat']);
            %             load([dest_sub '_cell_resp_amp.mat']);
            load([dest_sub '_cell_categories_LG.mat']);
            load([dest_sub '_release_movies.mat']);
            load([dest_sub 'parse_behavior_LG.mat']);
            load([dest_sub '_spike_variance.mat']);
            load([dest_sub '_dff_lick_LG.mat']);
            load([dest_sub '_lick_stats_LG.mat']);
            load([dest_sub '_lick_movies_LG.mat']);
            %             load([dest_sub '_omitR_previousTrial.mat'])
            
            ncells{id} = size(press_resp,2);
            
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
            lickandresp_cells{id} = intersect(RS_cells{id}, lick_cells);
            respnotlick_cells{id} = intersect(RS_cells{id}, notlick_cells);
            lick_resp_cells{id} = lick_cells;
            notlick_cells_all{id} = notlick_cells;
            
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

            success_TC_notlick = [success_TC_notlick; avg_success(respnotlick_cells{id},:)];
            success_TC_lick = [success_TC_lick; avg_success(lickandresp_cells{id}, :)];
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
           
            success_OR_resp_RS_lateWin{id} = mean((success_resp2(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            success_OR_resp_RS_earlyWin{id} = mean((success_resp1(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            
            fail_TC_notlick = [fail_TC_notlick; avg_fail(respnotlick_cells{id},:)];
            fail_TC_lick = [fail_TC_lick; avg_fail(lickandresp_cells{id}, :)];
            early_fail_TC_notlick = [early_fail_TC_notlick; avg_early_fail(respnotlick_cells{id},:)];
            early_fail_TC_lick = [early_fail_TC_lick; avg_early_fail(lickandresp_cells{id}, :)];
            late_fail_TC_notlick = [late_fail_TC_notlick; avg_late_fail(respnotlick_cells{id},:)];
            late_fail_TC_lick = [late_fail_TC_lick; avg_late_fail(lickandresp_cells{id}, :)];
            [fail_lick_freq, ~] = binLicks(fail_lick_long, double(ifi));
            faillong_lick_freq_all = [faillong_lick_freq_all; fail_lick_freq];
            fail_TC_long = [fail_TC_long; avg_fail_long(RS_cells{id},:)];
            
            if ~isempty(avg_OR)
                omitR_TC_RS = [omitR_TC; avg_OR(RS_cells{id},:)];
                
                [omitR_lick_freq, ~] = binLicks(omitR_lick_long, double(ifi));
                omitRlong_lick_freq_all = [omitRlong_lick_freq_all; omitR_lick_freq];
                omitR_TC_long = [omitR_TC_long; avg_OR_long(RS_cells{id},:)];
                omitR_resp_RS_lateWin{id} = mean((OR_resp2(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
                omitR_resp_RS_earlyWin{id} = mean((OR_resp1(:,RS_cells{id})-OR_base(:,RS_cells{id})),1);
            end
            
            if ~isempty(avg_IR)
                itiR_TC_notlick = [itiR_TC_notlick; avg_IR(respnotlick_cells{id},:)];
                itiR_TC_lick = [itiR_TC_lick; avg_IR(lickandresp_cells{id}, :)];
                itiR_TC = [itiR_TC; avg_IR(RS_cells{id}, :)];
                [itiR_lick_freq, ~] = binLicks(itiR_lick_long, double(ifi));
                itiRlong_lick_freq_all = [itiRlong_lick_freq_all; itiR_lick_freq];
                itiR_TC_long = [itiR_TC_long; avg_IR_long(RS_cells{id},:)];
                
                licks_TC_lick_IR = [licks_TC_lick_IR; avg_single_lick(lickandresp_cells{id},:)];
                licks_TC_notlick_IR = [licks_TC_notlick_IR; avg_single_lick(respnotlick_cells{id},:)];
                licks_TC_IR = [licks_TC_IR; avg_single_lick(RS_cells{id},:)];
                lickb_TC_lick_IR = [lickb_TC_lick_IR; avg_lick(lickandresp_cells{id},:)];
                lickb_TC_notlick_IR = [lickb_TC_notlick_IR; avg_lick(respnotlick_cells{id},:)];
                lickb_TC_IR = [lickb_TC_IR; avg_lick(RS_cells{id},:)];
                
                [success_lick_freq, ~] = binLicks(success_lick, double(ifi));
                success_lick_freq_IR = [success_lick_freq_IR; success_lick_freq];
                
                success_TC_IR_lick = [success_TC_IR_lick; avg_success(lickandresp_cells{id},:)];
                success_TC_IR_notlick = [success_TC_IR_notlick; avg_success(respnotlick_cells{id},:)];
                success_TC_IR = [success_TC_IR; avg_success(RS_cells{id},:)];
                
                resp_success_IR = [resp_success_IR mean(success_resp(:,RS_cells{id})-success_base(:,RS_cells{id}),1)];
                resp_IR = [resp_IR mean(IR_resp(:,RS_cells{id})-IR_base(:,RS_cells{id}),1)];
            end
            
            early_fail_TC_RS = [early_fail_TC_RS; avg_early_fail(RS_cells{id},:)];
            late_fail_TC_RS = [late_fail_TC_RS; avg_late_fail(RS_cells{id},:)];
            
            [early_fail_lick_freq, ~] = binLicks(early_fail_lick, double(ifi));
            early_fail_lick_freq_all = [early_fail_lick_freq_all; early_fail_lick_freq];
            
            [late_fail_lick_freq, ~] = binLicks(late_fail_lick, double(ifi));
            late_fail_lick_freq_all = [late_fail_lick_freq_all; late_fail_lick_freq];
            
            
            licks_TC = [licks_TC; avg_single_lick(RS_cells{id},:)];
            lickb_TC = [lickb_TC; avg_lick(RS_cells{id},:)];
            licks_TC_lick = [licks_TC_lick; avg_single_lick(lickandresp_cells{id},:)];
            lickb_TC_lick = [lickb_TC_lick; avg_lick(lickandresp_cells{id},:)];
            lickb_TC_notlick = [lickb_TC_notlick; avg_lick(respnotlick_cells{id},:)];
            licks_TC_notlick = [licks_TC_notlick; avg_single_lick(respnotlick_cells{id},:)];

            licks_movie{id} = single_lick_movie(:,lick_resp_cells{id},:,:);
            licks_movie_notlick{id} = single_lick_movie(:,notlick_cells_all{id},:,:);

            lickb_resp_lick = [lickb_resp_lick mean(lickb_resp(:,lickandresp_cells{id})-lickb_base(:,lickandresp_cells{id}),1)];
            licks_resp_lick = [licks_resp_lick mean(licks_resp(:,lickandresp_cells{id})-licks_base(:,lickandresp_cells{id}),1)];
            lickb_resp_notlick = [lickb_resp_notlick mean(lickb_resp(:,respnotlick_cells{id})-lickb_base(:,respnotlick_cells{id}),1)];
            licks_resp_notlick = [licks_resp_notlick mean(licks_resp(:,respnotlick_cells{id})-licks_base(:,respnotlick_cells{id}),1)];
            
            resp_success_lick = [resp_success_lick mean(success_resp(:,lickandresp_cells{id})-success_base(:,lickandresp_cells{id}),1)];
            resp_success_notlick = [resp_success_notlick mean(success_resp(:,respnotlick_cells{id})-success_base(:,respnotlick_cells{id}),1)];
            
            resp_fail_lick = [resp_fail_lick mean(fail_resp(:,lickandresp_cells{id})-fail_base(:,lickandresp_cells{id}),1)];
            resp_fail_notlick = [resp_fail_notlick mean(fail_resp(:,respnotlick_cells{id})-fail_base(:,respnotlick_cells{id}),1)];
            
            resp_early_fail_lick = [resp_early_fail_lick mean(early_fail_resp(:,lickandresp_cells{id})-early_fail_base(:,lickandresp_cells{id}),1)];
            resp_early_fail_notlick = [resp_early_fail_notlick mean(early_fail_resp(:,respnotlick_cells{id})-early_fail_base(:,respnotlick_cells{id}),1)];
            
            resp_late_fail_lick = [resp_late_fail_lick mean(late_fail_resp(:,lickandresp_cells{id})-late_fail_base(:,lickandresp_cells{id}),1)];
            resp_late_fail_notlick = [resp_late_fail_notlick mean(late_fail_resp(:,respnotlick_cells{id})-late_fail_base(:,respnotlick_cells{id}),1)];
           
        end
    end
end
%% figures
total_resp = sum(cell2mat(tot_resp));

tt = (-pre_release_frames:post_release_frames)*double(ifi);
tt_lick_ind = (1:3:45); 
tt_lick_bin = tt(tt_lick_ind);

%scatter of correct vs early 
figure;
subplot(3,2,1)
scatter(resp_success_notlick, resp_fail_notlick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_notlick,2), mean(resp_fail_notlick,2), std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2)), std(resp_fail_notlick,[],2)./sqrt(size(resp_fail_notlick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Early dF/F')
axis square
[h p] = ttest(resp_success_notlick,resp_fail_notlick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_notlick,2)) '   '   num2str(mean(resp_fail_notlick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2))) '   '   num2str(std(resp_fail_notlick,[],2)./sqrt(size(resp_fail_notlick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Release not lick resp cells-', ['n = ' num2str(size(resp_success_notlick,2)) ' out of total resp cells n = ' num2str(total_resp)]})

subplot(3,2,2)
scatter(resp_success_lick, resp_fail_lick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_lick,2), mean(resp_fail_lick,2), std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2)), std(resp_fail_lick,[],2)./sqrt(size(resp_fail_lick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Early dF/F')
axis square
[h p] = ttest(resp_success_lick,resp_fail_lick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_lick,2)) '   '   num2str(mean(resp_fail_lick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2))) '   '   num2str(std(resp_fail_lick,[],2)./sqrt(size(resp_fail_lick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Lick and release resp cells-', ['n = ' num2str(size(resp_success_lick,2)) ' out of total resp cells n = ' num2str(total_resp)]})

subplot(3,2,3)
scatter(resp_success_notlick, resp_fail_notlick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_notlick,2), mean(resp_early_fail_notlick,2), std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2)), std(resp_early_fail_notlick,[],2)./sqrt(size(resp_early_fail_notlick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Early Early dF/F')
axis square
[h p] = ttest(resp_success_notlick,resp_early_fail_notlick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_notlick,2)) '   '   num2str(mean(resp_early_fail_notlick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2))) '   '   num2str(std(resp_early_fail_notlick,[],2)./sqrt(size(resp_early_fail_notlick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Release not lick resp cells-', ['n = ' num2str(size(resp_success_notlick,2)) ' out of total resp cells n = ' num2str(total_resp)]})

subplot(3,2,4)
scatter(resp_success_lick, resp_early_fail_lick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_lick,2), mean(resp_early_fail_lick,2), std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2)), std(resp_early_fail_lick,[],2)./sqrt(size(resp_early_fail_lick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Early Early dF/F')
axis square
[h p] = ttest(resp_success_lick,resp_early_fail_lick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_lick,2)) '   '   num2str(mean(resp_early_fail_lick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2))) '   '   num2str(std(resp_early_fail_lick,[],2)./sqrt(size(resp_early_fail_lick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Lick and release resp cells-', ['n = ' num2str(size(resp_success_lick,2)) ' out of total resp cells n = ' num2str(total_resp)]})

subplot(3,2,5)
scatter(resp_success_notlick, resp_late_fail_notlick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_notlick,2), mean(resp_late_fail_notlick,2), std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2)), std(resp_late_fail_notlick,[],2)./sqrt(size(resp_late_fail_notlick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Late Early dF/F')
axis square
[h p] = ttest(resp_success_notlick,resp_late_fail_notlick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_notlick,2)) '   '   num2str(mean(resp_late_fail_notlick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_notlick,[],2)./sqrt(size(resp_success_notlick,2))) '   '   num2str(std(resp_late_fail_notlick,[],2)./sqrt(size(resp_late_fail_notlick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Release not lick resp cells-', ['n = ' num2str(size(resp_success_notlick,2)) ' out of total resp cells n = ' num2str(total_resp)]})

subplot(3,2,6)
scatter(resp_success_lick, resp_late_fail_lick,10,[0.5 0.5 0.5])
hold on
errorbarxy(mean(resp_success_lick,2), mean(resp_late_fail_lick,2), std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2)), std(resp_late_fail_lick,[],2)./sqrt(size(resp_late_fail_lick,2)),{'r', 'r', 'r'})
xlim([-0.2 0.4])
ylim([-0.2 0.4])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Late Early dF/F')
axis square
[h p] = ttest(resp_success_lick,resp_late_fail_lick);
text(-0.15,0.35, ['mean = ' num2str(mean(resp_success_lick,2)) '   '   num2str(mean(resp_late_fail_lick,2))])
text(-0.15,0.30, ['sem = ' num2str(std(resp_success_lick,[],2)./sqrt(size(resp_success_lick,2))) '   '   num2str(std(resp_late_fail_lick,[],2)./sqrt(size(resp_late_fail_lick,2)))])
text(-0.15,0.25, ['p = ' num2str(chop(p,3))])
title({'Lick and release resp cells-', ['n = ' num2str(size(resp_success_lick,2)) ' out of total resp cells n = ' num2str(total_resp)]})
print([out_base '\Summary_SuccessFail_LickNoLick_scatter_LG.pdf'], '-dpdf','-fillpage');

%summary of ITI reward
figure;
subplot(2,3,1)
shadedErrorBar(tt, mean(success_TC_IR_notlick-mean(success_TC_IR_notlick(:,1:5),2),1),std(success_TC_IR_notlick,[],1)./sqrt(size(success_TC_IR_notlick,1)),'k');
hold on
shadedErrorBar(tt, mean(itiR_TC_notlick-mean(itiR_TC_notlick(:,1:5),2),1),std(itiR_TC_notlick,[],1)./sqrt(size(itiR_TC_notlick,1)),'g');
xlim([-500 500])
ylim([-0.01 0.1])
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title({'Not lick resp-', ['n = ' num2str(size(itiR_TC_notlick,1))]})

subplot(2,3,2)
shadedErrorBar(tt, mean(success_TC_IR_lick-mean(success_TC_IR_lick(:,1:5),2),1),std(success_TC_IR_lick,[],1)./sqrt(size(success_TC_IR_lick,1)),'k');
hold on
shadedErrorBar(tt, mean(itiR_TC_lick-mean(itiR_TC_lick(:,1:5),2),1),std(itiR_TC_lick,[],1)./sqrt(size(itiR_TC_lick,1)),'g');
xlim([-500 500])
ylim([-0.01 0.1])
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title({'Lick resp-', ['n = ' num2str(size(itiR_TC_lick,1))]})

subplot(2,3,3)
shadedErrorBar(tt, mean(success_TC_IR-mean(success_TC_IR(:,1:5),2),1),std(success_TC_IR,[],1)./sqrt(size(success_TC_IR,1)),'k');
hold on
shadedErrorBar(tt, mean(itiR_TC-mean(itiR_TC(:,1:5),2),1),std(itiR_TC,[],1)./sqrt(size(itiR_TC,1)),'g');
xlim([-500 500])
ylim([-0.01 0.1])
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title({'All resp-', ['n = ' num2str(size(itiR_TC,1))]})

subplot(2,3,4)
scatter(resp_success_IR, resp_IR,10,[0.5 0.5 0.5]);
hold on
errorbarxy(mean(resp_success_IR,2), mean(resp_IR,2), std(resp_success_IR,[],2)./sqrt(size(resp_success_IR,2)), std(resp_IR,[],2)./sqrt(size(resp_IR,2)),{'r', 'r', 'r'});
xlim([-0.05 0.3])
ylim([-0.05 0.3])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Correct dF/F')
ylabel('Unexpected Reward dF/F')
axis square
[h p] = ttest(resp_success_IR, resp_IR);
title(['p = ' num2str(chop(p,3))])

subplot(2,3,5)
bar(tt_lick_bin,mean(success_lick_freq_IR,1),'k')
hold on
errorbar(tt_lick,mean(success_lick_freq_IR,1), std(success_lick_freq_IR,[],1)./sqrt(size(success_lick_freq_IR,1)),'-k')
xlim([-500 500])
ylim([-2 15])
axis square
xlabel('Time (ms)')
ylabel('Lick rate (Hz)')

subplot(2,3,6)
bar(tt_lick_bin(1:11),mean(itiRlong_lick_freq_all(:,1:11),1),'g')
hold on
errorbar(tt_lick(1:11),mean(itiRlong_lick_freq_all(:,1:11),1), std(itiRlong_lick_freq_all(:,1:11),[],1)./sqrt(size(itiRlong_lick_freq_all,1)),'-g')
xlim([-500 500])
ylim([-2 15])
axis square
xlabel('Time (ms)')
ylabel('Lick rate (Hz)')

print([out_base '\Summary_itiRvsSpontLick_LickNoLick_TC_LG.pdf'], '-dpdf','-fillpage');

%single lick vs lick bout
ttl = (-round(2500/double(ifi)):round(2500/double(ifi)))*double(ifi);
[tts ttl_ind] = intersect(ttl,tt);
figure;
subplot(2,2,1)
shadedErrorBar(tts, mean(licks_TC_lick(:,ttl_ind)-mean(licks_TC_lick(:,ttl_ind(1:5)),2),1),std(licks_TC_lick(:,ttl_ind),[],1)./sqrt(size(licks_TC_lick,1)),'k');
hold on
shadedErrorBar(tts, mean(lickb_TC_lick(:,ttl_ind)-mean(lickb_TC_lick(:,ttl_ind(1:5)),2),1),std(lickb_TC_lick(:,ttl_ind),[],1)./sqrt(size(lickb_TC_lick,1)),'b');
xlim([-500 500])
ylim([-0.03 0.1])
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title({'Lick resp-', ['n = ' num2str(size(lickb_TC_lick,1))]})

subplot(2,2,2)
shadedErrorBar(tts, mean(licks_TC_notlick(:,ttl_ind)-mean(licks_TC_notlick(:,ttl_ind(1:5)),2),1),std(licks_TC_notlick(:,ttl_ind),[],1)./sqrt(size(licks_TC_notlick,1)),'k');
hold on
shadedErrorBar(tts, mean(lickb_TC_notlick(:,ttl_ind)-mean(lickb_TC_notlick(:,ttl_ind(1:5)),2),1),std(lickb_TC_notlick(:,ttl_ind),[],1)./sqrt(size(lickb_TC_notlick,1)),'b');
xlim([-500 500])
ylim([-0.01 0.1])
xlabel('Time (ms)')
ylabel('dF/F')
axis square
title({'Not Lick resp-', ['n = ' num2str(size(lickb_TC_notlick,1))]})

subplot(2,2,3)
scatter(licks_resp_lick, lickb_resp_lick,10,[0.5 0.5 0.5]);
hold on
errorbarxy(mean(licks_resp_lick,2), mean(lickb_resp_lick,2), std(licks_resp_lick,[],2)./sqrt(size(licks_resp_lick,2)), std(lickb_resp_lick,[],2)./sqrt(size(lickb_resp_lick,2)),{'r', 'r', 'r'});
xlim([-0.1 0.3])
ylim([-0.1 0.3])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Single lick dF/F')
ylabel('Lick bout dF/F')
axis square
[h p] = ttest(licks_resp_lick, lickb_resp_lick);
title(['Lick resp - p = ' num2str(chop(p,3))])

subplot(2,2,4)
scatter(licks_resp_notlick, lickb_resp_notlick,10,[0.5 0.5 0.5]);
hold on
errorbarxy(mean(licks_resp_notlick,2), mean(lickb_resp_notlick,2), std(licks_resp_notlick,[],2)./sqrt(size(licks_resp_notlick,2)), std(lickb_resp_notlick,[],2)./sqrt(size(lickb_resp_notlick,2)),{'r', 'r', 'r'});
xlim([-0.1 0.3])
ylim([-0.1 0.3])
hline(0,'--k')
vline(0,'--k')
refline(1,0)
xlabel('Single lick dF/F')
ylabel('Lick bout dF/F')
axis square
[h p] = ttest(licks_resp_notlick, lickb_resp_notlick);
title(['Not lick resp - p = ' num2str(chop(p,3))])

print([out_base '\Summary_SingleVsLickBout_LG.pdf'], '-dpdf','-fillpage');



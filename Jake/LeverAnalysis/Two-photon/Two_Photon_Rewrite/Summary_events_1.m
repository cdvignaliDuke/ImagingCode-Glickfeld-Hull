%compares event magnitudes, event rates(PSTHs), event synchrony across conditions

clear
file_info
% out_base = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis\');
%out_base = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis\');
out_base = fullfile('Z:\Analysis\2P Analysis\Lever\');
count_base = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis\');
mouseID = mouseID(1:30);
success_nevents_tot = 0; fail_nevents_tot = 0; spont_nevents_tot = 0;
fail_ntrials = 0; success_ntrials = 0;
for id = 1:30%size(mouseID,2)
    id
    for rID  = 1:2
%         dest_sub  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
%         dest_sub = ['Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure\', date{id}, '_', runID{rID}, '_', mouseID{id}, '\'];
        dest_sub  = fullfile('\\crash.dhe.duke.edu\data\home\ziye\2P_Analysis\2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        
        if exist(dest_sub)
            load([dest_sub 'parse_behavior']);
            load([dest_sub '_cell_categories2.mat']);           
            load([dest_sub '_event_peaks.mat']);
            load([dest_sub '_spont_events.mat']);
            load([dest_sub '_event_hist.mat']);
            load([dest_sub '_evoked_events.mat']);
            load([dest_sub '_norm2spont.mat']);
            load([count_base, 'cell_count.mat']);
            load([dest_sub '_event_summary.mat']);
            
            ncells(id) = size(press,2);
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
            
            TC_length(id) = size(all_success_event,2);
            TC_ifi(id) = ifi;
            
            press_event_TC{id} = all_press_event;
            success_event_TC{id} = all_success_event;
            fail_event_TC{id} = all_fail_event;
%             spont_event_TC{id} = all_spont_event;
            press_event_peak{id} = all_press_peak;
            release_event_peak{id} = all_release_peak;
            success_event_peak{id} = all_success_peak;
            fail_event_peak{id} = all_fail_peak;
            
            press_event_TC_RS{id} = all_press_event(RS_cells{id},:);
            success_event_TC_RS{id} = all_success_event(RS_cells{id},:);
            fail_event_TC_RS{id} = all_fail_event(RS_cells{id},:);
%             spont_event_TC_RS{id} = all_spont_event(RS_cells{id},:);
            press_event_peak_RS{id} = all_press_peak(:,RS_cells{id});
            release_event_peak_RS{id} = all_release_peak(:,RS_cells{id});
            success_event_peak_RS{id} = all_success_peak(:,RS_cells{id});
            fail_event_peak_RS{id} = all_fail_peak(:,RS_cells{id});
            
            success_hist{id} = all_success_hist.*(1000/double(ifi)); %convert to Hz
            fail_hist{id} = all_fail_hist.*(1000/double(ifi));
            success_hist_RS{id} = resp_success_hist.*(1000/double(ifi));
            fail_hist_RS{id} = resp_fail_hist.*(1000/double(ifi));
            
%             psuccess_hist{id} = all_psuccess_hist.*(1000/double(ifi));
%             pfail_hist{id} = all_pfail_hist.*(1000/double(ifi));
            psuccess_EndF = trial_outcome.success_ptimeEndF;
            pfail_EndF = trial_outcome.early_ptimeEndF;
            stop_frame = max([psuccess_EndF pfail_EndF]);
%             for ip = 1:length(psuccess_EndF)
%                 resp_psuccess_hist(psuccess_EndF(ip)+pre_buffer+2:end,ip) = nan;
%             end
%             for ip = 1:length(pfail_EndF)
%                 resp_pfail_hist(pfail_EndF(ip)+pre_buffer+2:end,ip) = nan;
%             end
%             resp_psuccess_hist(stop_frame+1 :end,:) = [];
%             resp_pfail_hist(stop_frame+1 :end,:) = [];
%             psuccess_hist_RS{id} = resp_psuccess_hist.*(1000/double(ifi));
%             pfail_hist_RS{id} = resp_pfail_hist.*(1000/double(ifi));
            
            success_rate{id} = success_rate_all.*(1000/double(ifi));
            fail_rate{id} = fail_rate_all.*(1000/double(ifi));
            success_RS_rate{id} = success_rate_RS.*(1000/double(ifi));
            fail_RS_rate{id} = fail_rate_RS.*(1000/double(ifi));
            
            success_ind{id} = (success_ind_all-1).*double(ifi);
            fail_ind{id} = (fail_ind_all-1).*double(ifi);
            success_RS_ind{id} = (success_ind_RS-1).*double(ifi);
            fail_RS_ind{id} = (fail_ind_RS-1).*double(ifi);
            
            release_norm_all{id} = release_norm_avg;
            press_norm_all{id} = press_norm_avg;
            spont_norm_all{id} = spont_norm_avg;
            success_norm_all{id} = success_norm_avg;
            fail_norm_all{id} = fail_norm_avg;
%             lickb_norm_all{id} = lickb_norm_avg;
%             licks_norm_all{id} = licks_norm_avg;
            
            release_norm_RS{id} = release_norm_avg(:,RS_cells{id});
            press_norm_RS{id} = press_norm_avg(:,RS_cells{id});
            spont_norm_RS{id} = spont_norm_avg(:,RS_cells{id});
            success_norm_RS{id} = success_norm_avg(:,RS_cells{id});
            fail_norm_RS{id} = fail_norm_avg(:,RS_cells{id});
%             lickb_norm_RS{id} = lickb_norm_avg(:,RS_cells{id});
%             licks_norm_RS{id} = licks_norm_avg(:,RS_cells{id});
            
            rate{id} = events_rate;
            
            success_syn_all_RS{id} = success_syn_RS;
            fail_syn_all_RS{id}  = fail_syn_RS;
            success_syn_c_RS{id} = success_syn_cell_RS;
            fail_syn_c_RS{id}  = fail_syn_cell_RS;
            cue_syn_c_RS{id} = cue_syn_cell_RS;
            
            %get all cells 
            success_syn_c{id} = success_syn_cell;
            fail_syn_c{id} = fail_syn_cell;
%             cue_syn_cell

            %under construction ==============================================================================
%             success_syn_all_spike_ind_all{id} = success_syn_all_spike_ind.spike_latency;
%             success_syn_cell_spike_ind_all{id} = success_syn_cell_spike_ind.spike_latency;
%             cue_syn_cell_spike_ind_all{id} = cue_syn_cell_spike_ind.spike_latency; 

            %RT_offset_ind_frames_all{id} = RT_offset_ind_frames;
            %success_syn_RS_spike_ind_all{id} = success_syn_RS_spike_ind; 
            success_syn_cell_RS_spike_ind_all{id} = success_syn_cell_RS_spike_ind; 
            fail_syn_cell_RS_spike_ind_all{id} = fail_syn_cell_RS_spike_ind;
            cue_syn_cell_RS_spike_ind_all{id} = cue_syn_cell_RS_spike_ind;
            
            %===================================================================================================
%             fail_syn_all_spike_ind
%             fail_syn_RS_spike_ind
%             fail_syn_cell_spike_ind
%             fail_syn_cell_RS_spike_ind
            
            ts{id} = [-data_start:data_end].*double(ifi);
            data_start_frame{id} = data_start;
            data_end_frame{id} = data_end;
            sz1 = size(all_success_hist,1);
            th{id} = [1-(sz1/2):(sz1/2)].*double(ifi);
            
            success_nevents_tot = success_nevents_tot + sum(success_nevents);
            fail_nevents_tot = fail_nevents_tot + sum(fail_nevents);
            spont_nevents_tot = spont_nevents_tot + sum(spont_nevents);
            
            success_ntrials = size(success_syn_all,2) + success_ntrials;
            fail_ntrials = size(fail_syn_all,2) + fail_ntrials;
            
            acor = [];
            
            acor_all{id} = iti_cor_mean';
        end
    end
end


num_trials_succ = 0;
num_cells_fail = 0; 
num_cells_succ = 0;
num_trials_fail = 0;
num_cells_succ_all = 0;
num_cells_fail_all = 0;

length(find([fail_syn_c{:}] == 0));
num_trials_succ = sum(~isnan([success_syn_all_RS{:}]));
num_trials_fail = sum( ~isnan([fail_syn_all_RS{:}]));
num_cells_succ =  sum( ~isnan([success_syn_c_RS{:}]));
num_cells_fail = sum( ~isnan([fail_syn_c_RS{:}]));
num_cells_succ_all =  sum( ~isnan([success_syn_c{:}]));
num_cells_fail_all =  sum( ~isnan([fail_syn_c{:}]));

num_cells_succ =  length([success_syn_c_RS{:}]);
num_cells_fail = length([fail_syn_c_RS{:}]);
num_cells_succ_all =  length([success_syn_c{:}]);
num_cells_fail_all =  length([fail_syn_c{:}]);


%%plotting
% ['peach' 'army green' 'yellow' 'purple' 'black']
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
    0.5 0.4 0; 0.5 0.5 0.5; 0.3 0.5 1; 0.1 0.5 0.7; 0 0.6 0.2;0.8 0.8 0.4;0.1 0.1 0.1;0.3 0.7 0; 0 0 0; 0.1 0.5 0.5];

% summary of autocorrelation
stframe_size = cell2mat(cellfun(@size, acor_all, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

acor_all_all = interp_frame(acor_all, max_frame);
a = acor_all_all';
tta = -2:2/(size(a,1)-1):2;
a = [fliplr(mean(a,2)')'; mean(a,2)];
a(size(a,1)/2) = [];
fig = figure;bar(tta, a);hold on;errorbar(tta, a, std(a,[],2)/sqrt(size(a,2)),'.')
xlim([-2 2])
xlabel('Lag (s)')
ylabel('Auto correlation');
title('auto correlation for spont events for all cells');
saveas(fig, [out_base 'Summary_autocorrelation.fig']);
print([out_base 'Summary_autocorrelation.eps'], '-depsc');
print([out_base 'Summary_autocorrelation.pdf'], '-dpdf');

b = acor_all{end}';
ttb = -2:2/(size(b,1)-1):2;
b1 = b(:,55); b2 = b(:,70); b3 = b(:,30); b4 = b(:,20);
b1 = [fliplr(b1')'; b1]; b1(size(b,1))= [];
b2 = [fliplr(b2')'; b2]; b2(size(b,1))= [];
b3 = [fliplr(b3')'; b3]; b3(size(b,1))= [];
b4 = [fliplr(b4')'; b4]; b4(size(b,1))= [];
fig = figure;subplot(4,1,1); bar(ttb, b1);xlim([-2 2])
subplot(4,1,2); bar(ttb, b2); xlim([-2 2])
subplot(4,1,3); bar(ttb, b3);xlim([-2 2])
subplot(4,1,4); bar(ttb, b4);xlim([-2 2])
xlabel('Lag (s)')
ylabel('Auto correlation');
supertitle('example auto correlation for cells from img 59');
saveas(fig, [out_base 'img59_autocorrelation.fig']);
print([out_base 'imd59_autocorrelation.eps'], '-depsc');
print([out_base 'img59_autocorrelation.pdf'], '-dpdf');


% hold time and event rate
sResp = cell2mat(mean_success_rate')';
sResp(isnan(sResp)) = [];
sHold = cell2mat(succ_hold_dur);
sHold(isnan(sHold)) = [];
[LinearCoeffS1, fitS1] = polyfit(sHold, sResp, 1);
CorrSfit1 = polyval(LinearCoeffS1, sHold);
[BS,idx] = histc(sHold,0:0.25:round(max(sHold)));
sHold_bins = accumarray(idx(:),sHold,[],@mean);
sResp_std = accumarray(idx(:),sResp,[],@sem);
sResp_bins = accumarray(idx(:),sResp,[],@mean);
[LinearCoeffS2, fitS2] = polyfit(sHold_bins(BS>=10), sResp_bins(BS>=10), 1);
CorrSfit2 = polyval(LinearCoeffS2, sHold_bins(BS>=10));

fResp = cell2mat(mean_fail_rate')';
fResp(isnan(fResp)) = [];
fHold = cell2mat(fail_hold_dur);
fHold(isnan(fHold)) = [];
[LinearCoeffF1, fitF1] = polyfit(fHold, fResp, 1);
CorrFfit1 = polyval(LinearCoeffF1, fHold);
[BF,idx] = histc(fHold,0:0.25:round(max(fHold)));
fHold_bins = accumarray(idx(:),fHold,[],@mean);
fResp_std = accumarray(idx(:),fHold,[],@sem);
fResp_bins = accumarray(idx(:),fResp,[],@mean);
[LinearCoeffF2, fitF2] = polyfit(fHold_bins(BF>=10), fResp_bins(BF>=10), 1);
CorrFfit2 = polyval(LinearCoeffF2, fHold_bins(BF>=10));


% scatter_plot(mouseID, succ_hold_dur, mean_success_rate, col_mat_s);
% scatter(sHold_bins(BS>=10), sResp_bins(BS>=10), 4, 'MarkerEdgeColor', [0.5 0.5 0.5] );
fig = figure;
hold on
errorbar(sHold_bins(BS>=10), sResp_bins(BS>=10), sResp_std(BS>=10), '.', 'MarkerSize', 10, 'color', [0.5,0.5,0.5])
plot(sHold_bins(BS>=10), CorrSfit2, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);


% scatter_plot(mouseID, fail_hold_dur, mean_fail_rate, col_mat_s);
% scatter(fHold_bins(BF>=10), fResp_bins(BF>=10), 4, 'MarkerEdgeColor', [0.9 0 0] );
errorbar(fHold_bins(BF>=10), fResp_bins(BF>=10), fResp_std(BF>=10), '.', 'MarkerSize', 10, 'color', [0.9,0,0])
plot(fHold_bins(BF>=10), CorrFfit2, 'color', [0.9, 0,0], 'Linewidth', 1.5);

%plot(x,y,'-k')
title('event rate and hold time for responsive cells, each point is one trial, black-correct, red-early');
xlabel('hold time(s)');
ylabel('event rate(Hz)');
saveas(fig, [out_base 'Summary_eventRate_holdtime_trial_avg.fig']);
print([out_base 'Summary_eventRate_holdtime_trial_avg.eps'], '-depsc');
print([out_base 'Summary_eventRate_holdtime_trial_avg.pdf'], '-dpdf');

fig = figure;
col_mat_s = [0.5, 0.5, 0.5];
scatter_plot(mouseID, succ_hold_dur, mean_success_rate, col_mat_s);
% scatter(sHold_bins(3:end-1), sResp_bins(3:end-1), 4, 'MarkerEdgeColor', [0.5 0.5 0.5] );
hold on
plot(sHold, CorrSfit1, 'color', [0.5, 0.5,0.5], 'Linewidth', 1.5);

col_mat_s = [0.9 0 0];
scatter_plot(mouseID, fail_hold_dur, mean_fail_rate, col_mat_s);
plot(fHold, CorrFfit1, 'color', [0.9, 0,0], 'Linewidth', 1.5);

%plot(x,y,'-k')
title('event rate and hold time for responsive cells, each point is one trial, black-correct, red-early');
xlabel('hold time(s)');
ylabel('event rate(Hz)');
saveas(fig, [out_base 'Summary_eventRate_holdtime_trial_scatter.fig']);
print([out_base 'Summary_eventRate_holdtime_trial_scatter.eps'], '-depsc');
print([out_base 'Summary_eventRate_holdtime_trial_scatter.pdf'], '-dpdf');
% col_mat_s = repmat([0.5,0.5,0.5],11,1);

%summary of event rate
fig=figure;
rate_all = [];
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
        plot(id*ones(size(rate{id})), rate{id}, 'o', 'color', col_mat(cid,:));
    
    rate_all = [rate_all rate{id}];
    hold on
end
ylim([0 2])
xlim([0.5 size(mouseID,2)+.5])
ylabel('Spike rate (Hz)')
xlabel('session')
title(['Avg rate: ' num2str(chop(mean(rate_all,2),3)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
saveas(fig, [out_base 'Summary_event_rate.fig']);
print([out_base 'Summary_event_rate.eps'], '-depsc');
print([out_base 'Summary_event_rate.pdf'], '-dpdf');

fig = figure;
histogram(rate_all)
xlabel('Firing rate (HZ)')
ylabel('#Cell')

title(['Avg rate: ' num2str(chop(mean(rate_all,2),3)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
saveas(fig, [out_base 'Summary_event_rate_hist.fig']);
print([out_base 'Summary_event_rate_hist.eps'], '-depsc');
print([out_base 'Summary_event_rate_hist.pdf'], '-dpdf');


%summary of average events
fig=figure;
for id = 1:size(mouseID,2)
    press_events = press_event_TC{id};
    success_events = success_event_TC{id};
    fail_events = fail_event_TC{id};
    cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
    subplot(4,5,id)
    errorbar(ts{id}, nanmean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
    hold on
    errorbar(ts{id}, nanmean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
    hold on
    errorbar(ts{id}, nanmean(fail_events(cell_use,:),1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
    title([date{id} ' ' mouseID{id} '- ' num2str(length(cell_use)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
supertitle(['Good events- All cells- Black: success; Red: early; Cyan: press'])
saveas(fig, [out_base 'Summary_allcells_event_TC.fig']);
print([out_base 'Summary_allcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_allcells_event_TC.pdf'], '-dpdf');
fig=figure;
for id = 1:size(mouseID,2)
    press_events = press_event_TC_RS{id};
    success_events = success_event_TC_RS{id};
    fail_events = fail_event_TC_RS{id};
    cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
    subplot(4,5,id)
    errorbar(ts{id}, nanmean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
    hold on
    errorbar(ts{id}, nanmean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
    hold on
    errorbar(ts{id}, nanmean(fail_events,1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
    title([date{id} ' ' mouseID{id} '- ' num2str(length(cell_use)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
supertitle(['Good events- Resp cells- Black: success; Red: early; Cyan: press'])
saveas(fig, [out_base 'Summary_respcells_event_TC.fig']);
print([out_base 'Summary_respcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_respcells_event_TC.pdf'], '-dpdf');


%% commented to fix acquisition rates
%average across days
fig=figure;
ra = [];
pa = [];
rr = [];
pr = [];
sa = [];
fa = [];
sr = [];
fr = [];
cid = 1;
for id = 1:size(mouseID,2)
    ra = [ra release_event_peak{id}];
    pa = [pa press_event_peak{id}];
    rr = [rr release_event_peak_RS{id}];
    pr = [pr press_event_peak_RS{id}];
    sa = [sa success_event_peak{id}];
    fa = [fa fail_event_peak{id}];
    sr = [sr success_event_peak_RS{id}];
    fr = [fr fail_event_peak_RS{id}];
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
    subplot(2,2,1)
    scatter(release_event_peak{id}, press_event_peak{id}, 'MarkerEdgeColor', col_mat(cid,:));
    hold on
    subplot(2,2,2)
    scatter(success_event_peak{id}, fail_event_peak{id}, 'MarkerEdgeColor', col_mat(cid,:));
    hold on
    subplot(2,2,3)
    scatter(release_event_peak_RS{id}, press_event_peak_RS{id}, 'MarkerEdgeColor', col_mat(cid,:));
    hold on
    subplot(2,2,4)
    scatter(success_event_peak_RS{id}, fail_event_peak_RS{id}, 'MarkerEdgeColor', col_mat(cid,:));
    hold on
    
end
x = 0:.1:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release event peak dF/F')
ylabel('Press event peak dF/F')
[h_rpa p_rpa]= ttest(ra,pa);
title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_rpa,2))])
subplot(2,2,2)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
[h_sfa p_sfa]= ttest(sa,fa);
title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_sfa,2))])
subplot(2,2,3)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release event peak dF/F')
ylabel('Press event peak dF/F')
[h_rpr p_rpr]= ttest(rr,pr);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_rpr,2))])
subplot(2,2,4)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
[h_sfr p_sfr]= ttest(sr,fr);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sfr,2))])
suptitle('Peak event amplitude')
saveas(fig, [out_base 'Summary_event_amp_scatter.fig']);
print([out_base 'Summary_event_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_event_amp_scatter.pdf'], '-dpdf');

%average event amplitude by mouse
fig=figure;
subplot(2,2,1)
scatter_plot(mouseID, release_event_peak, press_event_peak, col_mat);

subplot(2,2,2)
scatter_plot(mouseID, success_event_peak, fail_event_peak, col_mat);

subplot(2,2,3)
scatter_plot(mouseID, release_event_peak_RS, press_event_peak_RS, col_mat);

subplot(2,2,4)
scatter_plot(mouseID, success_event_peak_RS, fail_event_peak_RS, col_mat);

% for id = 1:size(mouseID,2)
%     subplot(2,2,1)
%     rp_peak = [release_event_peak{id}; press_event_peak{id}];
%     rp_use = find(~isnan(mean(rp_peak,1)));
%     errorbarxy(mean(rp_peak(1,rp_use),2), nanmean(rp_peak(2,rp_use),2), std(rp_peak(1,rp_use),[],2)./sqrt(length(rp_use)), std(rp_peak(2,rp_use),[],2)./sqrt(length(rp_use)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(rp_peak(1,rp_use),2), nanmean(rp_peak(2,rp_use),2),col_mat(id,:));
%     subplot(2,2,2)
%     sf_peak = [success_event_peak{id}; fail_event_peak{id}];
%     sf_use = find(~isnan(mean(sf_peak,1)));
%     errorbarxy(mean(sf_peak(1,sf_use),2), nanmean(sf_peak(2,sf_use),2),std(sf_peak(1,sf_use),[],2)./sqrt(length(sf_use)), std(sf_peak(2,sf_use),[],2)./sqrt(length(sf_use)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(sf_peak(1,sf_use),2), nanmean(sf_peak(2,sf_use),2),col_mat(id,:));
%     subplot(2,2,3)
%     rp_peak_RS = [release_event_peak_RS{id}; press_event_peak_RS{id}];
%     rp_use_RS = find(~isnan(mean(rp_peak_RS,1)));
%     errorbarxy(mean(rp_peak_RS(1,rp_use_RS),2), nanmean(rp_peak_RS(2,rp_use_RS),2), std(rp_peak_RS(1,rp_use_RS),[],2)./sqrt(length(rp_use_RS)), std(rp_peak_RS(2,rp_use_RS),[],2)./sqrt(length(rp_use_RS)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(rp_peak_RS(1,rp_use_RS),2), nanmean(rp_peak_RS(2,rp_use_RS),2), col_mat(id,:));
%     subplot(2,2,4)
%     sf_peak_RS = [success_event_peak_RS{id}; fail_event_peak_RS{id}];
%     sf_use_RS = find(~isnan(mean(sf_peak_RS,1)));
%     errorbarxy(mean(sf_peak_RS(1,sf_use_RS),2), nanmean(sf_peak_RS(2,sf_use_RS),2), std(sf_peak_RS(1,sf_use_RS),[],2)./sqrt(length(sf_use_RS)), std(sf_peak_RS(2,sf_use_RS),[],2)./sqrt(length(sf_use_RS)),{['o' col_mat(id,:)], col_mat(id,:),col_mat(id,:)});
%     hold on
%     scatter(mean(sf_peak_RS(1,sf_use_RS),2), nanmean(sf_peak_RS(2,sf_use_RS),2), col_mat(id,:));
% end
x = 0:.1:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release event peak dF/F')
ylabel('Press event peak dF/F')
title(['All cells'])
subplot(2,2,2)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
title(['All cells'])
subplot(2,2,3)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release event peak dF/F')
ylabel('Press event peak dF/F')
title(['Responsive cells'])
subplot(2,2,4)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
title(['Responsive cells'])
supertitle('Peak event amplitude')
saveas(fig, [out_base 'Summary_event_amp_avg_scatter.fig']);
print([out_base 'Summary_event_amp_avg_scatter.eps'], '-depsc');
print([out_base 'Summary_event_amp_avg_scatter.pdf'], '-dpdf');

%summary of PSTH- average
fig=figure;
for id = 1:size(mouseID,2)
    subplot(4,5,id)
    shadedErrorBar(th{id}, nanmean(success_hist{id},2), std(success_hist{id},[],2)./sqrt(size(success_hist{id},2)), 'k');
    hold on;
    shadedErrorBar(th{id}, nanmean(fail_hist{id},2), std(fail_hist{id},[],2)./sqrt(size(fail_hist{id},2)), 'r');
    title([date{id} ' ' mouseID{id}])
end
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
supertitle(['PSTH- all cells- Black: success; Red: Failure'])
saveas(fig, [out_base 'Summary_PSTH_allcells.fig']);
print([out_base 'Summary_PSTH_allcells.eps'], '-depsc');
print([out_base 'Summary_PSTH_allcells.pdf'], '-dpdf');

fig=figure;
for id = 1:size(mouseID,2)
    subplot(4,5,id)
    shadedErrorBar(th{id}, nanmean(success_hist_RS{id},2), std(success_hist_RS{id},[],2)./sqrt(size(success_hist_RS{id},2)), 'k');
    hold on;
    shadedErrorBar(th{id}, nanmean(fail_hist_RS{id},2), std(fail_hist_RS{id},[],2)./sqrt(size(fail_hist_RS{id},2)), 'r');
    title([date{id} ' ' mouseID{id}])
end
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
supertitle(['PSTH- responsive cells- Black: success; Red: Failure'])
saveas(fig, [out_base 'Summary_PSTH_respcells.fig']);
print([out_base 'Summary_PSTH_respcells.eps'], '-depsc');
print([out_base 'Summary_PSTH_respcells.pdf'], '-dpdf');


%% needs to be updated for different acquisition rates
%PSTH PLOTTED HERE
success_hist = cellfun(@transpose, success_hist, 'UniformOutput', 0);
fail_hist = cellfun(@transpose, fail_hist, 'UniformOutput', 0);
success_hist_RS = cellfun(@transpose, success_hist_RS, 'UniformOutput', 0);
fail_hist_RS = cellfun(@transpose, fail_hist_RS, 'UniformOutput', 0);

frame_size = cell2mat(cellfun(@size, success_hist, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

success_hist_all = interp_frame(success_hist, max_frame);
fail_hist_all    = interp_frame(fail_hist, max_frame);
success_hist_RS_all   = interp_frame(success_hist_RS, max_frame);
fail_hist_RS_all = interp_frame(fail_hist_RS, max_frame);

% tth = [1-(sz1/2):(sz1/2)].*double(min((TC_ifi)));
tth = [1-(size(fail_hist_RS_all,2)/2):(size(fail_hist_RS_all,2)/2)].*double(min((TC_ifi)));
fig=figure;
subplot(2,1,1)
shadedErrorBar(tth, nanmean(success_hist_all,1), std(success_hist_all,[],1)./sqrt(size(success_hist_all,1)), 'k');
hold on;
shadedErrorBar(tth, nanmean(fail_hist_all,1), std(fail_hist_all,[],1)./sqrt(size(fail_hist_all,1)), 'r');
xlim([-1000 1000])
ylim([0 3])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['PSTH- all cells- n = ' num2str(total_cells)])
subplot(2,1,2)
shadedErrorBar(tth, nanmean(success_hist_RS_all,1), std(success_hist_RS_all,[],1)./sqrt(size(success_hist_RS_all,1)), 'k');
hold on;
shadedErrorBar(tth, nanmean(fail_hist_RS_all,1), std(fail_hist_RS_all,[],1)./sqrt(size(fail_hist_RS_all,1)), 'r');
xlim([-1000 1000])
ylim([0 3])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['PSTH- responsive cells- n = ' num2str(total_resp)])
supertitle(['PSTH- responsive cells- Black: success; Red: Failure'])
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz.pdf'], '-dpdf');

% second y axis
fig = figure;
shadedErrorBaryy(tth,nanmean(success_hist_RS_all,1),std(success_hist_RS_all,[],1)./sqrt(size(success_hist_RS_all,1)),'k', ...
    tth,nanmean(success_hist_RS_all,1)*double(min(TC_ifi))/1000,std(success_hist_RS_all*double(min(TC_ifi))/1000,[],1)./sqrt(size(success_hist_RS_all,1)),'k', ...
    tth,nanmean(fail_hist_RS_all,1),std(fail_hist_RS_all,[],1)./sqrt(size(fail_hist_RS_all,1)),'r')
shadedErrorBar(tth, nanmean(success_hist_RS_all,1), std(success_hist_RS_all,[],1)./sqrt(size(success_hist_RS_all,1)), 'k');
hold on;
shadedErrorBar(tth, nanmean(fail_hist_RS_all,1), std(fail_hist_RS_all,[],1)./sqrt(size(fail_hist_RS_all,1)), 'r');
xlim([-1000 1000]); ylim([0 3])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['PSTH- responsive cells- n = ' num2str(total_resp)])
supertitle(['PSTH- responsive cells- Black: success; Red: Failure'])
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_resp.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_resp.pdf'], '-dpdf');

% PSTH align to press
psuccess_hist_RS = cellfun(@transpose, psuccess_hist_RS, 'UniformOutput', 0);
pfail_hist_RS = cellfun(@transpose, pfail_hist_RS, 'UniformOutput', 0);

frame_size = cell2mat(cellfun(@size, psuccess_hist, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

psuccess_hist_RS_all   = interp_frame(psuccess_hist_RS, max_frame);
pfail_hist_RS_all = interp_frame(pfail_hist_RS, max_frame);

% tth = [1-(sz1/2):(sz1/2)].*double(min((TC_ifi)));
tth = [1-(size(pfail_hist_RS_all,2)/2):(size(pfail_hist_RS_all,2)/2)].*double(min((TC_ifi)));
ttp = [1-pre_buffer:size(pfail_hist_RS_all,2) - pre_buffer].*double(min((TC_ifi)));
fig=figure;
shadedErrorBar(ttp, nanmean(psuccess_hist_RS_all,1), nanstd(psuccess_hist_RS_all,[],1)./sqrt(size(psuccess_hist_RS_all,1)), 'k');
hold on;
shadedErrorBar(ttp, nanmean(pfail_hist_RS_all,1), nanstd(pfail_hist_RS_all,[],1)./sqrt(size(pfail_hist_RS_all,1)), 'r');
% xlim([-1000 1000])
ylim([0 3])
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
% title(['PSTH- responsive cells- n = ' num2str(total_resp)])
supertitle(['PSTH- responsive cells- Black: success; Red: Failure'])
saveas(fig, [out_base 'Summary_PSTH_avgexpt_30Hz_press.fig']);
print([out_base 'Summary_PSTH_avgexpt_30Hz_press.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_30Hz_press.pdf'], '-dpdf');

%scatter of peak rate
fig=figure; 
all_sr = [];
all_si = [];
RS_sr = [];
RS_si = [];
all_fr = [];
all_fi = [];
RS_fr = [];
RS_fi = [];
cid = 1;
for id = 1:size(mouseID,2)
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(2,2,1)
%     scatter(success_rate{id}, fail_rate{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
%     subplot(2,2,2)
%     scatter(success_RS_rate{id}, fail_RS_rate{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
%     subplot(2,2,3)
%     scatter(success_ind{id}, fail_ind{id}, 'MarkerEdgeColor', col_mat(cid,:))
%     hold on
%     subplot(2,2,4)
%     scatter(success_RS_ind{id}, fail_RS_ind{id},'MarkerEdgeColor',  col_mat(cid,:))
%     hold on
    all_sr = [all_sr success_rate{id}];
    all_si = [all_si success_ind{id}];
    RS_sr = [RS_sr success_RS_rate{id}];
    RS_si = [RS_si success_RS_ind{id}];
    all_fr = [all_fr fail_rate{id}];
    all_fi = [all_fi fail_ind{id}];
    RS_fr = [RS_fr fail_RS_rate{id}];
    RS_fi = [RS_fi fail_RS_ind{id}];
end
col_mat_s = 0.5*ones(size(col_mat));
    subplot(2,2,1)
    scatter_plot(mouseID, success_rate, fail_rate,  col_mat_s)
    hold on
    subplot(2,2,2)
    scatter_plot(mouseID, success_RS_rate, fail_RS_rate,  col_mat_s)
    hold on
    subplot(2,2,3)
    scatter_plot(mouseID, success_ind, fail_ind, col_mat_s)
    hold on
    subplot(2,2,4)
    scatter_plot(mouseID, success_RS_ind, fail_RS_ind,  col_mat_s)
x = 0:.1:10;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlabel('Correct rate (Hz)')
ylabel('Early rate (Hz)')
xlim([-.5 8])
ylim([-.5 8])
[h_ra p_ra] = ttest(all_sr, all_fr);
title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_ra,2))])
subplot(2,2,2)
plot(x,y,'k')
xlim([-.5 8])
ylim([-.5 8])
xlabel('Correct rate (Hz)')
ylabel('Early rate (Hz)')
[h_rr p_rr] = ttest(RS_sr, RS_fr);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_rr,2))])
subplot(2,2,3)
x = 0:.1:150;
y = x;
plot(x,y,'k');
xlabel('Correct latency (ms)')
ylabel('Early latency (ms)')
[h_la p_la] = ttest(all_si, all_fi);
title(['All cells- n = ' num2str(total_cells) '; p = ' num2str(chop(p_la,2))])
subplot(2,2,4); hold on
plot(x,y,'k')
xlabel('Correct latency (ms)')
ylabel('Early latency (ms)')
[h_lr p_lr] = ttest(RS_si, RS_fi);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_lr,2))])
supertitle(['Event rate and latency on releases- all cells'])
saveas(fig, [out_base 'Summary_rate_latency_scatter.fig']);
print([out_base 'Summary_rate_latency_scatter.eps'], '-depsc');
print([out_base 'Summary_rate_latency_scatter.pdf'], '-dpdf');

fig=figure;
subplot(2,2,1)
scatter_plot(mouseID, success_rate, fail_rate, col_mat);

subplot(2,2,2)
scatter_plot(mouseID, success_RS_rate, fail_RS_rate, col_mat);

subplot(2,2,3)
scatter_plot(mouseID, success_ind, fail_ind, col_mat);

subplot(2,2,4)
scatter_plot(mouseID, success_RS_ind, fail_RS_ind, col_mat);

x = 0:.1:10;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlim([0 6])
ylim([0 6])
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
title(['All cells- rate'])
subplot(2,2,2)
plot(x,y,'k')
xlim([0 6])
ylim([0 6])
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
title(['Responsive cells- rate'])
subplot(2,2,3)
x = 0:1:150;
y = x;
plot(x,y,'k')
xlabel('Success latency (ms)')
ylabel('Fail latency (ms)')
title(['All cells- latency'])
subplot(2,2,4)
plot(x,y,'k')
xlabel('Success latency (ms)')
ylabel('Fail latency (ms)')
title(['Responsive cells- latency'])
supertitle(['Event rate and latency on releases- average within expts'])
saveas(fig, [out_base 'Summary_rate_latency_scatter_avg.fig']);
print([out_base 'Summary_rate_latency_scatter_avg.eps'], '-depsc');
print([out_base 'Summary_rate_latency_scatter_avg.pdf'], '-dpdf');

%==================================================================================================================================
col_mat = [ 0  0  0;
    0 0  0;
    0  0  0;
    0  0  0;
    0  0  0;
    0  0  0;
    0  0  0;
    0  0  0;
    0  0 0 ;
    0  0 0;
    0 0 0; 0.5 0.5 0.5; 0.3 0.5 1; 0.1 0.5 0.7; 0 0.6 0.2;0.8 0.8 0.4;0.1 0.1 0.1;0.3 0.7 0; 0 0 0; 0.1 0.5 0.5];
fig=figure;
subplot(1,2,1)
%scatter_plot(mouseID, success_syn_all_RS, fail_syn_all_RS, col_mat);
scatter_plot3(mouseID, success_syn_all_RS, fail_syn_all_RS);
hold on
x = 0:1:250;
y = x;
plot(x,y,'k'); xlim([0 215]); ylim([0 215]);
axis square
xlabel('Success standard deviation of latency')
ylabel('Fail standard deviation of latency')
title(['All RS cells- st.dev of latency across cells. Averaged across trials'])

subplot(1,2,2)
%scatter_plot(mouseID, success_syn_c_RS, fail_syn_c_RS, col_mat);
scatter_plot3(mouseID, success_syn_c_RS, fail_syn_c_RS);
hold on
x = 0:1:250;
y = x;
plot(x,y,'k'); xlim([0 215]); ylim([0 215]);
axis square
xlabel('Success standard deviation of latency')
ylabel('Fail standard deviation of latency')
title(['All RS cells- st.dev of latency across trials. Avg across cells.'])   

subplot(1,3,3);  
scatter_plot(mouseID, success_syn_c_RS, cue_syn_c_RS, col_mat);
%scatter_plot3(mouseID, success_syn_c_RS, cue_syn_c_RS);
hold on
x = 0:1:250;
y = x;
plot(x,y,'k'); xlim([0 215]); ylim([0 215]);
xlabel('lever aligned standard deviation of latency')
ylabel('cue aligned standard deviation of latency')
title(['All RS cells- st.dev of latency across trials. avg across cells.'])
axis square
supertitle('Syncrony');
saveas(fig, [out_base 'Summary_std_latency_scatter_avg_abs.fig']);
print([out_base 'Summary_std_latency_scatter_avg_abs.eps'], '-depsc');
print([out_base 'Summary_std_latency_scatter_avg_abs.pdf'], '-dpdf');
%==========================================================================================================================================================

%summary of average event waveform relative to spont
fig=figure;
for id = 1:size(mouseID,2)
    subplot(4,5,id)
    release_norm = bsxfun(@rdivide, release_norm_all{id}, max(spont_norm_all{id},[],1));
    press_norm = bsxfun(@rdivide,press_norm_all{id}, max(spont_norm_all{id},[],1));
    spont_norm = bsxfun(@rdivide,spont_norm_all{id}, max(spont_norm_all{id},[],1));
    success_norm = bsxfun(@rdivide,success_norm_all{id}, max(spont_norm_all{id},[],1));
    fail_norm = bsxfun(@rdivide,fail_norm_all{id}, max(spont_norm_all{id},[],1));
    release_norm_peak{id} = max(release_norm,[],1);
    press_norm_peak{id} = max(press_norm,[],1);
    spont_norm_peak{id} = max(spont_norm,[],1);
    success_norm_peak{id} = max(success_norm,[],1);
    fail_norm_peak{id} = max(fail_norm,[],1);
    
    cell_use = find(~isnan(mean([release_norm; press_norm; spont_norm],1)));
    errorbar(ts{id}, nanmean(release_norm(:,cell_use),2), std(release_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
    hold on;
    errorbar(ts{id}, nanmean(press_norm(:,cell_use),2), std(press_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
    hold on;
    errorbar(ts{id}, nanmean(spont_norm(:,cell_use),2), std(spont_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
    title([date{id} ' ' mouseID{id} '-'  num2str(length(cell_use)) ' cells'])   
    xlabel('Time (ms)')
    ylabel('Norm dF/F')
end
supertitle('All cells- Event waveform normalized to spont')
saveas(fig, [out_base 'Summary_avgevent_norm2spont.fig']);
print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');

fig=figure;
for id = 1:size(mouseID,2)
%     subplot(4,5,id)
    normScale = max(spont_norm_RS{id},[],1);
    normScale = 1;
    release_renorm_RS = bsxfun(@rdivide,release_norm_RS{id}, normScale);
    press_renorm_RS = bsxfun(@rdivide,press_norm_RS{id}, normScale);
    spont_renorm_RS = bsxfun(@rdivide,spont_norm_RS{id}, normScale);
    success_renorm_RS = bsxfun(@rdivide,success_norm_RS{id}, normScale);
    fail_renorm_RS = bsxfun(@rdivide,fail_norm_RS{id}, normScale);
    release_norm_peak_RS{id} = max(release_renorm_RS,[],1);
    press_norm_peak_RS{id} = max(press_renorm_RS,[],1);
     spont_norm_peak_RS{id} = max(spont_renorm_RS,[],1);
    success_norm_peak_RS{id} = max(success_renorm_RS,[],1);
    fail_norm_peak_RS{id} = max(fail_renorm_RS,[],1);
    
%     cell_use = find(~isnan(mean([release_renorm_RS; press_renorm_RS; spont_renorm_RS],1)));
%     errorbar(ts{id}, nanmean(release_renorm_RS(:,cell_use),2), std(release_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
%     hold on;
%     errorbar(ts{id}, nanmean(press_renorm_RS(:,cell_use),2), std(press_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
%     hold on;
%     errorbar(ts{id}, nanmean(spont_renorm_RS(:,cell_use),2), std(spont_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
%     title([date{id} ' ' mouseID{id} '-'  num2str(length(cell_use)) ' cells'])   
%     xlabel('Time (ms)')
%     ylabel('Norm dF/F')
end
supertitle('Responsive cells- Event waveform normalized to spont')
saveas(fig, [out_base 'Summary_avgevent_norm2spont.fig']);
print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');

%peak scatters
r_all = [];
s_all = [];  % success
f_all = [];  % fail
sp_all = []; % spont
p_all = [];
r_RS = [];
s_RS = [];
f_RS = [];
sp_RS = [];
p_RS = [];
fig=figure;
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
    r_all = [r_all release_norm_peak{id}];
    p_all = [p_all press_norm_peak{id}];
    r_RS = [r_RS release_norm_peak_RS{id}];
    p_RS = [p_RS press_norm_peak_RS{id}];   
    subplot(1,2,1)
    scatter(release_norm_peak{id}, press_norm_peak{id}, 'MarkerEdgeColor', col_mat(cid,:))
    hold on;
    subplot(1,2,2)
    scatter(release_norm_peak_RS{id}, press_norm_peak_RS{id},'MarkerEdgeColor',  col_mat(cid,:))
    hold on;
end
subplot(1,2,1)
xlim([0 2.5])
ylim([0 2])
hline(1, '--k')
hold on;
vline(1, '--k')
use_cells_all = find(~isnan(mean([p_all; r_all],1)));
[h_rall p_rall] = ttest(r_all(1,use_cells_all),1);
[h_pall p_pall] = ttest(p_all(1,use_cells_all),1);
xlabel('Release amplitude')
ylabel('Press amplitude')
ncells = sum(~isnan(mean([p_all; r_all],1)),2);
title(['All cells- n = ' num2str(total_cells) '; Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
subplot(1,2,2)
xlim([0 2.5])
ylim([0 2])
hline(1, '--k')
hold on;
vline(1, '--k')
use_cells_RS = find(~isnan(mean([p_RS; r_RS],1)),2);
[h_rRS p_rRS] = ttest(r_RS(1,use_cells_RS),1);
[h_pRS p_pRS] = ttest(p_RS(1,use_cells_RS),1);
xlabel('Release amplitude')
ylabel('Press amplitude')
ncells = sum(~isnan(mean([p_RS; r_RS],1)),2);
title(['Responsive cells- n = ' num2str(total_resp) '; Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
supertitle('Event waveform normalized to spont')
saveas(fig, [out_base 'Summary_pr_scatter_norm2spont.fig']);
print([out_base 'Summary_pr_scatter_norm2spont.eps'], '-depsc');
print([out_base 'Summary_pr_scatter_norm2spont.pdf'], '-dpdf');


%peak scatters for spont, fail, success and press
s_all = [];  % success
f_all = [];  % fail
sp_all = []; % spont
p_all = [];
s_RS = [];
f_RS = [];
sp_RS = [];
p_RS = [];
x = [0:.01:1.5];
y = x;
fig=figure;
col_mat_s = [0.5, 0.5, 0.5];
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
%     s_all = [s_all success_norm_peak{id}];
%     f_all = [f_all fail_norm_peak{id}];
%     p_all = [p_all press_norm_peak{id}];
%     sp_all = [sp_all spont_norm_peak{id}];
    s_RS = [s_RS success_norm_peak_RS{id}];
    f_RS = [f_RS fail_norm_peak_RS{id}];
    p_RS = [p_RS press_norm_peak_RS{id}]; 
%     sp_RS = [sp_RS spont_norm_peak_RS{id}];
    
%     subplot(1,3,1)
%     scatter(press_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor', col_mat_s)
%     hold on;
%     subplot(1,3,2)
%     scatter(success_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor',  col_mat_s)
%     hold on;
%     subplot(1,3,3)
%     scatter(fail_norm_peak_RS{id}, spont_norm_peak_RS{id}, 4, 'MarkerEdgeColor',  col_mat_s)
%     hold on;
end
subplot(1,3,1)
hold on; plot(x,y,'-k');
xlim([0 1.5])
ylim([0 1.5])
use_cells_RS = find(~isnan(mean([p_RS; sp_RS],1)));
[h_pspRS p_pspRS] = ttest(sp_RS(1,use_cells_RS),p_RS(1,use_cells_RS));
% [h_pRS p_pRS] = ttest(p_RS(1,use_cells_RS),1);
xlabel('Press amplitude')
ylabel('Spontaneous amplitude')
ncells = sum(~isnan(mean([p_RS; sp_RS],1)),2);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_pspRS,2))])

subplot(1,3,2)
hold on; plot(x,y,'-k');
xlim([0 1.5])
ylim([0 1.5])
use_cells_RS = find(~isnan(mean([s_RS; sp_RS],1)),2);
[h_sRS p_sspRS] = ttest(s_RS(1,use_cells_RS),sp_RS(1,use_cells_RS));
% [h_spRS p_spRS] = ttest(sp_RS(1,use_cells_RS),1);
xlabel('Correct amplitude')
ylabel('Spontaneous amplitude')
ncells = sum(~isnan(mean([s_RS; sp_RS],1)),2);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sspRS,2))])

subplot(1,3,3)
hold on; plot(x,y,'-k');
xlim([0 1.5])
ylim([0 1.5])
use_cells_RS = find(~isnan(mean([f_RS; sp_RS],1)),2);
[h_fRS p_fspRS] = ttest(f_RS(1,use_cells_RS),sp_RS(1,use_cells_RS));
% [h_spRS p_spRS] = ttest(sp_RS(1,use_cells_RS),1);
xlabel('Early amplitude')
ylabel('Spontaneous amplitude')
ncells = sum(~isnan(mean([f_RS; sp_RS],1)),2);
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_fspRS,2))])

supertitle('Peak of event waveform')
% saveas(fig, [out_base 'Summary_event_peak_scatter.fig']);
% print([out_base 'Summary_event_peak_scatter.eps'], '-depsc');
% print([out_base 'Summary_event_peak_scatter.pdf'], '-dpdf');

figure;
% turn on line 619 to scale to max
cid =1;
for id = 1:size(mouseID,2)
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
scatter(success_norm_peak_RS{id}, fail_norm_peak_RS{id}, 4, 'MarkerEdgeColor',  col_mat_s)
hold on;
end

set(gca, 'XScale',' log', 'YScale', 'log');
xlim([0.1 10])
ylim([0.1 10])
hline(1, '--k')
hold on;
vline(1, '--k')
xlabel('Correct amplitude')
ylabel('Early amplitude')
use_cells_RS = find(~isnan(mean([s_RS; f_RS],1)),2);
[h_sfRS p_sfRS] = ttest(s_RS(1,use_cells_RS),f_RS(1,use_cells_RS));
title(['Responsive cells- n = ' num2str(total_resp) '; p = ' num2str(chop(p_sfRS,2))])
saveas(fig, [out_base 'Summary_event_peak_scatter_CE_norm.fig']);
print([out_base 'Summary_event_peak_scatter_CE_norm.eps'], '-depsc');
print([out_base 'Summary_event_peak_scatter_CE_norm.pdf'], '-dpdf');

%average across mice
fig=figure;
cid = 1;
for id = 1:size(mouseID,2) 
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
    subplot(1,2,1)
    temp_release = release_norm_peak{id};
    temp_press = press_norm_peak{id};
    use_cells_all = find(~isnan(mean([temp_release; temp_press],1)));
    temp_release_RS = release_norm_peak_RS{id};
    temp_press_RS = press_norm_peak_RS{id};
    use_cells_RS = find(~isnan(mean([temp_release_RS; temp_press_RS],1)));
    errorbarxy(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), std(temp_release(:,use_cells_all),[],2)./sqrt(length(use_cells_all)), std(temp_press(:,use_cells_all),[],2)./sqrt(length(use_cells_all)),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
    hold on;
    scatter(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), 'MarkerEdgeColor', col_mat(cid,:))
    subplot(1,2,2)
    errorbarxy(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2), std(temp_release_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)), std(temp_press_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
    hold on;
    scatter(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2),'MarkerEdgeColor', col_mat(cid,:))
end
subplot(1,2,1)
xlim([0.5 1.5])
ylim([0.5 1.5])
hline(1, '--k')
hold on;
vline(1, '--k')
xlabel('Release amplitude')
ylabel('Press amplitude')
title('All cells')
subplot(1,2,2)
xlim([0 0.6])
ylim([0 0.6])
hline(1, '--k')
hold on;
vline(1, '--k')
xlabel('Release amplitude')
ylabel('Press amplitude')
title('Responsive cells')
supertitle('Avg event waveform normalized to spont')
saveas(fig, [out_base 'Summary_avgpr_scatter_norm2spont.fig']);
print([out_base 'Summary_avgpr_scatter_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgpr_scatter_norm2spont.pdf'], '-dpdf');


%press and release alone
fig=figure;
lever = strvcat('release', 'press');
% cid = 1;
% for id = 1:size(mouseID,2) 
%     if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
%             
%             cid = cid + 1;
%              
%     end
%     subplot(2,2,1)
%     scatter(ones(size(release_norm_peak{id})), release_norm_peak{id}, 'MarkerEdgeColor',col_mat(cid,:))
%     hold on;
%     scatter(ones(size(press_norm_peak{id})).*2, press_norm_peak{id}, 'MarkerEdgeColor',col_mat(cid,:))
%     hold on;
%     subplot(2,2,2)
%     scatter(ones(size(release_norm_peak_RS{id})), release_norm_peak_RS{id}, 'MarkerEdgeColor',col_mat(cid,:))
%     hold on;
%     scatter(ones(size(press_norm_peak_RS{id})).*2, press_norm_peak_RS{id}, 'MarkerEdgeColor',col_mat(cid,:))
%     hold on;
% end
% subplot(2,2,1)
% xlim([0.5 2.5])
% ylim([0 2.5])
% hline(1, '--k')
% set(gca, 'XTick', 1:2, 'XTickLabels', lever);
% [h_rall p_rall] = ttest(r_all,1);
% [h_pall p_pall] = ttest(p_all,1);
% ylabel('Normalized amplitude')
% title(['Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
% subplot(2,2,2)
% xlim([0.5 2.5])
% ylim([0 2.5])
% hline(1, '--k')
% set(gca, 'XTick', 1:2,'XTickLabels', lever);
% [h_rRS p_rRS] = ttest(r_RS,1);
% [h_pRS p_pRS] = ttest(p_RS,1);
% ylabel('Normalized amplitude')
% title(['Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
cid = 1;
%average across mice
for id = 1:size(mouseID,2) 
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
    end
    subplot(1,2,1)
    temp_release = release_norm_peak{id};
    temp_press = press_norm_peak{id};
    temp_release_RS = release_norm_peak_RS{id};
    temp_press_RS = press_norm_peak_RS{id};
    errorbar(1, nanmean(temp_release,2),nanstd(temp_release,[],2)./sqrt(size(temp_release,2)), 'o', 'color', col_mat(cid,:))
    hold on;
    errorbar(2, nanmean(temp_press,2), nanstd(temp_press,[],2)./sqrt(size(temp_press,2)), 'o', 'color', col_mat(cid,:))
    hold on;
    subplot(1,2,2)
    errorbar(1, nanmean(temp_release_RS,2), nanstd(temp_release_RS,[],2)./sqrt(size(temp_release_RS,2)), 'o', 'color', col_mat(cid,:))
    hold on;
    errorbar(2, nanmean(temp_press_RS,2), nanstd(temp_press_RS,[],2)./sqrt(size(temp_press_RS,2)), 'o', 'color',  col_mat(cid,:))
    hold on;
end
subplot(1,2,1)
xlim([0.5 2.5])
ylim([0 1.5])
hline(1, '--k')
set(gca, 'XTick', 1:2, 'XTickLabel', lever);
ylabel('Normalized amplitude')
title(['All cells- Release n = ' num2str(total_cells) ' cells; Press n = ' num2str(sum(~isnan(p_all),2))])
subplot(1,2,2)
xlim([0.5 2.5])
ylim([0 1.5])
hline(1, '--k')
set(gca, 'XTick', 1:2, 'XTickLabel', lever);
ylabel('Normalized amplitude')
title(['Responsive cells- Release n = ' num2str(total_resp) ' cells; Press n = ' num2str(sum(~isnan(p_RS),2))])
supertitle('Avg event waveform normalized to spont')
saveas(fig, [out_base 'Summary_avg_event_norm2spont.fig']);
print([out_base 'Summary_avg_event_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avg_event_norm2spont.pdf'], '-dpdf');

% %% need to correct for acquisition rate
% %spont waveform comparison
% press_norm_all = cellfun(@transpose, press_norm_all, 'UniformOutput', 0);
% press_norm_RS = cellfun(@transpose, press_norm_RS, 'UniformOutput', 0);
% success_norm_all = cellfun(@transpose, success_norm_all, 'UniformOutput', 0);
% success_norm_RS = cellfun(@transpose, success_norm_RS, 'UniformOutput', 0);
% spont_norm_all = cellfun(@transpose, spont_norm_all, 'UniformOutput', 0);
% spont_norm_RS = cellfun(@transpose, spont_norm_RS, 'UniformOutput', 0);
% 
% frame_size = cell2mat(cellfun(@size, success_norm_all, 'UniformOutput', 0));
% max_frame = max(frame_size(2:2:end));
% 
% all_press_norm = interp_frame(press_norm_all, max_frame);
% RS_press_norm    = interp_frame(press_norm_RS, max_frame);
% all_release_norm   = interp_frame(success_norm_all, max_frame);
% RS_success_norm = interp_frame(success_norm_RS, max_frame);
% all_spont_norm    = interp_frame(spont_norm_all, max_frame);
% RS_spont_norm   = interp_frame(spont_norm_RS, max_frame);
% 
% fig=figure;
% all_spont_renorm = bsxfun(@rdivide, all_spont_norm, max(all_spont_norm,[],2));
% all_success_renorm = bsxfun(@rdivide, all_release_norm, max(all_spont_norm,[],2));
% all_press_renorm = bsxfun(@rdivide, all_press_norm, max(all_spont_norm,[],2));
% RS_spont_renorm = bsxfun(@rdivide, RS_spont_norm, max(RS_spont_norm,[],2));
% RS_success_renorm = bsxfun(@rdivide, RS_success_norm, max(RS_spont_norm,[],2));
% RS_press_renorm = bsxfun(@rdivide, RS_press_norm, max(RS_spont_norm,[],2));
% cell_use = find(~isnan(mean([all_spont_renorm all_success_renorm all_press_renorm],2)));
% tth = [1-(size(all_spont_renorm,2)/2):(size(all_spont_renorm,2)/2)].*double(min((TC_ifi)));
% subplot(1,2,1)
% errorbar(tth, nanmean(all_spont_renorm(cell_use,:),1), std(all_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(all_success_renorm(cell_use,:),1), std(all_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(all_press_renorm(cell_use,:),1), std(all_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- all cells- ' num2str(total_cells) ' cells'])
% 
% subplot(1,2,2)
% cell_use = find(~isnan(mean([RS_spont_renorm RS_success_renorm RS_press_renorm],2)));
% errorbar(tth, nanmean(RS_spont_renorm(cell_use,:),1), std(RS_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(RS_success_renorm(cell_use,:),1), std(RS_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(RS_press_renorm(cell_use,:),1), std(RS_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- responsive cells- ' num2str(total_resp) ' cells'])
% supertitle('Average event - Red: spont; Black: release; Cyan: press')
% saveas(fig, [out_base 'Summary_event_TC_norm2spont.fig']);
% print([out_base 'Summary_event_TC_norm2spont.eps'], '-depsc');
% print([out_base 'Summary_event_TC_norm2spont.pdf'], '-dpdf');


%spont waveform comparison
press_norm_all = cellfun(@transpose, press_norm_all, 'UniformOutput', 0);
press_norm_RS = cellfun(@transpose, press_norm_RS, 'UniformOutput', 0);
success_norm_all = cellfun(@transpose, success_norm_all, 'UniformOutput', 0);
success_norm_RS = cellfun(@transpose, success_norm_RS, 'UniformOutput', 0);
fail_norm_all = cellfun(@transpose, fail_norm_all, 'UniformOutput', 0);
fail_norm_RS = cellfun(@transpose, fail_norm_RS, 'UniformOutput', 0);
spont_norm_all = cellfun(@transpose, spont_norm_all, 'UniformOutput', 0);
spont_norm_RS = cellfun(@transpose, spont_norm_RS, 'UniformOutput', 0);
% lickb_norm_all = cellfun(@transpose, lickb_norm_all, 'UniformOutput', 0);
% lickb_norm_RS = cellfun(@transpose, lickb_norm_RS, 'UniformOutput', 0);
% licks_norm_all = cellfun(@transpose, licks_norm_all, 'UniformOutput', 0);
% licks_norm_RS = cellfun(@transpose, licks_norm_RS, 'UniformOutput', 0);

frame_size = cell2mat(cellfun(@size, success_norm_all, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

all_press_norm = interp_frame(press_norm_all, max_frame);
RS_press_norm    = interp_frame(press_norm_RS, max_frame);
all_success_norm   = interp_frame(success_norm_all, max_frame);
RS_success_norm = interp_frame(success_norm_RS, max_frame);
all_fail_norm   = interp_frame(fail_norm_all, max_frame);
RS_fail_norm = interp_frame(fail_norm_RS, max_frame);
all_spont_norm    = interp_frame(spont_norm_all, max_frame);
RS_spont_norm   = interp_frame(spont_norm_RS, max_frame);

fig=figure;
normScale = max(all_spont_norm,[],2);
normScale = 1;
all_spont_renorm = bsxfun(@rdivide, all_spont_norm, normScale);
all_success_renorm = bsxfun(@rdivide, all_success_norm, normScale);
all_fail_renorm = bsxfun(@rdivide, all_fail_norm, normScale);
all_press_renorm = bsxfun(@rdivide, all_press_norm, normScale);
RS_spont_renorm = bsxfun(@rdivide, RS_spont_norm, normScale);
RS_success_renorm = bsxfun(@rdivide, RS_success_norm, normScale);
RS_fail_renorm = bsxfun(@rdivide, RS_fail_norm, normScale);
RS_press_renorm = bsxfun(@rdivide, RS_press_norm, normScale);
cell_use = find(~isnan(mean([all_spont_renorm all_success_renorm all_fail_renorm all_press_renorm],2)));
tth = [1-(size(all_spont_renorm,2)/2):(size(all_spont_renorm,2)/2)].*double(min((TC_ifi)));

subplot(1,2,1); 
errorbar(tth, nanmean(all_spont_renorm(cell_use,:),1), nanstd(all_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-m')
hold on;
errorbar(tth, nanmean(all_success_renorm(cell_use,:),1), nanstd(all_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tth, nanmean(all_press_renorm(cell_use,:),1), nanstd(all_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
errorbar(tth, nanmean(all_fail_renorm(cell_use,:),1), nanstd(all_fail_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
xlim([ -500 500]);axis square
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- all cells- ' num2str(total_cells) ' cells'])

subplot(1,2,2);
cell_use = find(~isnan(mean([RS_spont_renorm RS_success_renorm RS_fail_renorm RS_press_renorm],2)));
errorbar(tth, nanmean(RS_spont_renorm(cell_use,:),1), nanstd(RS_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-m')
hold on;
errorbar(tth, nanmean(RS_success_renorm(cell_use,:),1), nanstd(RS_success_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tth, nanmean(RS_press_renorm(cell_use,:),1), nanstd(RS_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
errorbar(tth, nanmean(RS_fail_renorm(cell_use,:),1), nanstd(RS_fail_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
xlim([ -500 500]); ylim([-0.5 1.5]); axis square
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- responsive cells- ' num2str(total_resp) ' cells'])
suptitle('Average event- Black: success; Red: early; Cyan: press; Magenta: spont')
saveas(fig, [out_base 'Summary_event_TC_exptavg.fig']);
print([out_base 'Summary_event_TC_exptavg.eps'], '-depsc');
print([out_base 'Summary_event_TC_exptavg.pdf'], '-dpdf');

% %% plot events amplitude timecourse
% all_success_event = cellfun(@transpose, success_event_TC, 'UniformOutput', 0);
% all_fail_event = cellfun(@transpose, fail_event_TC, 'UniformOutput', 0);
% all_press_event = cellfun(@transpose, press_event_TC, 'UniformOutput', 0);
% all_spont_event = cellfun(@transpose, spont_event_TC, 'UniformOutput', 0);
% RS_success_event = cellfun(@transpose, success_event_TC_RS, 'UniformOutput', 0);
% RS_fail_event = cellfun(@transpose, fail_event_TC_RS, 'UniformOutput', 0);
% RS_press_event = cellfun(@transpose, press_event_TC_RS, 'UniformOutput', 0);
% RS_spont_event = cellfun(@transpose, spont_event_TC_RS, 'UniformOutput', 0);
% 
% frame_size = cell2mat(cellfun(@size, success_event_TC, 'UniformOutput', 0));
% max_frame = max(frame_size(2:2:end));
% 
% all_success_event = interp_frame(all_success_event, max_frame);
% all_fail_event    = interp_frame(all_fail_event, max_frame);
% all_spont_event    = interp_frame(all_spont_event, max_frame);
% all_press_event    = interp_frame(all_press_event, max_frame);
% RS_success_event   = interp_frame(RS_success_event, max_frame);
% RS_fail_event = interp_frame(RS_fail_event, max_frame);
% RS_spont_event = interp_frame(RS_spont_event, max_frame);
% RS_press_event = interp_frame(RS_press_event, max_frame);
% 
% cell_use = find(~isnan(mean([all_success_event all_fail_event all_spont_event all_press_event],1)));
% 
% % tth = [1-(sz1/2):(sz1/2)].*double(min((TC_ifi)));
% tth = [1-(size(all_success_event,2)/2):(size(all_success_event,2)/2)].*double(min((TC_ifi)));
% fig=figure;
% subplot(2,1,1)
% errorbar(tth, nanmean(all_success_event(cell_use,:),1), std(all_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(all_fail_event(cell_use,:),1), std(all_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(all_press_event(cell_use,:),1), std(all_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% errorbar(tth, nanmean(all_spont_event(cell_use,:),1), std(all_spont_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-m')
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- all cells- ' num2str(length(cell_use)) ' cells'])
% 
% subplot(1,2,2)
% cell_use = find(~isnan(mean([RS_success_event RS_fail_event RS_press_event],2)));
% errorbar(tth, nanmean(RS_success_event(cell_use,:),1), std(RS_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
% hold on;
% errorbar(tth, nanmean(RS_fail_event(cell_use,:),1), std(RS_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
% hold on;
% errorbar(tth, nanmean(RS_press_event(cell_use,:),1), std(RS_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% hold on;
% errorbar(tth, nanmean(RS_spont_event(cell_use,:),1), std(RS_spont_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
% xlabel('Time (ms)')
% ylabel('dF/F')
% title(['Average event- responsive cells- ' num2str(length(cell_use)) ' cells'])
% suptitle('Average event- Black: success; Red: early; Cyan: press; Megneta: spont')
% print([out_base 'Summary_event_TC_exptavg.eps'], '-depsc');
% print([out_base 'Summary_event_TC_exptavg.pdf'], '-dpdf');
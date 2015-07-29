clear all
Jake_2P_exptlist
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    run = run_mat(id,:,1);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    nrun = nrun_mat(id,:);
    if nrun == 1
        run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    else
        run_name = [date '_' mouse '_run' run(length(run)-2:end) '-00' num2str(nrun-1)];
    end
    
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);
    
    dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
    load([dest '_parse_behavior']);
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub '_event_summary.mat']);
    load([dest_sub '_event_peaks.mat']);
    load([dest_sub '_spont_events.mat']);
    load([dest_sub '_event_hist.mat']);
    
    ncells(id) = size(press,2);
    RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells]);
    
    TC_length(id) = size(all_success_event,2);
    TC_ifi(id) = ifi;

    press_event_TC{id} = all_press_event;
    success_event_TC{id} = all_success_event;
    fail_event_TC{id} = all_fail_event;
    press_event_peak{id} = all_press_peak;
    release_event_peak{id} = all_release_peak;
    success_event_peak{id} = all_success_peak;
    fail_event_peak{id} = all_fail_peak;
    
    press_event_TC_RS{id} = all_press_event(RS_cells{id},:);
    success_event_TC_RS{id} = all_success_event(RS_cells{id},:);
    fail_event_TC_RS{id} = all_fail_event(RS_cells{id},:);
    press_event_peak_RS{id} = all_press_peak(:,RS_cells{id});
    release_event_peak_RS{id} = all_release_peak(:,RS_cells{id});
    success_event_peak_RS{id} = all_success_peak(:,RS_cells{id});
    fail_event_peak_RS{id} = all_fail_peak(:,RS_cells{id});
    
    success_hist{id} = all_success_hist.*(1000/double(ifi));
    fail_hist{id} = all_fail_hist.*(1000/double(ifi));
    success_hist_RS{id} = resp_success_hist.*(1000/double(ifi));
    fail_hist_RS{id} = resp_fail_hist.*(1000/double(ifi));
    
    success_rate{id} = success_rate_all;
    fail_rate{id} = fail_rate_all;
    success_RS_rate{id} = success_rate_RS;
    fail_RS_rate{id} = fail_rate_RS;
    
    success_ind{id} = (success_ind_all-1).*double(ifi);
    fail_ind{id} = (fail_ind_all-1).*double(ifi);
    success_RS_ind{id} = (success_ind_RS-1).*double(ifi);
    fail_RS_ind{id} = (fail_ind_RS-1).*double(ifi);
    
    rate{id} = events_rate;
    
    ts{id} = [-data_start:data_end].*double(ifi);
    sz1 = size(all_success_hist,1);
    th{id} = [1-(sz1/2):(sz1/2)].*double(ifi);
end
col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm');
%summary of event rate
figure
rate_all = [];
for id = 1:length(date_mat)
    plot(id*ones(size(rate{id})), rate{id}, ['o' col_mat(id,:)]);
    rate_all = [rate_all rate{id}];
    hold on
end
ylim([0 2])
xlim([0.5 length(date_mat)+.5])
ylabel('Spike rate (Hz)')
xlabel('session')
title(['Avg rate: ' num2str(chop(mean(rate_all,2),2)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
print([out_base 'Summary_event_rate.eps'], '-depsc');
print([out_base 'Summary_event_rate.pdf'], '-dpdf');

%summary of average events
figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    errorbar(ts{id}, nanmean(press_event_TC{id},1), nanstd(press_event_TC{id},[],1)./sqrt(ncells(id)), 'c');
    hold on
    errorbar(ts{id}, nanmean(success_event_TC{id},1), nanstd(success_event_TC{id},[],1)./sqrt(ncells(id)), 'k');
    hold on
    errorbar(ts{id}, nanmean(fail_event_TC{id},1), nanstd(fail_event_TC{id},[],1)./sqrt(ncells(id)), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:) '- ' num2str(ncells(id)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
suptitle(['Good events- All cells- Black: success; Red: early; Cyan: press'])
print([out_base 'Summary_allcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_allcells_event_TC.pdf'], '-dpdf');
figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    errorbar(ts{id}, nanmean(press_event_TC_RS{id},1), nanstd(press_event_TC_RS{id},[],1)./sqrt(length(RS_cells{id})), 'c');
    hold on
    errorbar(ts{id}, nanmean(success_event_TC_RS{id},1), nanstd(success_event_TC_RS{id},[],1)./sqrt(length(RS_cells{id})), 'k');
    hold on
    errorbar(ts{id}, nanmean(fail_event_TC_RS{id},1), nanstd(fail_event_TC_RS{id},[],1)./sqrt(length(RS_cells{id})), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:) '- ' num2str(ncells(id)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
suptitle(['Good events- Resp cells- Black: success; Red: early; Cyan: press'])
print([out_base 'Summary_respcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_respcells_event_TC.pdf'], '-dpdf');

all_press_event = [];
RS_press_event = [];
all_fail_event = [];
RS_fail_event = [];
all_success_event = [];
RS_success_event = [];
for id = 1:4
    all_press_event = [all_press_event; press_event_TC{id}];
    RS_press_event = [RS_press_event; press_event_TC_RS{id}];
    all_fail_event = [all_fail_event; fail_event_TC{id}];
    RS_fail_event = [RS_fail_event; fail_event_TC_RS{id}];
    all_success_event = [all_success_event; success_event_TC{id}];
    RS_success_event = [RS_success_event; success_event_TC_RS{id}];
end
figure;
subplot(1,2,1)
errorbar(ts{1}, nanmean(all_success_event,1), nanstd(all_success_event,[],1)./sqrt(sum(~isnan(all_success_event(:,1)),1)), '-k')
hold on;
errorbar(ts{1}, nanmean(all_fail_event,1), nanstd(all_fail_event,[],1)./sqrt(sum(~isnan(all_fail_event(:,1)),1)), '-r')
hold on;
errorbar(ts{1}, nanmean(all_press_event,1), nanstd(all_press_event,[],1)./sqrt(sum(~isnan(all_press_event(:,1)),1)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title('Average event- all cells')

subplot(1,2,2)
errorbar(ts{1}, nanmean(RS_success_event,1), nanstd(RS_success_event,[],1)./sqrt(sum(~isnan(RS_success_event(:,1)),1)), '-k')
hold on;
errorbar(ts{1}, nanmean(RS_fail_event,1), nanstd(RS_fail_event,[],1)./sqrt(sum(~isnan(RS_fail_event(:,1)),1)), '-r')
hold on;
errorbar(ts{1}, nanmean(RS_press_event,1), nanstd(RS_press_event,[],1)./sqrt(sum(~isnan(RS_press_event(:,1)),1)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title('Average event- responsive cells')
suptitle('Average event - 15Hz expts - Black: success; Red: early; Cyan: press')
print([out_base 'Summary_event_TC_15Hz_exptavg.eps'], '-depsc');
print([out_base 'Summary_event_TC_15Hz_exptavg.pdf'], '-dpdf');

%summary of average event peaks
figure;
ra = [];
pa = [];
rr = [];
pr = [];
sa = [];
fa = [];
sr = [];
fr = [];
for id = 1:6
    ra = [ra release_event_peak{id}];
    pa = [pa press_event_peak{id}];
    rr = [rr release_event_peak_RS{id}];
    pr = [pr press_event_peak_RS{id}];
    sa = [sa success_event_peak{id}];
    fa = [fa fail_event_peak{id}];
    sr = [sr success_event_peak_RS{id}];
    fr = [fr fail_event_peak_RS{id}];
    subplot(2,2,1)
    scatter(release_event_peak{id}, press_event_peak{id}, col_mat(id,:));
    hold on
    subplot(2,2,2)
    scatter(success_event_peak{id}, fail_event_peak{id}, col_mat(id,:));
    hold on
    subplot(2,2,3)
    scatter(release_event_peak_RS{id}, press_event_peak_RS{id}, col_mat(id,:));
    hold on
    subplot(2,2,4)
    scatter(success_event_peak_RS{id}, fail_event_peak_RS{id}, col_mat(id,:));
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
title(['All cells- n = ' num2str(sum(~isnan(ra+pa),2)) '; p = ' num2str(chop(p_rpa,2))])
subplot(2,2,2)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
[h_sfa p_sfa]= ttest(sa,fa);
title(['All cells- n = ' num2str(sum(~isnan(sa+fa),2)) '; p = ' num2str(chop(p_sfa,2))])
subplot(2,2,3)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release event peak dF/F')
ylabel('Press event peak dF/F')
[h_rpr p_rpr]= ttest(rr,pr);
title(['Responsive cells- n = ' num2str(sum(~isnan(rr+pr),2)) '; p = ' num2str(chop(p_rpr,2))])
subplot(2,2,4)
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success event peak dF/F')
ylabel('Fail event peak dF/F') 
[h_sfr p_sfr]= ttest(sr,fr);
title(['Responsive cells- n = ' num2str(sum(~isnan(sr+fr),2)) '; p = ' num2str(chop(p_sfr,2))])
suptitle('Peak event amplitude- img24 and img25')
print([out_base 'Summary_event_amp.eps'], '-depsc');
print([out_base 'Summary_event_amp.pdf'], '-dpdf');
% suptitle('Peak event amplitude- img24 and img25')
% print([out_base 'Summary_event_amp_img24_25.eps'], '-depsc');
% print([out_base 'Summary_event_amp_img24_25.pdf'], '-dpdf');

%summary of PSTH- average
figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    shadedErrorBar(th{id}, mean(success_hist{id},2), std(success_hist{id},[],2)./sqrt(size(success_hist{id},2)), 'k');
    hold on;
    shadedErrorBar(th{id}, mean(fail_hist{id},2), std(fail_hist{id},[],2)./sqrt(size(fail_hist{id},2)), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:)])
end
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
suptitle(['PSTH- all cells- Black: success; Red: Failure'])
print([out_base 'Summary_PSTH_allcells.eps'], '-depsc');
print([out_base 'Summary_PSTH_allcells.pdf'], '-dpdf');

figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    shadedErrorBar(th{id}, mean(success_hist_RS{id},2), std(success_hist_RS{id},[],2)./sqrt(size(success_hist_RS{id},2)), 'k');
    hold on;
    shadedErrorBar(th{id}, mean(fail_hist_RS{id},2), std(fail_hist_RS{id},[],2)./sqrt(size(fail_hist_RS{id},2)), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:)])
end
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
suptitle(['PSTH- responsive cells- Black: success; Red: Failure'])
print([out_base 'Summary_PSTH_respcells.eps'], '-depsc');
print([out_base 'Summary_PSTH_respcells.pdf'], '-dpdf');

success_hist_all = [];
fail_hist_all = [];
success_hist_RS_all = [];
fail_hist_RS_all = [];
for id = 1:4
    success_hist_all = [success_hist_all success_hist{id}];
    fail_hist_all = [fail_hist_all fail_hist{id}];
    success_hist_RS_all = [success_hist_RS_all success_hist_RS{id}];
    fail_hist_RS_all = [fail_hist_RS_all fail_hist_RS{id}];
end
figure;
subplot(2,1,1)
shadedErrorBar(th{1}, mean(success_hist_all,2), std(success_hist_all,[],2)./sqrt(size(success_hist_all,2)), 'k');
hold on;
shadedErrorBar(th{1}, mean(fail_hist_all,2), std(fail_hist_all,[],2)./sqrt(size(fail_hist_all,2)), 'r');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['PSTH- all cells- n = ' num2str(size(fail_hist_all,2))])
subplot(2,1,2)
shadedErrorBar(th{1}, mean(success_hist_RS_all,2), std(success_hist_RS_all,[],2)./sqrt(size(success_hist_RS_all,2)), 'k');
hold on;
shadedErrorBar(th{1}, mean(fail_hist_RS_all,2), std(fail_hist_RS_all,[],2)./sqrt(size(fail_hist_RS_all,2)), 'r');
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['PSTH- responsive cells- n = ' num2str(size(fail_hist_RS_all,2))])
suptitle(['PSTH- responsive cells- Black: success; Red: Failure'])
print([out_base 'Summary_PSTH_avgexpt_15Hz.eps'], '-depsc');
print([out_base 'Summary_PSTH_avgexpt_15Hz.pdf'], '-dpdf');

%scatter of peak rate
figure; 
all_sr = [];
all_si = [];
RS_sr = [];
RS_si = [];
all_fr = [];
all_fi = [];
RS_fr = [];
RS_fi = [];
for id = 1:length(date_mat)
    subplot(2,2,1)
    scatter(success_rate{id}, fail_rate{id}, col_mat(id,:))
    hold on
    subplot(2,2,2)
    scatter(success_RS_rate{id}, fail_RS_rate{id}, col_mat(id,:))
    hold on
    subplot(2,2,3)
    scatter(success_ind{id}, fail_ind{id}, col_mat(id,:))
    hold on
    subplot(2,2,4)
    scatter(success_RS_ind{id}, fail_RS_ind{id}, col_mat(id,:))
    hold on
    all_sr = [all_sr success_rate{id}];
    all_si = [all_si success_ind{id}];
    RS_sr = [RS_sr success_RS_rate{id}];
    RS_si = [RS_si success_RS_ind{id}];
    all_fr = [all_fr fail_rate{id}];
    all_fi = [all_fi fail_ind{id}];
    RS_fr = [RS_fr fail_RS_rate{id}];
    RS_fi = [RS_fi fail_RS_ind{id}];
end
x = 0:.1:.5;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
xlim([-.05 0.3])
ylim([-.05 0.3])
[h_ra p_ra] = ttest(all_sr, all_fr);
title(['All cells- n = ' num2str(sum(~isnan(all_sr+all_fr),2)) '; p = ' num2str(chop(p_ra,2))])
subplot(2,2,2)
plot(x,y,'k')
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
xlim([-.05 0.3])
ylim([-.05 0.3])
[h_rr p_rr] = ttest(RS_sr, RS_fr);
title(['Responsive cells- n = ' num2str(sum(~isnan(RS_sr+RS_fr),2)) '; p = ' num2str(chop(p_rr,2))])
subplot(2,2,3)
xlabel('Success latency (ms)')
ylabel('Fail latency (ms)')
[h_la p_la] = ttest(all_si, all_fi);
title(['All cells- n = ' num2str(sum(~isnan(all_si+all_fi),2)) '; p = ' num2str(chop(p_la,2))])
subplot(2,2,4)
plot(x,y,'k')
xlabel('Success latency (ms)')
ylabel('Fail latency (ms)')
[h_lr p_lr] = ttest(RS_si, RS_fi);
title(['Responsive cells- n = ' num2str(sum(~isnan(RS_si+RS_fi),2)) '; p = ' num2str(chop(p_lr,2))])
suptitle(['Event rate and latency on releases- all cells'])
print([out_base 'Summary_rate_latency_scatter.eps'], '-depsc');
print([out_base 'Summary_rate_latency_scatter.pdf'], '-dpdf');

for id = 1:length(date_mat)
    subplot(2,2,1)
    scatter(mean(success_rate{id},2), mean(fail_rate{id},2), col_mat(id,:))
    hold on
    subplot(2,2,2)
    scatter(mean(success_RS_rate{id},2), mean(fail_RS_rate{id},2), col_mat(id,:))
    hold on
    subplot(2,2,3)
    scatter(mean(success_ind{id},2), mean(fail_ind{id},2), col_mat(id,:))
    hold on
    subplot(2,2,4)
    scatter(mean(success_RS_ind{id},2), mean(fail_RS_ind{id},2), col_mat(id,:))
    hold on
end
x = 0:.1:.5;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
xlim([-.05 0.3])
ylim([-.05 0.3])
title(['All cells- rate'])
subplot(2,2,2)
plot(x,y,'k')
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
xlim([-.05 0.3])
ylim([-.05 0.3])
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
suptitle(['Event rate and latency on releases- average within expts'])
print([out_base 'Summary_rate_latency_scatter_avg.eps'], '-depsc');
print([out_base 'Summary_rate_latency_scatter_avg.pdf'], '-dpdf');
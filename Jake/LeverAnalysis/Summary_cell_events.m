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
    load([dest_sub '_norm2spont.mat']);
    
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
    
    release_norm_RS{id} = release_norm_avg(:,RS_cells{id});
    press_norm_RS{id} = press_norm_avg(:,RS_cells{id});
    spont_norm_RS{id} = spont_norm_avg(:,RS_cells{id});
    
    rate{id} = events_rate;
    
    ts{id} = [-data_start:data_end].*double(ifi);
    sz1 = size(all_success_hist,1);
    th{id} = [1-(sz1/2):(sz1/2)].*double(ifi);
end
%%plotting
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
title(['Avg rate: ' num2str(chop(mean(rate_all,2),3)) ' +/-' num2str(chop(std(rate_all,[],2)./(sqrt(size(rate_all,2))),2)) ' n = ' num2str(size(rate_all,2)) ' cells'])
print([out_base 'Summary_event_rate.eps'], '-depsc');
print([out_base 'Summary_event_rate.pdf'], '-dpdf');

%summary of average events
figure;
for id = 1:length(date_mat)
    press_events = press_event_TC{id};
    success_events = success_event_TC{id};
    fail_events = fail_event_TC{id};
    cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
    subplot(2,3,id)
    errorbar(ts{id}, mean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
    hold on
    errorbar(ts{id}, mean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
    hold on
    errorbar(ts{id}, mean(fail_events(cell_use,:),1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:) '- ' num2str(length(cell_use)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
suptitle(['Good events- All cells- Black: success; Red: early; Cyan: press'])
print([out_base 'Summary_allcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_allcells_event_TC.pdf'], '-dpdf');
figure;
for id = 1:length(date_mat)
    press_events = press_event_TC_RS{id};
    success_events = success_event_TC_RS{id};
    fail_events = fail_event_TC_RS{id};
    cell_use = find(~isnan(mean([press_events success_events fail_events],2)));
    subplot(2,3,id)
    errorbar(ts{id}, mean(press_events(cell_use,:),1), std(press_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'c');
    hold on
    errorbar(ts{id}, mean(success_events(cell_use,:),1), std(success_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'k');
    hold on
    errorbar(ts{id}, mean(fail_events(cell_use,:),1), std(fail_events(cell_use,:),[],1)./sqrt(length(cell_use)), 'r');
    title([date_mat(id,:) ' ' mouse_mat(id,:) '- ' num2str(length(cell_use)) ' cells'])
end
xlabel('Time (ms)')
ylabel('dF/F')
suptitle(['Good events- Resp cells- Black: success; Red: early; Cyan: press'])
print([out_base 'Summary_respcells_event_TC.eps'], '-depsc');
print([out_base 'Summary_respcells_event_TC.pdf'], '-dpdf');

%average across days- only including 15 and 30Hz animals (img24,25,27)
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
for id = 6
    press_TC_temp = press_event_TC{id};
    sz = size(press_TC_temp);
    all_press_event = [all_press_event(:,2:16); squeeze(mean(reshape(press_TC_temp(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    press_TC_temp_RS = press_event_TC_RS{id};
    sz = size(press_TC_temp_RS);
    RS_press_event = [RS_press_event(:,2:16); squeeze(mean(reshape(press_TC_temp_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    success_TC_temp = success_event_TC{id};
    sz = size(success_TC_temp);
    all_success_event = [all_success_event(:,2:16); squeeze(mean(reshape(success_TC_temp(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    success_TC_temp_RS = success_event_TC_RS{id};
    sz = size(success_TC_temp_RS);
    RS_success_event = [RS_success_event(:,2:16); squeeze(mean(reshape(success_TC_temp_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    fail_TC_temp = fail_event_TC{id};
    sz = size(fail_TC_temp);
    all_fail_event = [all_fail_event(:,2:16); squeeze(mean(reshape(fail_TC_temp(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    fail_TC_temp_RS = fail_event_TC_RS{id};
    sz = size(fail_TC_temp_RS);
    RS_fail_event = [RS_fail_event(:,2:16); squeeze(mean(reshape(fail_TC_temp_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
end
figure;
cell_use = find(~isnan(mean([all_success_event all_fail_event all_press_event],2)));
tts = ts{1};
subplot(1,2,1)
errorbar(tts(1,2:16), mean(all_success_event(cell_use,:),1), std(all_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tts(1,2:16), mean(all_fail_event(cell_use,:),1), std(all_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
hold on;
errorbar(tts(1,2:16), mean(all_press_event(cell_use,:),1), std(all_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- all cells- ' num2str(length(cell_use)) ' cells'])

subplot(1,2,2)
cell_use = find(~isnan(mean([RS_success_event RS_fail_event RS_press_event],2)));
errorbar(tts(1,2:16), mean(RS_success_event(cell_use,:),1), std(RS_success_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tts(1,2:16), mean(RS_fail_event(cell_use,:),1), std(RS_fail_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-r')
hold on;
errorbar(tts(1,2:16), mean(RS_press_event(cell_use,:),1), std(RS_press_event(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- responsive cells- ' num2str(length(cell_use)) ' cells'])
suptitle('Average event - img24, img25 and img27 - Black: success; Red: early; Cyan: press')
print([out_base 'Summary_event_TC_exptavg.eps'], '-depsc');
print([out_base 'Summary_event_TC_exptavg.pdf'], '-dpdf');

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
for id = [1:4 6]
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
suptitle('Peak event amplitude')
print([out_base 'Summary_event_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_event_amp_scatter.pdf'], '-dpdf');

%average event amplitude by mouse
figure;
for id = [1:4 6]
    subplot(2,2,1)
    rp_peak = [release_event_peak{id}; press_event_peak{id}];
    rp_use = find(~isnan(mean(rp_peak,1)));
    errorbarxy(mean(rp_peak(1,rp_use),2), mean(rp_peak(2,rp_use),2), std(rp_peak(1,rp_use),[],2)./sqrt(length(rp_use)), std(rp_peak(2,rp_use),[],2)./sqrt(length(rp_use)),[],[],[], col_mat(id,:));
    hold on
    scatter(mean(rp_peak(1,rp_use),2), mean(rp_peak(2,rp_use),2),col_mat(id,:));
    subplot(2,2,2)
    sf_peak = [success_event_peak{id}; fail_event_peak{id}];
    sf_use = find(~isnan(mean(sf_peak,1)));
    errorbarxy(mean(sf_peak(1,sf_use),2), mean(sf_peak(2,sf_use),2),std(sf_peak(1,sf_use),[],2)./sqrt(length(sf_use)), std(sf_peak(2,sf_use),[],2)./sqrt(length(sf_use)),[],[],[], col_mat(id,:));
    hold on
    scatter(mean(sf_peak(1,sf_use),2), mean(sf_peak(2,sf_use),2),col_mat(id,:));
    subplot(2,2,3)
    rp_peak_RS = [release_event_peak_RS{id}; press_event_peak_RS{id}];
    rp_use_RS = find(~isnan(mean(rp_peak_RS,1)));
    errorbarxy(mean(rp_peak_RS(1,rp_use_RS),2), mean(rp_peak_RS(2,rp_use_RS),2), std(rp_peak_RS(1,rp_use_RS),[],2)./sqrt(length(rp_use_RS)), std(rp_peak_RS(2,rp_use_RS),[],2)./sqrt(length(rp_use_RS)),[],[],[], col_mat(id,:));
    hold on
    scatter(mean(rp_peak_RS(1,rp_use_RS),2), mean(rp_peak_RS(2,rp_use_RS),2), col_mat(id,:));
    subplot(2,2,4)
    sf_peak_RS = [success_event_peak_RS{id}; fail_event_peak_RS{id}];
    sf_use_RS = find(~isnan(mean(sf_peak_RS,1)));
    errorbarxy(mean(sf_peak_RS(1,sf_use_RS),2), mean(sf_peak_RS(2,sf_use_RS),2), std(sf_peak_RS(1,sf_use_RS),[],2)./sqrt(length(sf_use_RS)), std(sf_peak_RS(2,sf_use_RS),[],2)./sqrt(length(sf_use_RS)),[],[],[], col_mat(id,:));
    hold on
    scatter(mean(sf_peak_RS(1,sf_use_RS),2), mean(sf_peak_RS(2,sf_use_RS),2), col_mat(id,:));
end
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
suptitle('Peak event amplitude')
print([out_base 'Summary_event_amp_avg_scatter.eps'], '-depsc');
print([out_base 'Summary_event_amp_avg_scatter.pdf'], '-dpdf');

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
x = 0:.1:10;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
xlim([-.5 8])
ylim([-.5 8])
[h_ra p_ra] = ttest(all_sr, all_fr);
title(['All cells- n = ' num2str(sum(~isnan(all_sr+all_fr),2)) '; p = ' num2str(chop(p_ra,2))])
subplot(2,2,2)
plot(x,y,'k')
xlim([-.5 8])
ylim([-.5 8])
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
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

figure;
for id = 1:length(date_mat)
    subplot(2,2,1)
    errorbarxy(mean(success_rate{id},2), mean(fail_rate{id},2), std(success_rate{id},[],2)./sqrt(size(success_rate{id},2)),std(fail_rate{id},[],2)./sqrt(size(fail_rate{id},2)),[],[],[], col_mat(id,:))
    hold on
    scatter(mean(success_rate{id},2), mean(fail_rate{id},2), col_mat(id,:))    
    subplot(2,2,2)
    errorbarxy(mean(success_RS_rate{id},2), mean(fail_RS_rate{id},2), std(success_RS_rate{id},[],2)./sqrt(size(success_RS_rate{id},2)),std(fail_RS_rate{id},[],2)./sqrt(size(fail_RS_rate{id},2)),[],[],[], col_mat(id,:))
    hold on
    scatter(mean(success_RS_rate{id},2), mean(fail_RS_rate{id},2), col_mat(id,:))
    subplot(2,2,3)
    errorbarxy(mean(success_ind{id},2), mean(fail_ind{id},2), std(success_ind{id},[],2)./sqrt(size(success_ind{id},2)),std(fail_ind{id},[],2)./sqrt(size(fail_ind{id},2)),[],[],[], col_mat(id,:))
    hold on
    scatter(mean(success_ind{id},2), mean(fail_ind{id},2), col_mat(id,:))
    subplot(2,2,4)
    errorbarxy(mean(success_RS_ind{id},2), mean(fail_RS_ind{id},2), std(success_RS_ind{id},[],2)./sqrt(size(success_RS_ind{id},2)),std(fail_RS_ind{id},[],2)./sqrt(size(fail_RS_ind{id},2)),[],[],[], col_mat(id,:))
    hold on
    scatter(mean(success_RS_ind{id},2), mean(fail_RS_ind{id},2), col_mat(id,:))
end
x = 0:.1:10;
y = x;
subplot(2,2,1)
plot(x,y,'k')
xlim([-.5 3])
ylim([-.5 3])
xlabel('Success rate (Hz)')
ylabel('Fail rate (Hz)')
title(['All cells- rate'])
subplot(2,2,2)
plot(x,y,'k')
xlim([-.5 3])
ylim([-.5 3])
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
suptitle(['Event rate and latency on releases- average within expts'])
print([out_base 'Summary_rate_latency_scatter_avg.eps'], '-depsc');
print([out_base 'Summary_rate_latency_scatter_avg.pdf'], '-dpdf');

%summary of average event waveform relative to spont
figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    release_norm = bsxfun(@rdivide, release_norm_all{id}, max(spont_norm_all{id},[],1));
    press_norm = bsxfun(@rdivide,press_norm_all{id}, max(spont_norm_all{id},[],1));
    spont_norm = bsxfun(@rdivide,spont_norm_all{id}, max(spont_norm_all{id},[],1));
    release_norm_peak{id} = max(release_norm,[],1);
    press_norm_peak{id} = max(press_norm,[],1);
    spont_norm_peak{id} = max(spont_norm,[],1);
    cell_use = find(~isnan(mean([release_norm; press_norm; spont_norm],1)));
    errorbar(ts{id}, mean(release_norm(:,cell_use),2), std(release_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
    hold on;
    errorbar(ts{id}, mean(press_norm(:,cell_use),2), std(press_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
    hold on;
    errorbar(ts{id}, mean(spont_norm(:,cell_use),2), std(spont_norm(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
    title([date_mat(id,:) ' ' mouse_mat(id,:) '-'  num2str(length(cell_use)) ' cells'])   
    xlabel('Time (ms)')
    ylabel('Norm dF/F')
end
suptitle('All cells- Event waveform normalized to spont')
print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');

figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    release_renorm_RS = bsxfun(@rdivide,release_norm_RS{id}, max(spont_norm_RS{id},[],1));
    press_renorm_RS = bsxfun(@rdivide,press_norm_RS{id}, max(spont_norm_RS{id},[],1));
    spont_renorm_RS = bsxfun(@rdivide,spont_norm_RS{id}, max(spont_norm_RS{id},[],1));
    release_norm_peak_RS{id} = max(release_renorm_RS,[],1);
    press_norm_peak_RS{id} = max(press_renorm_RS,[],1);
    spont_norm_peak_RS{id} = max(spont_renorm_RS,[],1);
    cell_use = find(~isnan(mean([release_renorm_RS; press_renorm_RS; spont_renorm_RS],1)));
    errorbar(ts{id}, mean(release_renorm_RS(:,cell_use),2), std(release_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'k')
    hold on;
    errorbar(ts{id}, mean(press_renorm_RS(:,cell_use),2), std(press_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'c')
    hold on;
    errorbar(ts{id}, mean(spont_renorm_RS(:,cell_use),2), std(spont_renorm_RS(:,cell_use),[],2)./sqrt(length(cell_use)), 'b')
    title([date_mat(id,:) ' ' mouse_mat(id,:) '-'  num2str(length(cell_use)) ' cells'])   
    xlabel('Time (ms)')
    ylabel('Norm dF/F')
end
suptitle('Responsive cells- Event waveform normalized to spont')
print([out_base 'Summary_avgevent_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgevent_norm2spont.pdf'], '-dpdf');

%peak scatters
r_all = [];
s_all = [];
p_all = [];
r_RS = [];
s_RS = [];
p_RS = [];
figure;
for id = 1:length(date_mat)
    r_all = [r_all release_norm_peak{id}];
    p_all = [p_all press_norm_peak{id}];
    r_RS = [r_RS release_norm_peak_RS{id}];
    p_RS = [p_RS press_norm_peak_RS{id}];   
    subplot(1,2,1)
    scatter(release_norm_peak{id}, press_norm_peak{id}, col_mat(id,:))
    hold on;
    subplot(1,2,2)
    scatter(release_norm_peak_RS{id}, press_norm_peak_RS{id}, col_mat(id,:))
    hold on;
end
subplot(1,2,1)
xlim([0 2.5])
ylim([0 2.5])
hline(1, '--k')
hold on;
vline(1, '--k')
use_cells_all = find(~isnan(mean([p_all; r_all],1)),2);
[h_rall p_rall] = ttest(r_all(1,use_cells_all),1);
[h_pall p_pall] = ttest(p_all(1,use_cells_all),1);
xlabel('Release amplitude')
ylabel('Press amplitude')
ncells = sum(~isnan(mean([p_all; r_all],1)),2);
title(['All cells- n = ' num2str(ncells) '; Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
subplot(1,2,2)
xlim([0 2.5])
ylim([0 2.5])
hline(1, '--k')
hold on;
vline(1, '--k')
use_cells_RS = find(~isnan(mean([p_all; r_all],1)),2);
[h_rRS p_rRS] = ttest(r_RS(1,use_cells_RS),1);
[h_pRS p_pRS] = ttest(p_RS(1,use_cells_RS),1);
xlabel('Release amplitude')
ylabel('Press amplitude')
ncells = sum(~isnan(mean([p_RS; r_RS],1)),2);
title(['Responsive cells- n = ' num2str(ncells) '; Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
suptitle('Event waveform normalized to spont')
print([out_base 'Summary_pr_scatter_norm2spont.eps'], '-depsc');
print([out_base 'Summary_pr_scatter_norm2spont.pdf'], '-dpdf');
    
%average across mice
figure;
for id = 1:length(date_mat) 
    subplot(1,2,1)
    temp_release = release_norm_peak{id};
    temp_press = press_norm_peak{id};
    use_cells_all = find(~isnan(mean([temp_release; temp_press],1)));
    temp_release_RS = release_norm_peak_RS{id};
    temp_press_RS = press_norm_peak_RS{id};
    use_cells_RS = find(~isnan(mean([temp_release_RS; temp_press_RS],1)));
    errorbarxy(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), std(temp_release(:,use_cells_all),[],2)./sqrt(length(use_cells_all)), std(temp_press(:,use_cells_all),[],2)./sqrt(length(use_cells_all)),[],[],[], col_mat(id,:))
    hold on;
    scatter(mean(temp_release(:,use_cells_all),2), mean(temp_press(:,use_cells_all),2), col_mat(id,:))
    subplot(1,2,2)
    errorbarxy(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2), std(temp_release_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)), std(temp_press_RS(:,use_cells_RS),[],2)./sqrt(length(use_cells_RS)),[],[],[], col_mat(id,:))
    hold on;
    scatter(mean(temp_release_RS(:,use_cells_RS),2), mean(temp_press_RS(:,use_cells_RS),2), col_mat(id,:))
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
xlim([0.5 1.5])
ylim([0.5 1.5])
hline(1, '--k')
hold on;
vline(1, '--k')
xlabel('Release amplitude')
ylabel('Press amplitude')
title('Responsive cells')
suptitle('Avg event waveform normalized to spont')
print([out_base 'Summary_avgpr_scatter_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avgpr_scatter_norm2spont.pdf'], '-dpdf');


%press and release alone
figure;
lever = strvcat('release', 'press');
for id = 1:length(date_mat) 
    subplot(2,2,1)
    scatter(ones(size(release_norm_peak{id})), release_norm_peak{id}, col_mat(id,:))
    hold on;
    scatter(ones(size(press_norm_peak{id})).*2, press_norm_peak{id}, col_mat(id,:))
    hold on;
    subplot(2,2,2)
    scatter(ones(size(release_norm_peak_RS{id})), release_norm_peak_RS{id}, col_mat(id,:))
    hold on;
    scatter(ones(size(press_norm_peak_RS{id})).*2, press_norm_peak_RS{id}, col_mat(id,:))
    hold on;
end
subplot(2,2,1)
xlim([0.5 2.5])
ylim([0 2.5])
hline(1, '--k')
set(gca, 'XTick', 1:2, 'XTickLabels', lever);
[h_rall p_rall] = ttest(r_all,1);
[h_pall p_pall] = ttest(p_all,1);
ylabel('Normalized amplitude')
title(['Release p = ' num2str(chop(p_rall,2)) '; Press p = ' num2str(chop(p_pall,2))])
subplot(2,2,2)
xlim([0.5 2.5])
ylim([0 2.5])
hline(1, '--k')
set(gca, 'XTick', 1:2,'XTickLabels', lever);
[h_rRS p_rRS] = ttest(r_RS,1);
[h_pRS p_pRS] = ttest(p_RS,1);
ylabel('Normalized amplitude')
title(['Release p = ' num2str(chop(p_rRS,2)) '; Press p = ' num2str(chop(p_pRS,2))])
    
%average across mice
for id = 1:length(date_mat) 
    subplot(2,2,3)
    temp_release = release_norm_peak{id};
    temp_press = press_norm_peak{id};
    temp_release_RS = release_norm_peak_RS{id};
    temp_press_RS = press_norm_peak_RS{id};
    errorbar(1, nanmean(temp_release,2),nanstd(temp_release,[],2)./sqrt(size(temp_release,2)), ['o' col_mat(id,:)])
    hold on;
    errorbar(2, nanmean(temp_press,2), nanstd(temp_press,[],2)./sqrt(size(temp_press,2)), ['o' col_mat(id,:)])
    hold on;
    subplot(2,2,4)
    errorbar(1, nanmean(temp_release_RS,2), nanstd(temp_release_RS,[],2)./sqrt(size(temp_release_RS,2)), ['o' col_mat(id,:)])
    hold on;
    errorbar(2, nanmean(temp_press_RS,2), nanstd(temp_press_RS,[],2)./sqrt(size(temp_press_RS,2)), ['o' col_mat(id,:)])
    hold on;
end
subplot(2,2,3)
xlim([0.5 2.5])
ylim([0.5 1.5])
hline(1, '--k')
set(gca, 'XTick', 1:2, 'XTickLabel', lever);
ylabel('Normalized amplitude')
title(['All cells- Release n = ' num2str(sum(~isnan(r_all),2)) ' cells; Press n = ' num2str(sum(~isnan(p_all),2))])
subplot(2,2,4)
xlim([0.5 2.5])
ylim([0.5 1.5])
hline(1, '--k')
set(gca, 'XTick', 1:2, 'XTickLabel', lever);
ylabel('Normalized amplitude')
title(['Responsive cells- Release n = ' num2str(sum(~isnan(r_RS),2)) ' cells; Press n = ' num2str(sum(~isnan(p_RS),2))])
suptitle('Avg event waveform normalized to spont')
print([out_base 'Summary_avg_event_norm2spont.eps'], '-depsc');
print([out_base 'Summary_avg_event_norm2spont.pdf'], '-dpdf');

%spont waveform comparison
all_press_norm = [];
RS_press_norm = [];
all_release_norm = [];
RS_release_norm = [];
all_spont_norm = [];
RS_spont_norm = [];
for id = 1:4
    release_renorm_all = release_norm_all{id}';
    press_renorm_all = press_norm_all{id}';
    spont_renorm_all = spont_norm_all{id}';
    release_renorm_RS = release_norm_RS{id}';
    press_renorm_RS = press_norm_RS{id}';
    spont_renorm_RS = spont_norm_RS{id}';
    all_press_norm = [all_press_norm; press_renorm_all];
    RS_press_norm = [RS_press_norm; press_renorm_RS];
    all_release_norm = [all_release_norm; release_renorm_all];
    RS_release_norm = [RS_release_norm; release_renorm_RS];
    all_spont_norm = [all_spont_norm; spont_renorm_all];
    RS_spont_norm = [RS_spont_norm; spont_renorm_RS];
end
for id = 6
    release_renorm_all = release_norm_all{id}';
    press_renorm_all = press_norm_all{id}';
    spont_renorm_all = spont_norm_all{id}';
    release_renorm_RS = release_norm_RS{id}';
    press_renorm_RS = press_norm_RS{id}';
    spont_renorm_RS = spont_norm_RS{id}';
    sz = size(press_renorm_all);
    all_press_norm = [all_press_norm(:,2:16); squeeze(mean(reshape(press_renorm_all(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    sz = size(press_renorm_RS);
    RS_press_norm = [RS_press_norm(:,2:16); squeeze(mean(reshape(press_renorm_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    sz = size(spont_renorm_all);
    all_spont_norm = [all_spont_norm(:,2:16); squeeze(mean(reshape(spont_renorm_all(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    sz = size(spont_renorm_RS);
    RS_spont_norm = [RS_spont_norm(:,2:16); squeeze(mean(reshape(spont_renorm_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    sz = size(release_renorm_all);
    all_release_norm = [all_release_norm(:,2:16); squeeze(mean(reshape(release_renorm_all(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
    sz = size(release_renorm_RS);
    RS_release_norm = [RS_release_norm(:,2:16); squeeze(mean(reshape(release_renorm_RS(:,1:sz(2)-1),[sz(1) 2 (sz(2)-1)./2]),2))];
end
figure;
all_spont_renorm = bsxfun(@rdivide, all_spont_norm, max(all_spont_norm,[],2));
all_release_renorm = bsxfun(@rdivide, all_release_norm, max(all_spont_norm,[],2));
all_press_renorm = bsxfun(@rdivide, all_press_norm, max(all_spont_norm,[],2));
RS_spont_renorm = bsxfun(@rdivide, RS_spont_norm, max(RS_spont_norm,[],2));
RS_release_renorm = bsxfun(@rdivide, RS_release_norm, max(RS_spont_norm,[],2));
RS_press_renorm = bsxfun(@rdivide, RS_press_norm, max(RS_spont_norm,[],2));
cell_use = find(~isnan(mean([all_spont_renorm all_release_renorm all_press_renorm],2)));
tts = ts{1};
subplot(1,2,1)
errorbar(tts(1,2:16), mean(all_spont_renorm(cell_use,:),1), std(all_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-b')
hold on;
errorbar(tts(1,2:16), mean(all_release_renorm(cell_use,:),1), std(all_release_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tts(1,2:16), mean(all_press_renorm(cell_use,:),1), std(all_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- all cells- ' num2str(length(cell_use)) ' cells'])

subplot(1,2,2)
cell_use = find(~isnan(mean([RS_spont_renorm RS_release_renorm RS_press_renorm],2)));
errorbar(tts(1,2:16), mean(RS_spont_renorm(cell_use,:),1), std(RS_spont_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-b')
hold on;
errorbar(tts(1,2:16), mean(RS_release_renorm(cell_use,:),1), std(RS_release_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-k')
hold on;
errorbar(tts(1,2:16), mean(RS_press_renorm(cell_use,:),1), std(RS_press_renorm(cell_use,:),[],1)./sqrt(length(cell_use)), '-c')
hold on;
xlabel('Time (ms)')
ylabel('dF/F')
title(['Average event- responsive cells- ' num2str(length(cell_use)) ' cells'])
suptitle('Average event - img24, img25 and img27 - Black: spont; Red: early; Cyan: press')
print([out_base 'Summary_event_TC_norm2spont.eps'], '-depsc');
print([out_base 'Summary_event_TC_norm2spont.pdf'], '-dpdf');
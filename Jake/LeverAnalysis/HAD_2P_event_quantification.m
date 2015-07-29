% quantifies amplitude and timing of events for each condition
% 1. calculate average timecourses for press/release events
% 2. compare timecourse of press/release events
% 3. compare PSTH for success/fail and press
% 4. compare amplitudes of evoked and spontaneous events


load([dest_sub '_spont_events.mat'])
load([dest_sub '_evoked_events.mat'])
load([dest '_parse_behavior.mat'])
load([dest_sub '_pvals.mat'])
nCells = size(events,2);
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n)
    n2= n-1;
else
    n2 = n;
end

%% 1- calculate average events for different conditions
data_start = round(300./double(ifi));
data_end = round(700./double(ifi));    

for ic = 1:nCells
    events(ic).dfoverf_avg = mean(events(ic).dfoverf_chunk,1);
    events(ic).dfoverf_sem = std(events(ic).dfoverf_chunk,[],1)./sqrt(size(events(ic).dfoverf_chunk,1));
    
    success(ic).good_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).good_event==1),:);
    success(ic).good_event_dfoverf_avg = mean(success(ic).good_event_dfoverf,1);
    success(ic).good_event_dfoverf_sem = std(success(ic).good_event_dfoverf,[],1)./sqrt(size(success(ic).good_event_dfoverf,1));
    success(ic).all_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).event==1),:);
    success(ic).all_event_dfoverf_avg = mean(success(ic).all_event_dfoverf,1);
    success(ic).all_event_dfoverf_sem = std(success(ic).all_event_dfoverf,[],1)./sqrt(size(success(ic).all_event_dfoverf,1));
    success(ic).no_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).event==0),:);
    success(ic).no_event_dfoverf_avg = mean(success(ic).no_event_dfoverf,1);
    success(ic).no_event_dfoverf_sem = std(success(ic).no_event_dfoverf,[],1)./sqrt(size(success(ic).no_event_dfoverf,1));
    success(ic).dfoverf_avg = mean(success(ic).dfoverf_chunk,1);
    success(ic).dfoverf_sem = std(success(ic).dfoverf_chunk,[],1)./sqrt(size(success(ic).dfoverf_chunk,1));

    fail(ic).good_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).good_event==1),:);
    fail(ic).good_event_dfoverf_avg = mean(fail(ic).good_event_dfoverf,1);
    fail(ic).good_event_dfoverf_sem = std(fail(ic).good_event_dfoverf,[],1)./sqrt(size(fail(ic).good_event_dfoverf,1));
    fail(ic).all_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).event==1),:);
    fail(ic).all_event_dfoverf_avg = mean(fail(ic).all_event_dfoverf,1);
    fail(ic).all_event_dfoverf_sem = std(fail(ic).all_event_dfoverf,[],1)./sqrt(size(fail(ic).all_event_dfoverf,1));
    fail(ic).no_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).event==0),:);
    fail(ic).no_event_dfoverf_avg = mean(fail(ic).no_event_dfoverf,1);
    fail(ic).no_event_dfoverf_sem = std(fail(ic).no_event_dfoverf,[],1)./sqrt(size(fail(ic).no_event_dfoverf,1));
    fail(ic).dfoverf_avg = mean(fail(ic).dfoverf_chunk,1);
    fail(ic).dfoverf_sem = std(fail(ic).dfoverf_chunk,[],1)./sqrt(size(fail(ic).dfoverf_chunk,1));

    press(ic).good_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).good_event==1),:);
    press(ic).good_event_dfoverf_avg = mean(press(ic).good_event_dfoverf,1);
    press(ic).good_event_dfoverf_sem = std(press(ic).good_event_dfoverf,[],1)./sqrt(size(press(ic).good_event_dfoverf,1));
    press(ic).all_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).event==1),:);
    press(ic).all_event_dfoverf_avg = mean(press(ic).all_event_dfoverf,1);
    press(ic).all_event_dfoverf_sem = std(press(ic).all_event_dfoverf,[],1)./sqrt(size(press(ic).all_event_dfoverf,1));
    press(ic).no_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).event==0),:);
    press(ic).no_event_dfoverf_avg = mean(press(ic).no_event_dfoverf,1);
    press(ic).no_event_dfoverf_sem = std(press(ic).no_event_dfoverf,[],1)./sqrt(size(press(ic).no_event_dfoverf,1));
    press(ic).dfoverf_avg = mean(press(ic).dfoverf_chunk,1);
    press(ic).dfoverf_sem = std(press(ic).dfoverf_chunk,[],1)./sqrt(size(press(ic).dfoverf_chunk,1));        
    
    release(ic).good_event_dfoverf_avg = mean([success(ic).good_event_dfoverf; fail(ic).good_event_dfoverf],1); 
    press(ic).good_event_dfoverf_peak = mean(press(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    success(ic).good_event_dfoverf_peak = mean(success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    fail(ic).good_event_dfoverf_peak = mean(fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    release(ic).good_event_dfoverf_peak = mean([success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))); fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)))],1);
end
save([dest_sub '_event_summary.mat'], 'data_start', 'data_end', 'press', 'release', 'success', 'fail', 'events');
%% plots
ts = [-data_start:data_end].*double(ifi);

%% compare timecourse of aligned good events
%by cells
figure; 
for ic = 1:nCells
    subplot(n,n2,ic)
    errorbar(ts, success(ic).good_event_dfoverf_avg, success(ic).good_event_dfoverf_sem, 'k')
    hold on
    errorbar(ts, fail(ic).good_event_dfoverf_avg, fail(ic).good_event_dfoverf_sem, 'r')
    hold on
    errorbar(ts, press(ic).good_event_dfoverf_avg, press(ic).good_event_dfoverf_sem, 'c')
    xlabel('Time (ms)')
    ylabel('dF/F')
    if release_h(ic)
        rel_col = '\color{red}';
    else
        rel_col = '\color{black}';
    end
    if press_h(ic)
        pr_col = '\color{red}';
    else
        pr_col = '\color{black}';
    end
    title([rel_col num2str(size(success(ic).good_event_dfoverf,1)) ' S; ' rel_col num2str(size(fail(ic).good_event_dfoverf,1)) ' F; ' pr_col num2str(size(press(ic).good_event_dfoverf,1)) ' P'])
end
suptitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
print([dest_sub '_good_event_dfoverf.eps'], '-depsc');
print([dest_sub '_good_event_dfoverf.pdf'], '-dpdf');

%summary of peaks
figure;
resp_h = zeros(size(press_h));
resp_h(find((press_h+release_h+success_h+fail_h))) = 1;
all_release_peak = zeros(1,nCells);
all_press_peak = zeros(1,nCells);
all_success_peak = zeros(1,nCells);
all_fail_peak = zeros(1,nCells);
all_release_peak = zeros(1,nCells);
lever = strvcat('press','success', 'fail');
subplot(2,2,1)
for ic = 1:nCells 
    plot(1:3, [press(ic).good_event_dfoverf_peak success(ic).good_event_dfoverf_peak fail(ic).good_event_dfoverf_peak], '-ok')
    hold on
    all_press_peak(1,ic) = press(ic).good_event_dfoverf_peak;
    all_success_peak(1,ic) = success(ic).good_event_dfoverf_peak;
    all_fail_peak(1,ic) = fail(ic).good_event_dfoverf_peak;
    all_release_peak(1,ic) = release(ic).good_event_dfoverf_peak;
    if resp_h(ic)
        plot(1:3, [press(ic).good_event_dfoverf_peak success(ic).good_event_dfoverf_peak fail(ic).good_event_dfoverf_peak], 'or')
        hold on
    end
end
set(gca, 'XTick', 1:3, 'XTickLabel', lever);
xlim([0.5 3.5])
title(['All cells- ' num2str(nCells) ' cells'])

subplot(2,2,2)
errorbar(1, nanmean(all_press_peak,2), nanstd(all_press_peak,[],2)./sqrt(sum(~isnan(all_press_peak(1,:)),2)), 'oc');
hold on;
errorbar(2, nanmean(all_success_peak,2), nanstd(all_success_peak,[],2)./sqrt(sum(~isnan(all_success_peak(1,:)),2)), 'ok');
hold on;
errorbar(3, nanmean(all_fail_peak,2), nanstd(all_fail_peak,[],2)./sqrt(sum(~isnan(all_fail_peak(1,:)),2)), 'or');
set(gca, 'XTick', 1:3, 'XTickLabel', lever);
xlim([0.5 3.5])
title(['All cells- ' num2str(nCells) ' cells'])

subplot(2,2,3)
resp_press_peak = zeros(1,sum(resp_h,2));
resp_success_peak = zeros(1,sum(resp_h,2));
resp_fail_peak = zeros(1,sum(resp_h,2));
start = 1;
for ic = find(resp_h)
    resp_press_peak(1,start) = press(ic).good_event_dfoverf_peak;
    resp_success_peak(1,start) = success(ic).good_event_dfoverf_peak;
    resp_fail_peak(1,start) = fail(ic).good_event_dfoverf_peak;
    start = start+1;
end
errorbar(1, nanmean(resp_press_peak,2), nanstd(resp_press_peak,[],2)./sqrt(sum(~isnan(resp_press_peak(1,:)),2)), 'oc');
hold on;
errorbar(2, nanmean(resp_success_peak,2), nanstd(resp_success_peak,[],2)./sqrt(sum(~isnan(resp_success_peak(1,:)),2)), 'ok');
hold on;
errorbar(3, nanmean(resp_fail_peak,2), nanstd(resp_fail_peak,[],2)./sqrt(sum(~isnan(resp_fail_peak(1,:)),2)), 'or');
set(gca, 'XTick', 1:3, 'XTickLabel', lever);
xlim([0.5 3.5])
title(['Responsive cells- ' num2str(sum(press_h,2)) ' cells'])
suptitle([mouse ' ' date '- Good events peak amp- Black: success; Red: early; Cyan: press'])
print([dest_sub '_good_event_peak_amp.eps'], '-depsc');
print([dest_sub '_good_event_peak_amp.pdf'], '-dpdf');

%average across cells
figure;
all_success_event = [];
all_fail_event = [];
all_press_event = [];
all_release_event = [];
for ic = 1:nCells
    all_success_event = [all_success_event;success(ic).good_event_dfoverf_avg];
    all_fail_event = [all_fail_event; fail(ic).good_event_dfoverf_avg];
    all_press_event = [all_press_event; press(ic).good_event_dfoverf_avg];
    all_release_event = [all_release_event; release(ic).good_event_dfoverf_avg];
end
subplot(1,2,1)
errorbar(ts, nanmean(all_success_event,1), nanstd(all_success_event,[],1)./sqrt(sum(~isnan(all_success_event(:,1)),1)), 'k');
hold on;
errorbar(ts, nanmean(all_fail_event,1), nanstd(all_fail_event,[],1)./sqrt(sum(~isnan(all_fail_event(:,1)),1)), 'r');
hold on;
errorbar(ts, nanmean(all_press_event,1), nanstd(all_press_event,[],1)./sqrt(sum(~isnan(all_press_event(:,1)),1)), 'c');
xlabel('Time (ms)')
ylabel('dF/F')
title(['All cells - Success (' num2str(sum(~isnan(all_success_event(:,1)),1)) '); Early (' num2str(sum(~isnan(all_fail_event(:,1)),1)) '); Press (' num2str(sum(~isnan(all_press_event(:,1)),1)) ')'])

subplot(1,2,2)
resp_success_event = [];
resp_fail_event = [];
resp_press_event = [];
for ic = find(resp_h)
    resp_success_event = [resp_success_event; success(ic).good_event_dfoverf_avg];
    resp_fail_event = [resp_fail_event; fail(ic).good_event_dfoverf_avg];
    resp_press_event = [resp_press_event; press(ic).good_event_dfoverf_avg];
end
errorbar(ts, nanmean(resp_success_event,1), nanstd(resp_success_event,[],1)./sqrt(sum(~isnan(resp_success_event(:,1)),1)), 'k');
hold on;
errorbar(ts, nanmean(resp_fail_event,1), nanstd(resp_fail_event,[],1)./sqrt(sum(~isnan(resp_fail_event(:,1)),1)), 'r');
hold on;
errorbar(ts, nanmean(resp_press_event,1), nanstd(resp_press_event,[],1)./sqrt(sum(~isnan(resp_press_event(:,1)),1)), 'c');
xlabel('Time (ms)')
ylabel('dF/F')
title(['Responsive cells - Success (' num2str(sum(~isnan(resp_success_event(:,1)),1)) '); Early (' num2str(sum(~isnan(resp_fail_event(:,1)),1)) '); Press (' num2str(sum(~isnan(resp_press_event(:,1)),1)) ')'])

suptitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
print([dest_sub '_good_event_dfoverf_avg.eps'], '-depsc');
print([dest_sub '_good_event_dfoverf_avg.pdf'], '-dpdf');

%compare with no events
all_success_noevent = [];
all_fail_noevent = [];
all_press_noevent = [];
resp_success_noevent = [];
resp_fail_noevent = [];
resp_press_noevent = [];

for ic = 1:nCells
    all_success_noevent = [all_success_noevent;success(ic).no_event_dfoverf_avg];
    all_fail_noevent = [all_fail_noevent; fail(ic).no_event_dfoverf_avg];
    all_press_noevent = [all_press_noevent; press(ic).no_event_dfoverf_avg];
end
for ic = find(resp_h)
    resp_success_noevent = [resp_success_noevent;success(ic).no_event_dfoverf_avg];
    resp_fail_noevent = [resp_fail_noevent; fail(ic).no_event_dfoverf_avg];
    resp_press_noevent = [resp_press_noevent; press(ic).no_event_dfoverf_avg];
end


subplot(3,3,1)
errorbar(ts, nanmean(all_success_event,1), nanstd(all_success_event,[],1)./sqrt(sum(~isnan(all_success_event(:,1)),1)), 'k');
hold on;
errorbar(ts, nanmean(all_success_noevent,1), nanstd(all_success_noevent,[],1)./sqrt(sum(~isnan(all_success_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['All cells: Event (' num2str(sum(~isnan(all_success_event(:,1)),1)) '); No event (' num2str(sum(~isnan(all_success_noevent(:,1)),1)) ')'])
subplot(3,3,2)
errorbar(ts, nanmean(all_fail_event,1), nanstd(all_fail_event,[],1)./sqrt(sum(~isnan(all_fail_event(:,1)),1)), 'r');
hold on;
errorbar(ts, nanmean(all_fail_noevent,1), nanstd(all_fail_noevent,[],1)./sqrt(sum(~isnan(all_fail_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['All cells: Event (' num2str(sum(~isnan(all_fail_event(:,1)),1)) '); No event (' num2str(sum(~isnan(all_fail_noevent(:,1)),1)) ')'])
subplot(3,3,3)
errorbar(ts, nanmean(all_press_event,1), nanstd(all_press_event,[],1)./sqrt(sum(~isnan(all_press_event(:,1)),1)), 'c');
hold on;
errorbar(ts, nanmean(all_press_noevent,1), nanstd(all_press_noevent,[],1)./sqrt(sum(~isnan(all_press_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['All cells: Event (' num2str(sum(~isnan(all_press_event(:,1)),1)) '); No event (' num2str(sum(~isnan(all_press_noevent(:,1)),1)) ')'])


subplot(3,3,4)
errorbar(ts, nanmean(resp_success_event,1), nanstd(resp_success_event,[],1)./sqrt(sum(~isnan(resp_success_event(:,1)),1)), 'k');
hold on;
errorbar(ts, nanmean(resp_success_noevent,1), nanstd(resp_success_noevent,[],1)./sqrt(sum(~isnan(resp_success_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['Responsive cells: Event (' num2str(sum(~isnan(resp_success_event(:,1)),1)) '); No event (' num2str(sum(~isnan(resp_success_noevent(:,1)),1)) ')'])
subplot(3,3,5)
errorbar(ts, nanmean(resp_fail_event,1), nanstd(resp_fail_event,[],1)./sqrt(sum(~isnan(resp_fail_event(:,1)),1)), 'r');
hold on;
errorbar(ts, nanmean(resp_fail_noevent,1), nanstd(resp_fail_noevent,[],1)./sqrt(sum(~isnan(resp_fail_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['Responsive cells: Event (' num2str(sum(~isnan(resp_fail_event(:,1)),1)) '); No event (' num2str(sum(~isnan(resp_fail_noevent(:,1)),1)) ')'])
subplot(3,3,6)
errorbar(ts, nanmean(resp_press_event,1), nanstd(resp_press_event,[],1)./sqrt(sum(~isnan(resp_press_event(:,1)),1)), 'c');
hold on;
errorbar(ts, nanmean(resp_press_noevent,1), nanstd(resp_press_noevent,[],1)./sqrt(sum(~isnan(resp_press_noevent(:,1)),1)), 'm');
xlabel('Time (ms)')
ylabel('dF/F')
title(['Responsive cells: Event (' num2str(sum(~isnan(resp_press_event(:,1)),1)) '); No event (' num2str(sum(~isnan(resp_press_noevent(:,1)),1)) ')'])

suptitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
print([dest_sub '_no_event_dfoverf_avg.eps'], '-depsc');
print([dest_sub '_no_event_dfoverf_avg.pdf'], '-dpdf');

save([dest_sub '_event_peaks'], 'all_success_event', 'all_fail_event', 'all_press_event', 'all_release_event', 'all_success_peak', 'all_fail_peak', 'all_press_peak', 'all_release_peak');

%% Create PSTHs
%figure;
% plot(repmat(spike_frames,1,2)', bsxfun(@times,[ie:ie+1],ones(length(spike_frames),2))', 'k')

sz1 = size(success_tc,1)-1;
sz_s = size(success(1).f_chunk,1);
sz_f = size(fail(1).f_chunk,1);
all_success_hist = zeros(sz1,nCells);
all_fail_hist = zeros(sz1,nCells);
for ic = 1:nCells
    success(ic).hist = zeros(sz1,sz_s);
    for is = 1:sz_s
        spike_frames = success(ic).ind{is};
        success(ic).hist(spike_frames,is) = 1;
    end
    fail(ic).hist = zeros(sz1,sz_f);
    for ie = 1:sz_f
        spike_frames = fail(ic).ind{ie};
        fail(ic).hist(spike_frames,ie) = 1;
    end
    all_success_hist(:,ic) = mean(success(ic).hist,2);
    all_fail_hist(:,ic) = mean(fail(ic).hist,2);
end

figure;
ts = [1-(sz1/2):(sz1/2)].*double(ifi);
h_ev_succ = zeros(1,nCells);
p_ev_succ = zeros(1,nCells);
h_ev_fail = zeros(1,nCells);
p_ev_fail = zeros(1,nCells);
for ic = 1:nCells
    subplot(n,n2,ic)
    errorbar(ts, mean(success(ic).hist.*(1000/double(ifi)),2), std(success(ic).hist.*(1000/double(ifi)),[],2)./sqrt(size(success(ic).hist,2)),'k')
    hold on
    errorbar(ts, mean(fail(ic).hist.*(1000/double(ifi)),2), std(fail(ic).hist.*(1000/double(ifi)),[],2)./sqrt(size(fail(ic).hist,2)),'r')
    xlim([ts(1) ts(end)])
    if success_h(ic)
        s_col = '\color{red}';
    else
        s_col = '\color{black}';
    end
    if fail_h(ic)
        f_col = '\color{red}';
    else
        f_col = '\color{black}';
    end
    [h_ev_succ(1,ic), p_ev_succ(1,ic)] = ttest(mean(success(ic).hist(1:ceil(500/double(ifi)),:),1), mean(success(ic).hist((sz1/2):(sz1/2)+ceil(65/double(ifi)),:),1),'tail', 'left');
    [h_ev_fail(1,ic), p_ev_fail(1,ic)] = ttest(mean(fail(ic).hist(1:ceil(500/double(ifi)),:),1), mean(fail(ic).hist((sz1/2):(sz1/2)+ceil(65/double(ifi)),:),1),'tail', 'left');
    if h_ev_succ(1,ic)
        s_col = ['\bf' s_col];
    end
    if h_ev_fail(1,ic)
        f_col = ['\bf' f_col];
    end
    title([s_col 'Success ' f_col ' Early'])
end
suptitle([date ' ' mouse ' PSTH all events- black: success (' num2str(size(success(ic).hist,2)) '); red: fail ('  num2str(size(fail(ic).hist,2)) ')'])
print([dest_sub '_all_event_PSTH.eps'], '-depsc');
print([dest_sub '_all_event_PSTH.pdf'], '-dpdf');

figure;
subplot(2,1,1)
errorbar(ts,mean(all_success_hist.*(1000/double(ifi)),2),std(all_success_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_success_hist,2)),'k')
hold on
errorbar(ts,mean(all_fail_hist.*(1000/double(ifi)),2),std(all_fail_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_fail_hist,2)),'r')
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['All cells- n = ' num2str(size(all_fail_hist,2))])

resp_success_hist = zeros(sz1,sum(resp_h,2));
resp_fail_hist = zeros(sz1,sum(resp_h,2));
start = 1;
for ic = find(resp_h)
    success_hist = zeros(sz1,sz_s);
    for is = 1:sz_s
        spike_frames = success(ic).ind{is};
        success_hist(spike_frames,is) = 1;
    end
    fail_hist = zeros(sz1,sz_f);
    for ie = 1:sz_f
        spike_frames = fail(ic).ind{ie};
        fail_hist(spike_frames,ie) = 1;
    end
    resp_success_hist(:,start) = mean(success_hist,2);
    resp_fail_hist(:,start) = mean(fail_hist,2);
    start= start+1;
end
subplot(2,1,2)
errorbar(ts,mean(resp_success_hist.*(1000/double(ifi)),2),std(resp_success_hist.*(1000/double(ifi)),[],2)./sqrt(size(resp_success_hist,2)),'k')
hold on
errorbar(ts,mean(resp_fail_hist.*(1000/double(ifi)),2),std(resp_fail_hist.*(1000/double(ifi)),[],2)./sqrt(size(resp_fail_hist,2)),'r')
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['Resp cells- n = ' num2str(size(resp_success_hist,2))])
suptitle([date ' ' mouse 'PSTH all events- black: success; red: fail'])

print([dest_sub '_all_event_avgPSTH.eps'], '-depsc');
print([dest_sub '_all_event_avgPSTH.pdf'], '-dpdf');

%find peak and latency of PSTH
figure
lever = strvcat('success', 'fail');
success_ind_all = [];
success_rate_all = [];
fail_ind_all = [];
fail_rate_all = [];
for ic = 1:nCells
    release(ic).hist = [success(ic).hist fail(ic).hist];
    bases = mean(mean(success(ic).hist(1:(sz1/2)-ceil(65/double(ifi)),:),1),2);
    [success_rate, success_ind_all(ic)] = max(mean(success(ic).hist((sz1/2):(sz1/2)+ceil(100/double(ifi)),:),2),[],1);
    basef = mean(mean(fail(ic).hist(1:sz1/2-ceil(65/double(ifi)),:),1),2);
    [fail_rate, fail_ind_all(ic)] = max(mean(fail(ic).hist((sz1/2):(sz1/2)+ceil(100/double(ifi)),:),2),[],1);
    subplot(4,2,1)
    success_rate_all(ic) = success_rate-bases;
    fail_rate_all(ic) = fail_rate-basef;
    plot(1:2, [success_rate-bases fail_rate-basef].*(1000/double(ifi)), '-ok')
    hold on
    if release_h(ic)
        plot(1:2, [success_rate-bases fail_rate-basef].*(1000/double(ifi)), 'or')
        hold on
    end
    subplot(4,2,2)
    plot(1:2, [success_ind_all(ic) fail_ind_all(ic)].*double(ifi), '-ok')
    hold on
end
subplot(4,2,1)
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Rate (Hz)')
title('Peak rate')
subplot(4,2,2)
for ic = 1:nCells
    if release_h(ic)
        plot(1:2, [success_ind_all(ic) fail_ind_all(ic)].*double(ifi), 'or')
        hold on
    end
end
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Time (ms)')
title('Latency of peak rate increase')
subplot(4,2,3)
errorbar(1:2, [mean(success_rate_all,2) mean(fail_rate_all,2)].*(1000/double(ifi)), [std(success_rate_all,[],2)./sqrt(nCells) std(fail_rate_all,[],2)./sqrt(nCells)].*(1000/double(ifi)), 'ok')
[h_amp, p_amp] = ttest(success_rate_all,fail_rate_all);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Rate (Hz)')
title(['Peak rate- ' num2str(nCells) ' cells- p = ' num2str(chop(p_amp,2))])
subplot(4,2,4)
errorbar(1:2, [mean(success_ind_all,2) mean(fail_ind_all,2)].*double(ifi), [std(success_ind_all,[],2)./sqrt(nCells) std(fail_ind_all,[],2)./sqrt(nCells)].*double(ifi), 'ok')
[h_lat, p_lat] = ttest(success_ind_all,fail_ind_all);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Time (ms)')
title(['Latency- ' num2str(nCells) ' cells- p = ' num2str(chop(p_lat,2))])


success_ind_RS = [];
success_rate_RS = [];
fail_ind_RS = [];
fail_rate_RS = [];
start = 1;
for ic = find(resp_h)
    bases = mean(mean(success(ic).hist(1:(sz1/2)-ceil(65/double(ifi)),:),1),2);
    [success_rate, success_ind_RS(start)] = max(mean(success(ic).hist((sz1/2):(sz1/2)+ceil(100/double(ifi)),:),2),[],1);
    basef = mean(mean(fail(ic).hist(1:sz1/2-ceil(65/double(ifi)),:),1),2);
    [fail_rate, fail_ind_RS(start)] = max(mean(fail(ic).hist((sz1/2):(sz1/2)+ceil(100/double(ifi)),:),2),[],1);
    success_rate_RS(start) = success_rate-bases;
    fail_rate_RS(start) = fail_rate-basef;
    start = start+1;
end

subplot(4,2,5)
errorbar(1:2, [mean(success_rate_RS,2) mean(fail_rate_RS,2)].*(1000/double(ifi)), [std(success_rate_RS,[],2)./sqrt(sum(resp_h,2)) std(fail_rate_RS,[],2)./sqrt(sum(resp_h,2))].*(1000/double(ifi)), 'ok')
[h_amp, p_amp] = ttest(success_rate_RS,fail_rate_RS);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Rate (Hz)')
title(['Release resp- Peak rate- ' num2str(sum(release_h,2)) ' cells- p = ' num2str(chop(p_amp,2))])
subplot(4,2,6)
errorbar(1:2, [mean(success_ind_all,2) mean(fail_ind_all,2)].*double(ifi), [std(success_ind_all,[],2)./sqrt(sum(resp_h,2)) std(fail_ind_all,[],2)./sqrt(sum(resp_h,2))].*double(ifi), 'ok')
[h_lat, p_lat] = ttest(success_ind_all,fail_ind_all);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Time (ms)')
title(['Responsive- Latency- ' num2str(sum(resp_h,2)) ' cells- p = ' num2str(chop(p_lat,2))])
suptitle([date ' ' mouse ' release evoked event frequency and latency'])

print([dest_sub '_PSTH_amp_latency.eps'], '-depsc');
print([dest_sub '_PSTH_amp_latency.pdf'], '-dpdf');

save([dest_sub '_event_hist.mat'], 'success_ind_RS', 'success_rate_RS', 'fail_ind_RS', 'fail_rate_RS', 'success_ind_all', 'success_rate_all', 'fail_ind_all', 'fail_rate_all', 'all_success_hist', 'all_fail_hist', 'resp_success_hist', 'resp_fail_hist');
% %plot at 3x lower samp rate
% figure;
% dx = squeeze(mean(reshape(x(1:30,:),[3,10,nCells]),1));
% dy = squeeze(mean(reshape(y(1:30,:),[3,10,nCells]),1));
% dts = squeeze(mean(reshape(ts(1,1:30),[3,10]),1));
% figure;
% errorbar(dts,mean(dx.*(1000/double(ifi)),2),std(dx.*(1000/double(ifi)),[],2)./sqrt(size(dx,2)),'k')
% hold on
% errorbar(dts,mean(dy.*(1000/double(ifi)),2),std(dy.*(1000/double(ifi)),[],2)./sqrt(size(dx,2)),'r')
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% suptitle('PSTH all events, all cells- black: success; red: fail')
% print([dest_sub '_all_event_avgPSTH_3bin.eps'], '-depsc');
% print([dest_sub '_all_event_avgPSTH_3bin.pdf'], '-dpdf');

%% PSTH of press events
% z = [];
% sz1 = size(success_tc,1)-1;
% sz_p = size(press(1).f_chunk,1);
% for ic = 1:nCells
%     press(ic).hist = zeros(sz1,sz_p);
%     for ip = 1:sz_p
%         spike_frames = press(ic).ind{ip};
%         press(ic).hist(spike_frames,ip) = 1;
%     end
%     if size(press(ic).hist,2)>1
%         z = [z mean(press(ic).hist,2)];
%     end
% end
% 
% figure;
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     errorbar(ts, mean(press(ic).hist.*(1000/double(ifi)),2), std(press(ic).hist.*(1000/double(ifi)),[],2)./sqrt(size(press(ic).hist,2)),'c')
%     xlim([ts(1) ts(end)])
%     if press_h(ic)
%         title('Resp')
%     end
% end
% suptitle('PSTH all events- Press')
% print([dest_sub '_all_event_press_PSTH.eps'], '-depsc');
% print([dest_sub '_all_event_press_PSTH.pdf'], '-dpdf');
% 
% figure;
% subplot(2,1,1)
% errorbar(ts,mean(z.*(1000/double(ifi)),2),std(z.*(1000/double(ifi)),[],2)./sqrt(size(z,2)),'c')
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% title(['All cells- n = ' num2str(size(z,2))])
% 
% z = [];
% sz1 = size(success_tc,1)-1;
% sz_p = size(press(1).f_chunk,1);
% for ic = find(press_h)
%     press(ic).hist = zeros(sz1,sz_p);
%     for ip = 1:sz_p
%         spike_frames = press(ic).ind{ip};
%         press(ic).hist(spike_frames,ip) = 1;
%     end
%     if size(press(ic).hist,2)>1
%         z = [z mean(press(ic).hist,2)];
%     end
% end
% subplot(2,1,2)
% errorbar(ts,mean(z.*(1000/double(ifi)),2),std(z.*(1000/double(ifi)),[],2)./sqrt(size(z,2)),'c')
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% title(['Press responsive cells- n = ' num2str(size(z,2))])
% suptitle('PSTH all events- Press')
% print([dest_sub '_press_avgPSTH.eps'], '-depsc');
% print([dest_sub '_press_avgPSTH.pdf'], '-dpdf');

% %plot at 3x lower samp rate
% figure;
% dz = squeeze(mean(reshape(z(1:30,:),[3,10,nCells]),1));
% dts = squeeze(mean(reshape(ts(1,1:30),[3,10]),1));
% figure;
% errorbar(dts,mean(dz.*(1000/double(ifi)),2),std(dz.*(1000/double(ifi)),[],2)./sqrt(size(dz,2)),'c')
% xlabel('Time (ms)')
% ylabel('Firing rate (Hz)')
% suptitle('PSTH all events, all cells- Press- binned')
% print([dest_sub '_all_event_press_avgPSTH_3bin.eps'], '-depsc');
% print([dest_sub '_all_event_press_avgPSTH_3bin.pdf'], '-dpdf');

%% compare spontaneous and evoked event amplitudes and waveforms
%compare amplitude distributions of triggered events
% figure
% a = [];
% b = [];
% c = [];
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     [Hsp] = cdfplot_LG(events(ic).dfoverf_peak);
%     hold on
%     set(Hsp, 'Color', 'b')
%     if sum(success(ic).event,2)
%         success_peak = success(ic).dfoverf_peak(find(success(ic).event==1),:);
%         [Hsu] = cdfplot_LG(success_peak);
%         hold on
%         set(Hsu, 'Color', 'k')
%     end
%     if sum(fail(ic).event,2)
%         fail_peak = fail(ic).dfoverf_peak(find(fail(ic).event==1),:);
%         [Hfa] = cdfplot_LG(fail_peak);
%         hold on
%         set(Hfa, 'Color', 'r')
%     end
%     title([num2str(size(events(ic).dfoverf_peak,1)) ' spont; ' num2str(sum(success(ic).event,2)) ' success;' num2str(sum(fail(ic).event,2)) ' fail; ' num2str(sum(press(ic).event,2)) ' press'])
%     a = [a; sum(success(ic).event,2)];
%     b = [b; sum(fail(ic).event,2)];
%     c = [c; size(events(ic).dfoverf_peak,1)];
% end
% suptitle('Amplitude distribution- all events')
% print([dest_sub '_all_event_spontVrelease_ampdist.eps'], '-depsc');
% print([dest_sub '_all_event_spontVrelease_ampdist.pdf'], '-dpdf');
% 
% d = [];
% figure;
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     d = [d; sum(press(ic).event,2)];
%     [Hsp] = cdfplot_LG(events(ic).dfoverf_peak);
%     hold on
%     set(Hsp, 'Color', 'b')
%     press_peak = press(ic).dfoverf_peak(find(press(ic).event==1),:);
%     [Hpr] = cdfplot_LG(press_peak);
%     hold on
%     set(Hpr, 'Color', 'c')
%     title([num2str(size(events(ic).dfoverf_peak,1)) ' spont; ' num2str(sum(press(ic).event,2)) ' press'])
% end
% suptitle('Amplitude distribution- all events')
% print([dest_sub 'all_event_spontVpress_ampdist.eps'], '-depsc');
% print([dest_sub 'all_event_spontVpress_ampdist.pdf'], '-dpdf');
% 
% 
% for ic = 1:nCells
%     max_event = max([success(ic).dfoverf_peak(find(success(ic).event==1)); fail(ic).dfoverf_peak(find(fail(ic).event==1)); events(ic).dfoverf_peak; press(ic).dfoverf_peak(find(press(ic).event==1))],[],1);
%     success(ic).peak_norm = success(ic).dfoverf_peak(find(success(ic).event==1))./max_event;
%     fail(ic).peak_norm = fail(ic).dfoverf_peak(find(fail(ic).event==1))./max_event;
%     event(ic).peak_norm = events(ic).dfoverf_peak./max_event;
%     press(ic).peak_norm = press(ic).dfoverf_peak(find(press(ic).event==1))./max_event;
% end
% 
% %chose n_min samples from success (fail is always 36, spont is 107) for
% %each ROI- overlay CDFs and perform KS test- title has median ks test for
% %each pairing.
% for n_min = [5 10]
%     figure;
%     nrep = 100;
%     kp_succVfail = zeros(1,nrep);
%     kp_succVspont = zeros(1,nrep);
%     kp_failVspont = zeros(1,nrep);
%     tp_succVfail = zeros(1,nrep);
%     tp_succVspont = zeros(1,nrep);
%     tp_failVspont = zeros(1,nrep);
%     for irep  = 1:nrep
%         success_peak_norm = [];
%         fail_peak_norm = [];
%         spont_peak_norm = [];
%         ncells = 0;
%         for ic = 1:nCells
%             if size(success(ic).peak_norm,1)>n_min-1
%                 samp1 = randsample(sum(success(ic).event,2),n_min);
%                 samp2 = randsample(sum(fail(ic).event,2),min(b,[],1));
%                 samp3 = randsample(size(events(ic).dfoverf_peak,1),min(c,[],1));
%                 success_peak_norm = [success_peak_norm; success(ic).peak_norm(samp1,:)];
%                 fail_peak_norm = [fail_peak_norm; fail(ic).peak_norm(samp2,:)];
%                 spont_peak_norm = [spont_peak_norm; event(ic).peak_norm(samp3,:)];
%                 ncells = ncells+1;
%             end
%         end
%         if and(size(success_peak_norm,1)>0,size(fail_peak_norm,1)>0)
%             [h_succ] = cdfplot(success_peak_norm);
%             set(h_succ,'Color','k')
%             hold on;
%             [h_fail] = cdfplot(fail_peak_norm);
%             set(h_fail,'Color','r')
%             hold on
%             [h_spont] = cdfplot(spont_peak_norm);
%             set(h_spont,'Color','b')
%             [kh_succVfail kp_succVfail(:,irep)] = kstest2(success_peak_norm,fail_peak_norm);
%             [kh_succVspont kp_succVspont(:,irep)] = kstest2(success_peak_norm,spont_peak_norm);
%             [kh_failVspont kp_failVspont(:,irep)] = kstest2(spont_peak_norm,fail_peak_norm);
%             [th_succVfail tp_succVfail(:,irep)] = ttest2(success_peak_norm,fail_peak_norm);
%             [th_succVspont tp_succVspont(:,irep)] = ttest2(success_peak_norm,spont_peak_norm);
%             [th_failVspont tp_failVspont(:,irep)] = ttest2(spont_peak_norm,fail_peak_norm);
%         end
%     end
%     xk = sort(kp_succVfail,2);
%     yk = sort(kp_failVspont,2);
%     zk = sort(kp_succVspont,2);
%     xt = sort(tp_succVfail,2);
%     yt = sort(tp_failVspont,2);
%     zt = sort(tp_succVspont,2);
%     title([num2str(n_min) ' events; ' num2str(ncells) ' cells; succVfail = ' num2str(xk(50)) '; succVspont = ' num2str(zk(50)) '; failVspont = ' num2str(yk(50))])
%     print([dest_sub '_all_event_ampdist_spontVrelease_resamp_' num2str(n_min) 'events.eps'], '-depsc');
%     print([dest_sub '_all_event_ampdist_spontVrelease_resamp_' num2str(n_min) 'events.pdf'], '-dpdf');
% end
% 
% 
% %press vs spont
% press_peak_norm = [];
% spont_peak_norm = [];
% for ic = 1:nCells
%     samp1 = randsample(sum(press(ic).event,2),min(d,[],1));
%     samp2 = randsample(size(events(ic).dfoverf_peak,1),min(c,[],1));
%     press_peak_norm = [press_peak_norm; press(ic).peak_norm(samp1,:)];
%     spont_peak_norm = [spont_peak_norm; event(ic).peak_norm(samp2,:)];
% end
% figure;
% [h_pr] = cdfplot(press_peak_norm);
% set(h_pr,'Color','c')
% hold on;
% [h_spont] = cdfplot(spont_peak_norm);
% set(h_spont,'Color','b')
% [h_prVsp p_prVsp] = kstest2(press_peak_norm, spont_peak_norm);
% title(['Event amplitude distribution: Press V. Spont - p = ' num2str(p_prVsp)])
% print([dest_sub '_all_event_ampdist_pressVspont.eps'], '-depsc');
% print([dest_sub '_all_event_ampdist_pressVspont.pdf'], '-dpdf');
% 
% %combine success and fails
% figure
% x = [];
% y = [];
% ab = a+b;
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     [Hsp] = cdfplot_LG(events(ic).dfoverf_peak);
%     hold on
%     set(Hsp, 'Color', 'b')
%     peaks = [success(ic).dfoverf_peak(find(success(ic).event==1)); fail(ic).dfoverf_peak(find(fail(ic).event==1))];
%     [Hev] = cdfplot_LG(peaks);
%     hold on
%     set(Hev, 'Color', 'k')
%     [h p] = kstest2(events(ic).dfoverf_peak,peaks);
%     title([num2str(p) ' ' num2str(size(peaks,1))])
%     samp1 = randsample(size(events(ic).dfoverf_peak,1), min(c,[],1));
%     samp2 = randsample(size(peaks,1),min(ab,[],1));
%     x = [x; events(ic).dfoverf_peak(samp1,:)./max(events(ic).dfoverf_peak,[],1)];
%     y = [y; peaks(samp2,:)./max(peaks,[],1)];
% end
% suptitle('Amplitude distribution- all events')
% print([dest_sub '_all_event_ampdist_spontVrelease.eps'], '-depsc');
% print([dest_sub '_all_event_ampdist_spontVrelease.pdf'], '-dpdf');
% 
% figure
% [Hsp_mean] = cdfplot_LG(x);
% hold on
% set(Hsp_mean, 'Color', 'b')
% [Hev_mean] = cdfplot_LG(y);
% hold on
% set(Hev_mean, 'Color', 'k')
% xlim([0 1])
% [h p] = kstest2(x,y);
% title(['Amplitude distribution- all events- Spont vs Evoked- p = ' num2str(p)])
% print([dest_sub '_all_event_ampdist_spontVrelease_avg.eps'], '-depsc');
% print([dest_sub '_all_event_ampdist_spontVrelease_avg.pdf'], '-dpdf');
% 
% %compare evoked, press and spont for good events only
% a = [];
% b = [];
% c = [];
% for ic = 1:nCells
%     a = [a; nansum([success(ic).good_event fail(ic).good_event],2)];
%     b = [b; size(events(ic).dfoverf_peak,1)];
%     c = [c; nansum(press(ic).good_event,2)];
% end
% 
% x = [];
% y = [];
% a = sort(a);
% b = sort(b);
% for s = 1:3:10
%     ncells = 0;
%     figure
%     for ic = 1:nCells
%         subplot(n,n2,ic)
%         [Hsp] = cdfplot_LG(events(ic).dfoverf_peak);
%         hold on
%         set(Hsp, 'Color', 'b')
%         peaks = [success(ic).dfoverf_peak(find(success(ic).good_event==1)); fail(ic).dfoverf_peak(find(fail(ic).good_event==1))];
%         [Hev] = cdfplot_LG(peaks);
%         hold on
%         set(Hev, 'Color', 'k')
%         [h p] = kstest2(events(ic).dfoverf_peak,peaks);
%         title([num2str(p) ' ' num2str(size(peaks,1))])
%         if size(peaks,1)>a(s)-1
%             samp1 = randsample(size(events(ic).dfoverf_peak,1),b(1));
%             samp2 = randsample(size(peaks,1),a(s));
%             x = [x; events(ic).dfoverf_peak(samp1,:)./max(events(ic).dfoverf_peak,[],1)];
%             y = [y; peaks(samp2,:)./max(peaks,[],1)];
%             ncells = ncells+1;
%         end
%     end
%     suptitle(['Amplitude distribution- good events- Spont (n = ' num2str(b(1)) ') vs Evoked (n = ' num2str(a(s)) ')'])
% 
%     print([dest_sub '_good_event_ampdist_releaseVspont_' num2str(a(s)) 'events.eps'], '-depsc');
%     print([dest_sub '_good_event_ampdist_releaseVspont_' num2str(a(s)) 'events.pdf'], '-dpdf');
% 
%     figure
%     [Hsp_mean] = cdfplot_LG(x);
%     hold on
%     set(Hsp_mean, 'Color', 'b')
%     [Hev_mean] = cdfplot_LG(y);
%     hold on
%     set(Hev_mean, 'Color', 'k')
%     xlim([0 1])
%     [h p] = kstest2(x,y);
%     title(['Amplitude distribution- good events- ' num2str(ncells) 'cells; Spont vs Evoked - P = ' num2str(p)])
%     print([dest_sub '_good_event_allcells_ampdist_releaseVspont_' num2str(a(s)) 'events.eps'], '-depsc');
%     print([dest_sub '_good_event_allcells_ampdist_releaseVspont_' num2str(a(s)) 'events.pdf'], '-dpdf');
% end
% 
% %compare spontaneous and press distributions
% x = [];
% y = [];
% b = sort(b);
% figure
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     [Hsp] = cdfplot_LG(events(ic).dfoverf_peak);
%     hold on
%     set(Hsp, 'Color', 'b')
%     peaks = press(ic).dfoverf_peak(find(press(ic).good_event==1));
%     [Hpr] = cdfplot_LG(peaks);
%     hold on
%     set(Hpr, 'Color', 'c')
%     if c(ic)>6
%         [h p] = kstest2(events(ic).dfoverf_peak,peaks);
%         title([num2str(p) ' ' num2str(size(peaks,1))])
%         samp1 = randsample(size(events(ic).dfoverf_peak,1),b(1));
%         samp2 = randsample(size(peaks,1),7);
%         x = [x; events(ic).dfoverf_peak(samp1,:)./max(events(ic).dfoverf_peak,[],1)];
%         y = [y; peaks(samp2,:)./max(peaks,[],1)];
%     end
% end
% suptitle(['Amplitude distribution- good events- Spont (n = ' num2str(b(1)) ') vs Press (n = ' num2str(c(1)) ')'])
% print([dest_sub '_good_event_ampdist_pressVspont.eps'], '-depsc');
% print([dest_sub '_good_event_ampdist_pressVspont.pdf'], '-dpdf');
% 
% figure
% [Hsp_mean] = cdfplot_LG(x);
% hold on
% set(Hsp_mean, 'Color', 'b')
% [Hpr_mean] = cdfplot_LG(y);
% hold on
% set(Hpr_mean, 'Color', 'c')
% xlim([0 1])
% [h p] = kstest2(x,y);
% title(['Amplitude distribution- good events; Spont vs Press - P = ' num2str(p)])
% print([dest_sub '_good_event_allcells_ampdist_pressVspont.eps'], '-depsc');
% print([dest_sub '_good_event_allcells_ampdist_pressVspont.pdf'], '-dpdf');
% 
% %compare evoked and spontaneous waveforms
% figure;
% spont_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
% evoked_norm_avg = zeros(size(success(1).dfoverf_chunk,2),nCells);
% spont_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
% evoked_norm_sem = zeros(size(success(1).dfoverf_chunk,2),nCells);
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     max_amp = mean(events(ic).dfoverf_peak,1); %normalizing to the average spont, not max
%     spont_norm = events(ic).dfoverf_chunk./max_amp;
%     evoked_norm = [success(ic).dfoverf_chunk(find(success(ic).good_event==1),:)./max_amp;  fail(ic).dfoverf_chunk(find(fail(ic).good_event==1),:)./max_amp];
%     spont_norm_avg(:,ic) = squeeze(mean(spont_norm,1));
%     evoked_norm_avg(:,ic) = squeeze(mean(evoked_norm,1));
%     spont_norm_sem(:,ic) = squeeze(std(spont_norm,[],1))./sqrt(size(spont_norm,1));
%     evoked_norm_sem(:,ic) = squeeze(std(evoked_norm,[],1))./sqrt(size(evoked_norm,1));
%     errorbar(spont_norm_avg(:,ic), spont_norm_sem(:,ic), 'b')
%     hold on
%     errorbar(evoked_norm_avg(:,ic), evoked_norm_sem(:,ic), 'k')
%     title([num2str(size(spont_norm,1)) ' spont; ' num2str(size(evoked_norm,1)) ' evoked'])
% end
% suptitle(['Average good events- Spont and Evoked- norm to average Spont'])
% print([dest_sub '_good_event_releaseVspont.eps'], '-depsc');
% print([dest_sub '_good_event_releaseVspont.pdf'], '-dpdf');
% 
% evoked_renorm_avg = bsxfun(@rdivide, evoked_norm_avg, max(spont_norm_avg,[],1));
% spont_renorm_avg = bsxfun(@rdivide, spont_norm_avg, max(spont_norm_avg,[],1));
% figure;
% errorbar(mean(evoked_renorm_avg,2), std(evoked_renorm_avg,[],2)./sqrt(nCells),'k')
% hold on
% errorbar(mean(spont_renorm_avg,2), std(spont_renorm_avg,[],2)./sqrt(nCells),'b')
% title(['Average good events- Spont and Evoked- norm to average Spont'])
% print([dest_sub '_good_event_avg_releaseVspont.eps'], '-depsc');
% print([dest_sub '_good_event_avg_releaseVspont.pdf'], '-dpdf');
% 
% %split evoked back into success and fail
% figure;
% spont_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
% success_norm_avg = zeros(size(success(1).dfoverf_chunk,2),nCells);
% fail_norm_avg = zeros(size(fail(1).dfoverf_chunk,2),nCells);
% spont_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
% success_norm_sem = zeros(size(success(1).dfoverf_chunk,2),nCells);
% fail_norm_sem = zeros(size(fail(1).dfoverf_chunk,2),nCells);
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     max_amp = mean(events(ic).dfoverf_peak,1); %normalizing to the average spont, not max
%     spont_norm = events(ic).dfoverf_chunk./max_amp;
%     success_norm = success(ic).dfoverf_chunk(find(success(ic).good_event==1),:)./max_amp;
%     fail_norm = fail(ic).dfoverf_chunk(find(fail(ic).good_event==1),:)./max_amp;
%     spont_norm_avg(:,ic) = squeeze(mean(spont_norm,1));
%     success_norm_avg(:,ic) = squeeze(mean(success_norm,1));
%     fail_norm_avg(:,ic) = squeeze(mean(fail_norm,1));
%     spont_norm_sem(:,ic) = squeeze(std(spont_norm,[],1))./sqrt(size(spont_norm,1));
%     success_norm_sem(:,ic) = squeeze(std(success_norm,[],1))./sqrt(size(success_norm,1));
%     fail_norm_sem(:,ic) = squeeze(std(fail_norm,[],1))./sqrt(size(fail_norm,1));
%     errorbar(spont_norm_avg(:,ic), spont_norm_sem(:,ic), 'b')
%     hold on
%     errorbar(success_norm_avg(:,ic), success_norm_sem(:,ic), 'k')
%     hold on
%     errorbar(fail_norm_avg(:,ic), fail_norm_sem(:,ic), 'r')
%     title([num2str(size(spont_norm,1)) ' spont; ' num2str(size(success_norm,1)) ' success; ' num2str(size(fail_norm,1)) ' fail'])
% end
% suptitle(['Average good events- Spont and Evoked- norm to average Spont'])
% print([dest_sub '_good_event_successVfailVspont.eps'], '-depsc');
% print([dest_sub '_good_event_successVfailVspont.pdf'], '-dpdf');
% 
% success_renorm_avg = bsxfun(@rdivide, success_norm_avg, max(spont_norm_avg,[],1));
% fail_renorm_avg = bsxfun(@rdivide, fail_norm_avg, max(spont_norm_avg,[],1));
% spont_renorm_avg = bsxfun(@rdivide, spont_norm_avg, max(spont_norm_avg,[],1));
% figure;
% errorbar(nanmean(success_renorm_avg,2), nanstd(success_renorm_avg,[],2)./sqrt(nCells),'k')
% hold on
% errorbar(nanmean(fail_renorm_avg,2), nanstd(fail_renorm_avg,[],2)./sqrt(nCells),'r')
% hold on
% errorbar(nanmean(spont_renorm_avg,2), nanstd(spont_renorm_avg,[],2)./sqrt(nCells),'b')
% title(['Average good events- Spont(' num2str(sum(~isnan(spont_norm_avg(1,:)),2)) ')/Success(' num2str(sum(~isnan(success_norm_avg(1,:)),2)) ')/Fail(' num2str(sum(~isnan(fail_norm_avg(1,:)),2)) ')- norm to average Spont'])
% print([dest_sub '_good_event_avg_successVfailVspont.eps'], '-depsc');
% print([dest_sub '_good_event_avg_successVfailVspont.pdf'], '-dpdf');
% 
% %compare press and spontaneous waveforms
% figure;
% spont_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
% press_norm_avg = zeros(size(press(1).dfoverf_chunk,2),nCells);
% spont_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
% press_norm_sem = zeros(size(press(1).dfoverf_chunk,2),nCells);
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     max_amp = mean(events(ic).dfoverf_peak,1); %normalizing to the average spont, not max
%     spont_norm = events(ic).dfoverf_chunk./max_amp;
%     press_norm = press(ic).dfoverf_chunk(find(press(ic).good_event==1),:)./max_amp;
%     spont_norm_avg(:,ic) = squeeze(mean(spont_norm,1));
%     press_norm_avg(:,ic) = squeeze(mean(press_norm,1));
%     spont_norm_sem(:,ic) = squeeze(std(spont_norm,[],1))./sqrt(size(spont_norm,1));
%     press_norm_sem(:,ic) = squeeze(std(press_norm,[],1))./sqrt(size(press_norm,1));
%     errorbar(spont_norm_avg(:,ic), spont_norm_sem(:,ic), 'b')
%     hold on
%     errorbar(press_norm_avg(:,ic), press_norm_sem(:,ic), 'c')
%     title([num2str(size(spont_norm,1)) ' spont; ' num2str(size(press_norm,1)) ' press'])
% end
% suptitle(['Average good events- Spont and Press- norm to average Spont'])
% print([dest_sub '_good_event_pressVspont.eps'], '-depsc');
% print([dest_sub '_good_event_pressVspont.pdf'], '-dpdf');
% 
% press_renorm_avg = bsxfun(@rdivide, press_norm_avg, max(spont_norm_avg,[],1));
% spont_renorm_avg = bsxfun(@rdivide, spont_norm_avg, max(spont_norm_avg,[],1));
% figure;
% errorbar(nanmean(press_renorm_avg,2), nanstd(press_renorm_avg,[],2)./sqrt(nCells),'c')
% hold on
% errorbar(nanmean(spont_renorm_avg,2), nanstd(spont_renorm_avg,[],2)./sqrt(nCells),'b')
% title(['Average good events- Spont(' num2str(sum(~isnan(spont_norm_avg(1,:)),2)) ') and Press(' num2str(sum(~isnan(press_norm_avg(1,:)),2)) ')- norm to average Spont'])
% print([dest_sub '_good_event_avg_pressVspont.eps'], '-depsc');
% print([dest_sub '_good_event_avg_pressVspont.pdf'], '-dpdf');

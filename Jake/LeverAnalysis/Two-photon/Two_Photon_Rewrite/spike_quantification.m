% quantifies amplitude and timing of events for each condition
% 1. calculate average timecourses for press/release events
% 2. compare timecourse of press/release events
% 3. compare PSTH for success/fail and press
% 4. compare amplitudes of evoked and spontaneous events


load([dest_sub '_spont_events.mat'])
load([dest_sub '_evoked_events.mat'])
load([dest 'parse_behavior.mat'])
load([dest_sub '_pvals.mat'])
load([dest_sub '_cell_categories.mat']);

nCells = size(events,2);
RS_cells = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);

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
    events(ic).dfoverf_avg = nanmean(events(ic).dfoverf_chunk,1);
    events(ic).dfoverf_sem = std(events(ic).dfoverf_chunk,[],1)./sqrt(size(events(ic).dfoverf_chunk,1));
    
    success(ic).good_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).good_event==1),:);
    if isempty( success(ic).dfoverf_chunk(find(success(ic).good_event==1),:) )
        success(ic).good_event_dfoverf_avg = zeros(size(nanmean(success(ic).good_event_dfoverf,1)));
    else
        success(ic).good_event_dfoverf_avg = nanmean(success(ic).good_event_dfoverf,1);
    end
    success(ic).good_event_dfoverf_sem = std(success(ic).good_event_dfoverf,[],1)./sqrt(size(success(ic).good_event_dfoverf,1));
    success(ic).all_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).event==1),:);
    success(ic).all_event_dfoverf_avg = nanmean(success(ic).all_event_dfoverf,1);
    success(ic).all_event_dfoverf_sem = std(success(ic).all_event_dfoverf,[],1)./sqrt(size(success(ic).all_event_dfoverf,1));
    success(ic).no_event_dfoverf = success(ic).dfoverf_chunk(find(success(ic).event==0),:);
    success(ic).no_event_dfoverf_avg = nanmean(success(ic).no_event_dfoverf,1);
    success(ic).no_event_dfoverf_sem = std(success(ic).no_event_dfoverf,[],1)./sqrt(size(success(ic).no_event_dfoverf,1));
    success(ic).dfoverf_avg = nanmean(success(ic).dfoverf_chunk,1);
    success(ic).dfoverf_sem = std(success(ic).dfoverf_chunk,[],1)./sqrt(size(success(ic).dfoverf_chunk,1));

    fail(ic).good_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).good_event==1),:);
    if isempty( fail(ic).dfoverf_chunk(find(fail(ic).good_event==1),:) )
        fail(ic).good_event_dfoverf_avg = zeros(size(nanmean(fail(ic).good_event_dfoverf,1)));
    else
        fail(ic).good_event_dfoverf_avg = nanmean(fail(ic).good_event_dfoverf,1);
    end
    fail(ic).good_event_dfoverf_sem = std(fail(ic).good_event_dfoverf,[],1)./sqrt(size(fail(ic).good_event_dfoverf,1));
    fail(ic).all_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).event==1),:);
    fail(ic).all_event_dfoverf_avg = nanmean(fail(ic).all_event_dfoverf,1);
    fail(ic).all_event_dfoverf_sem = std(fail(ic).all_event_dfoverf,[],1)./sqrt(size(fail(ic).all_event_dfoverf,1));
    fail(ic).no_event_dfoverf = fail(ic).dfoverf_chunk(find(fail(ic).event==0),:);
    fail(ic).no_event_dfoverf_avg = nanmean(fail(ic).no_event_dfoverf,1);
    fail(ic).no_event_dfoverf_sem = std(fail(ic).no_event_dfoverf,[],1)./sqrt(size(fail(ic).no_event_dfoverf,1));
    fail(ic).dfoverf_avg = nanmean(fail(ic).dfoverf_chunk,1);
    fail(ic).dfoverf_sem = std(fail(ic).dfoverf_chunk,[],1)./sqrt(size(fail(ic).dfoverf_chunk,1));

    press(ic).good_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).good_event==1),:);
    if isempty( press(ic).dfoverf_chunk(find(press(ic).good_event==1),:) )
        press(ic).good_event_dfoverf_avg = zeros(size(nanmean(press(ic).good_event_dfoverf,1)));
    else
        press(ic).good_event_dfoverf_avg = nanmean(press(ic).good_event_dfoverf,1);
    end
    press(ic).good_event_dfoverf_sem = std(press(ic).good_event_dfoverf,[],1)./sqrt(size(press(ic).good_event_dfoverf,1));
    press(ic).all_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).event==1),:);
    press(ic).all_event_dfoverf_avg = nanmean(press(ic).all_event_dfoverf,1);
    press(ic).all_event_dfoverf_sem = std(press(ic).all_event_dfoverf,[],1)./sqrt(size(press(ic).all_event_dfoverf,1));
    press(ic).no_event_dfoverf = press(ic).dfoverf_chunk(find(press(ic).event==0),:);
    press(ic).no_event_dfoverf_avg = nanmean(press(ic).no_event_dfoverf,1);
    press(ic).no_event_dfoverf_sem = std(press(ic).no_event_dfoverf,[],1)./sqrt(size(press(ic).no_event_dfoverf,1));
    press(ic).dfoverf_avg = nanmean(press(ic).dfoverf_chunk,1);
    press(ic).dfoverf_sem = std(press(ic).dfoverf_chunk,[],1)./sqrt(size(press(ic).dfoverf_chunk,1));        
    
    release(ic).good_event_dfoverf = [success(ic).good_event_dfoverf; fail(ic).good_event_dfoverf]; 
    
    release(ic).good_event_dfoverf_avg = nanmean(release(ic).good_event_dfoverf,1);
    
    press(ic).good_event_dfoverf_peak = nanmean(press(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    success(ic).good_event_dfoverf_peak = nanmean(success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    fail(ic).good_event_dfoverf_peak = nanmean(fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))),1);
    release(ic).good_event_dfoverf_peak = nanmean([success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))); fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)))],1);
    events(ic).good_event_dfoverf_peak = nanmean(events(ic).dfoverf_chunk(:,data_start+ceil(100/double(ifi))),1);

    press(ic).good_event_dfoverf_peaks = press(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)));
    success(ic).good_event_dfoverf_peaks = success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)));
    fail(ic).good_event_dfoverf_peaks = fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)));
    release(ic).good_event_dfoverf_peaks = [success(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi))); fail(ic).good_event_dfoverf(:,data_start+ceil(100/double(ifi)))];
    events(ic).good_event_dfoverf_peaks = events(ic).dfoverf_chunk(:,data_start+ceil(100/double(ifi)));    
end
save([dest_sub '_event_summary.mat'], 'data_start', 'data_end', 'press', 'release', 'success', 'fail', 'events');
%% plots
ts = [-data_start:data_end].*double(ifi);

%% compare timecourse of aligned good events- evoked
%overlay individual events
fig=figure
for ic = 1:nCells
    subplot(n,n2,ic)
    ind = find(success(ic).good_event==1);
    if length(ind)>0
        plot(ts, success(ic).dfoverf_chunk(ind,:), 'k')
        hold on
    end
    ind2 = find(fail(ic).good_event==1);
    if length(ind2)>0
        plot(ts, fail(ic).dfoverf_chunk(ind2,:), 'r')
    end
end
supertitle([mouse ' ' date '- Good events- Black: success; Red: early'])
print([dest_sub '_good_event_overlay.eps'], '-depsc');
print([dest_sub '_good_event_overlay.pdf'], '-dpdf');

%by cells
fig=figure; 
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
supertitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
saveas(fig, [dest_sub '_good_event_dfoverf.fig']);
print([dest_sub '_good_event_dfoverf.eps'], '-depsc');
print([dest_sub '_good_event_dfoverf.pdf'], '-dpdf');

%summary of peaks
fig=figure;
resp_h = zeros(size(press_h));
resp_h(find((press_h+release_h+success_h+fail_h+tooFast_h))) = 1;
resp_h(find((press_h+release_h+success_h+fail_h))) = 1;
resp_h(tooFast_resp_cells) = 1;
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
supertitle([mouse ' ' date '- Good events peak amp- Black: success; Red: early; Cyan: press'])
saveas(fig, [dest_sub '_good_event_peak_amp.fig']);
print([dest_sub '_good_event_peak_amp.eps'], '-depsc');
print([dest_sub '_good_event_peak_amp.pdf'], '-dpdf');

%average across cells
fig=figure;
all_success_event = [];
all_fail_event = [];
all_press_event = [];
all_release_event = [];
all_spont_event = [];
for ic = 1:nCells
    all_success_event = [all_success_event;success(ic).good_event_dfoverf_avg];
    all_fail_event = [all_fail_event; fail(ic).good_event_dfoverf_avg];
    all_press_event = [all_press_event; press(ic).good_event_dfoverf_avg];
    all_release_event = [all_release_event; release(ic).good_event_dfoverf_avg];
    all_spont_event = [all_spont_event; events(ic).dfoverf_avg];
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

supertitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
saveas(fig, [dest_sub '_good_event_dfoverf_avg.fig']);
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

fig=figure;
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

supertitle([mouse ' ' date '- Good events- Black: success; Red: early; Cyan: press'])
saveas(fig, [dest_sub '_no_event_dfoverf_avg.fig']);
print([dest_sub '_no_event_dfoverf_avg.eps'], '-depsc');
print([dest_sub '_no_event_dfoverf_avg.pdf'], '-dpdf');

save([dest_sub '_event_peaks'], 'all_success_event', 'all_fail_event', 'all_press_event', 'all_release_event', 'all_success_peak', 'all_fail_peak', 'all_press_peak', 'all_release_peak', 'all_spont_event');

%% Create PSTHs
% fig=figure;
% plot(repmat(spike_frames,1,2)', bsxfun(@times,[ie:ie+1],ones(length(spike_frames),2))', 'k')

sz1 = size(success_tc,1);
sz1_half = floor(sz1/2);
sz_s = size(success(1).f_chunk,1);
sz_f = size(fail(1).f_chunk,1);
all_success_hist = zeros(sz1,nCells);
all_fail_hist = zeros(sz1,nCells);
all_success_syn = zeros(sz1,sz_s,nCells);
all_fail_syn = zeros(sz1,sz_f,nCells);


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
    all_success_syn(:,:,ic) = success(ic).hist;
    all_fail_syn(:,:,ic) = fail(ic).hist;
end

all_success_syn = permute(all_success_syn,[1 3 2]); % frames cells events
success_syn_all = zeros(1,sz_s);
for ie = 1:size(all_success_syn,3)
    spike_per_event_f = all_success_syn(sz1_half-round(300/ifi):sz1_half,:,ie);
    spike_per_event_l = all_success_syn(sz1_half+2:sz1_half+2+round(300/ifi),:,ie);
    zero_event = all_success_syn(sz1_half+1,:,ie);
    
    success_syn_all(ie) = calcEventStd(spike_per_event_f, spike_per_event_l, zero_event)*ifi;
end

RS_success_syn = all_success_syn(:,RS_cells,:);
success_syn_RS = zeros(1,sz_s);
for ie = 1:size(all_success_syn,3)
    spike_per_event_f = RS_success_syn(sz1_half-round(300/ifi):sz1_half,:,ie);
    spike_per_event_l = RS_success_syn(sz1_half+2:sz1_half+2+round(300/ifi),:,ie);
    zero_event = RS_success_syn(sz1_half+1,:,ie);
    
    success_syn_RS(ie) = calcEventStd(spike_per_event_f, spike_per_event_l, zero_event)*ifi;
end

success_syn_cell = zeros(1, nCells);
for ic = 1:size(all_success_syn,2)
    spike_per_cell_f = squeeze(all_success_syn(sz1_half-round(300/ifi):sz1_half,ic,:));
    spike_per_cell_l = squeeze(all_success_syn(sz1_half+2:sz1_half+round(300/ifi),ic,:));
    zero_event = squeeze(all_success_syn(sz1_half+1,ic,:));
    success_syn_cell(ic) = calcEventStd(spike_per_cell_f, spike_per_cell_l, zero_event)*ifi;
    
end
success_syn_cell_RS = success_syn_cell(RS_cells);

all_fail_syn = permute(all_fail_syn,[1 3 2]); % frames cells events
fail_syn_all = zeros(1,sz_f);
for ie = 1:size(all_fail_syn,3)
    spike_per_event_f = all_fail_syn(sz1_half-round(300/ifi):sz1_half,:,ie);
    spike_per_event_l = all_fail_syn(sz1_half+2:sz1_half+2+round(300/ifi),:,ie);
    zero_event = all_fail_syn(sz1_half+1,:,ie);
    fail_syn_all(ie) =  calcEventStd(spike_per_cell_f, spike_per_cell_l, zero_event)*ifi;
end

RS_fail_syn = all_fail_syn(:,RS_cells,:);
fail_syn_RS = zeros(1, sz_f);
for ie = 1:size(all_fail_syn,3)
    spike_per_event_f = RS_fail_syn(sz1_half-round(300/ifi):sz1_half,:,ie);
    spike_per_event_l = RS_fail_syn(sz1_half+2:sz1_half+2+round(300/ifi),:,ie);
    zero_event = RS_fail_syn(sz1_half+1,:,ie);
    fail_syn_RS(ie) =  calcEventStd(spike_per_cell_f, spike_per_cell_l, zero_event)*ifi;
end

fail_syn_cell = zeros(1, nCells);
for ic = 1:size(all_fail_syn,2)
    spike_per_cell_f = squeeze(all_fail_syn(sz1_half-round(300/ifi):sz1_half,ic,:));
    spike_per_cell_l = squeeze(all_fail_syn(sz1_half+2:sz1_half+2+round(300/ifi),ic,:));
    zero_event = squeeze(all_fail_syn(sz1_half+1,ic,:));
    fail_syn_cell(ic) = calcEventStd(spike_per_cell_f, spike_per_cell_l, zero_event)*ifi;
end

fail_syn_cell_RS = fail_syn_cell(RS_cells);

fig=figure;
ts = [-(sz1_half):(sz1_half)].*double(ifi);
h_ev_succ = zeros(1,nCells);
p_ev_succ = zeros(1,nCells);
h_ev_fail = zeros(1,nCells);
p_ev_fail = zeros(1,nCells);
for ic = 1:nCells
    subplot(n,n2,ic)
    errorbar(ts, nanmean(success(ic).hist.*(1000/double(ifi)),2), std(success(ic).hist.*(1000/double(ifi)),[],2)./sqrt(size(success(ic).hist,2)),'k')
    hold on
    errorbar(ts, nanmean(fail(ic).hist.*(1000/double(ifi)),2), std(fail(ic).hist.*(1000/double(ifi)),[],2)./sqrt(size(fail(ic).hist,2)),'r')
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
    [h_ev_succ(1,ic), p_ev_succ(1,ic)] = ttest(mean(success(ic).hist(1:ceil(500/double(ifi)),:),1), nanmean(success(ic).hist((sz1_half):(sz1_half)+ceil(65/double(ifi)),:),1),'tail', 'left');
    [h_ev_fail(1,ic), p_ev_fail(1,ic)] = ttest(mean(fail(ic).hist(1:ceil(500/double(ifi)),:),1), nanmean(fail(ic).hist((sz1_half):(sz1_half)+ceil(65/double(ifi)),:),1),'tail', 'left');
    if ~isnan(h_ev_succ(1,ic))
        if h_ev_succ(1,ic)
            s_col = ['\bf' s_col];   
        end
    else
        s_col = 'NaN';
    end
    if ~isnan(h_ev_fail(1,ic))
        if h_ev_fail(1,ic)
            f_col = ['\bf' f_col];
        end
    else
        s_col = 'NaN';
    end
    title([s_col 'Success ' f_col ' Early'])
end
supertitle([date ' ' mouse ' PSTH all events- black: success (' num2str(size(success(ic).hist,2)) '); red: fail ('  num2str(size(fail(ic).hist,2)) ')'])
saveas(fig, [dest_sub '_all_event_PSTH.fig']);
print([dest_sub '_all_event_PSTH.eps'], '-depsc');
print([dest_sub '_all_event_PSTH.pdf'], '-dpdf');

fig=figure;
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
    resp_success_hist(:,start) = nanmean(success_hist,2);
    resp_fail_hist(:,start) = nanmean(fail_hist,2);
    start= start+1;
end
subplot(2,1,2)
errorbar(ts,mean(resp_success_hist.*(1000/double(ifi)),2),std(resp_success_hist.*(1000/double(ifi)),[],2)./sqrt(size(resp_success_hist,2)),'k')
hold on
errorbar(ts,mean(resp_fail_hist.*(1000/double(ifi)),2),std(resp_fail_hist.*(1000/double(ifi)),[],2)./sqrt(size(resp_fail_hist,2)),'r')
xlabel('Time (ms)')
ylabel('Firing rate (Hz)')
title(['Resp cells- n = ' num2str(size(resp_success_hist,2))])
supertitle([date ' ' mouse 'PSTH all events- black: success; red: fail'])
saveas(fig, [dest_sub '_all_event_avgPSTH.fig']);
print([dest_sub '_all_event_avgPSTH.eps'], '-depsc');
print([dest_sub '_all_event_avgPSTH.pdf'], '-dpdf');

%find peak and latency of PSTH
fig=figure;
lever = strvcat('success', 'fail');
success_ind_all = [];
success_rate_all = [];
fail_ind_all = [];
fail_rate_all = [];
for ic = 1:nCells
    release(ic).hist = [success(ic).hist fail(ic).hist];
    bases = nanmean(mean(success(ic).hist(1:(sz1_half)-ceil(65/double(ifi)),:),1),2);
    [success_rate, success_ind_all(ic)] = max(mean(success(ic).hist((sz1_half):(sz1_half)+ceil(100/double(ifi)),:),2),[],1);
    basef = nanmean(mean(fail(ic).hist(1:sz1_half-ceil(65/double(ifi)),:),1),2);
    [fail_rate, fail_ind_all(ic)] = max(mean(fail(ic).hist((sz1_half):(sz1_half)+ceil(100/double(ifi)),:),2),[],1);
    subplot(4,2,1)
%     success_rate_all(ic) = success_rate-bases;
%     fail_rate_all(ic) = fail_rate-basef;
    success_rate_all(ic) = success_rate;
    fail_rate_all(ic) = fail_rate;
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
errorbar(1:2, [mean(success_rate_all,2) nanmean(fail_rate_all,2)].*(1000/double(ifi)), [std(success_rate_all,[],2)./sqrt(nCells) std(fail_rate_all,[],2)./sqrt(nCells)].*(1000/double(ifi)), 'ok')
[h_amp, p_amp] = ttest(success_rate_all,fail_rate_all);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Rate (Hz)')
title(['Peak rate- ' num2str(nCells) ' cells- p = ' num2str(chop(p_amp,2))])
subplot(4,2,4)
errorbar(1:2, [mean(success_ind_all,2) nanmean(fail_ind_all,2)].*double(ifi), [std(success_ind_all,[],2)./sqrt(nCells) std(fail_ind_all,[],2)./sqrt(nCells)].*double(ifi), 'ok')
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
    bases = nanmean(mean(success(ic).hist(1:(sz1_half)-ceil(65/double(ifi)),:),1),2);
    [success_rate, success_ind_RS(start)] = max(mean(success(ic).hist((sz1_half):(sz1_half)+ceil(100/double(ifi)),:),2),[],1);
    basef = nanmean(mean(fail(ic).hist(1:sz1_half-ceil(65/double(ifi)),:),1),2);
    [fail_rate, fail_ind_RS(start)] = max(mean(fail(ic).hist((sz1_half):(sz1_half)+ceil(100/double(ifi)),:),2),[],1);
%     success_rate_RS(start) = success_rate-bases;
%     fail_rate_RS(start) = fail_rate-basef;
    success_rate_RS(start) = success_rate;
    fail_rate_RS(start) = fail_rate;
    start = start+1;
end

subplot(4,2,5)
errorbar(1:2, [mean(success_rate_RS,2) nanmean(fail_rate_RS,2)].*(1000/double(ifi)), [std(success_rate_RS,[],2)./sqrt(sum(resp_h,2)) std(fail_rate_RS,[],2)./sqrt(sum(resp_h,2))].*(1000/double(ifi)), 'ok')
[h_amp, p_amp] = ttest(success_rate_RS,fail_rate_RS);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Rate (Hz)')
title(['Release resp- Peak rate- ' num2str(sum(release_h,2)) ' cells- p = ' num2str(chop(p_amp,2))])
subplot(4,2,6)
errorbar(1:2, [mean(success_ind_all,2) nanmean(fail_ind_all,2)].*double(ifi), [std(success_ind_all,[],2)./sqrt(sum(resp_h,2)) std(fail_ind_all,[],2)./sqrt(sum(resp_h,2))].*double(ifi), 'ok')
[h_lat, p_lat] = ttest(success_ind_all,fail_ind_all);
set(gca, 'XTick', 1:2, 'XTickLabel', lever)
xlim([0.5 2.5])
ylabel('Time (ms)')
title(['Responsive- Latency- ' num2str(sum(resp_h,2)) ' cells- p = ' num2str(chop(p_lat,2))])
supertitle([date ' ' mouse ' release evoked event frequency and latency'])
saveas(fig, [dest_sub '_PSTH_amp_latency.fig']);
print([dest_sub '_PSTH_amp_latency.eps'], '-depsc');
print([dest_sub '_PSTH_amp_latency.pdf'], '-dpdf');

save([dest_sub '_event_hist.mat'], 'success_ind_RS', 'success_rate_RS', 'fail_ind_RS', 'fail_rate_RS', 'success_ind_all', ...
    'success_rate_all', 'fail_ind_all', 'fail_rate_all', 'all_success_hist', 'all_fail_hist', 'resp_success_hist', 'resp_fail_hist',...
    'fail_syn_all', 'fail_syn_cell', 'success_syn_all', 'success_syn_cell', 'fail_syn_RS', 'fail_syn_cell_RS', 'success_syn_cell', ...
    'success_syn_cell_RS', 'success_syn_RS');

%% compare spontaneous and evoked event amplitudes and waveforms
%compare amplitude distributions of triggered events
fig=figure;
spontN = [];
releaseN = [];
pressN = [];
for ic = 1:nCells
    subplot(n,n2,ic)
    [Hsp] = cdfplot_LG(events(ic).good_event_dfoverf_peaks);
    hold on
    set(Hsp, 'Color', 'b')
    if size(release(ic).good_event_dfoverf_peaks,1)
        [Hsu] = cdfplot(release(ic).good_event_dfoverf_peaks);
        hold on
        set(Hsu, 'Color', 'k')
    end
    if size(press(ic).good_event_dfoverf_peaks,1)
        [Hsu] = cdfplot(press(ic).good_event_dfoverf_peaks);
        hold on
        set(Hsu, 'Color', 'c')
    end
    if resp_h(ic)
        sig_str = '\bf';
    else
        sig_str = '';
    end
    title([sig_str num2str(size(events(ic).good_event_dfoverf_peaks,1)) ' spont; ' num2str(size(release(ic).good_event_dfoverf_peaks,1)) ' release;' num2str(size(press(ic).good_event_dfoverf_peaks,1)) ' press'])
    spontN = [spontN; size(events(ic).good_event_dfoverf_peaks,1)];
    releaseN = [releaseN; size(release(ic).good_event_dfoverf_peaks,1)];
    pressN = [pressN; size(press(ic).good_event_dfoverf_peaks,1)];
end
supertitle('Amplitude distribution- all events')
saveas(fig, [dest_sub '_good_event_spontVevoked_ampdist.fig']);
print([dest_sub '_good_event_spontVevoked_ampdist.eps'], '-depsc');
print([dest_sub '_good_event_spontVevoked_ampdist.pdf'], '-dpdf');

for ic = 1:nCells
    max_event = max([events(ic).good_event_dfoverf_peaks; release(ic).good_event_dfoverf_peaks; press(ic).good_event_dfoverf_peaks],[],1);
    release(ic).peak_norm = release(ic).good_event_dfoverf_peaks./max_event;
    events(ic).peak_norm = events(ic).good_event_dfoverf_peaks./max_event;
    press(ic).peak_norm = press(ic).good_event_dfoverf_peaks./max_event;
end

%collect events across cells- not keeping track of event number
release_peak_norm = [];
press_peak_norm = [];
spont_peak_norm = [];
for ic = 1:nCells 
    release_peak_norm = [release_peak_norm; release(ic).peak_norm];
    press_peak_norm = [press_peak_norm; press(ic).peak_norm];
    spont_peak_norm = [spont_peak_norm; events(ic).peak_norm];
end
fig=figure;
subplot(1,2,1)
[h_rel] = cdfplot(release_peak_norm);
set(h_rel,'Color','k')
hold on;
[h_press] = cdfplot(press_peak_norm);
set(h_press,'Color','c')
hold on
[h_spont] = cdfplot(spont_peak_norm);
set(h_spont,'Color','b')
[kh_relVpr kp_relVpr] = kstest2(release_peak_norm,press_peak_norm);
[kh_relVspont kp_relVspont] = kstest2(release_peak_norm,spont_peak_norm);
[kh_prVspont kp_prVspont] = kstest2(spont_peak_norm,press_peak_norm);
title(['All cells- n = ' num2str(nCells) ' cells']) % rVp = ' num2str(chop(kp_relVpr,2)) '; rVs = ' num2str(chop(kp_relVspont,2)) '; pVs = ' num2str(chop(kp_prVspont,2))])

release_peak_norm = [];
press_peak_norm = [];
spont_peak_norm = [];
for ic = find(resp_h)
    release_peak_norm = [release_peak_norm; release(ic).peak_norm];
    press_peak_norm = [press_peak_norm; press(ic).peak_norm];
    spont_peak_norm = [spont_peak_norm; events(ic).peak_norm];
end
subplot(1,2,2)
[h_rel] = cdfplot(release_peak_norm);
set(h_rel,'Color','k')
hold on;
[h_press] = cdfplot(press_peak_norm);
set(h_press,'Color','c')
hold on
[h_spont] = cdfplot(spont_peak_norm);
set(h_spont,'Color','b')
[kh_relVpr kp_relVpr] = kstest2(release_peak_norm,press_peak_norm);
[kh_relVspont kp_relVspont] = kstest2(release_peak_norm,spont_peak_norm);
[kh_prVspont kp_prVspont] = kstest2(spont_peak_norm,press_peak_norm);
title(['Resp cells- n = ' num2str(sum(resp_h,2)) ' cells']) %rVp = ' num2str(chop(kp_relVpr,2)) '; rVs = ' num2str(chop(kp_relVspont,2)) '; pVs = ' num2str(chop(kp_prVspont,2))])
supertitle([date ' ' mouse 'peak event distribution- all good events'])
saveas(fig, [dest_sub '_ampdist_spontVevoked.fig']);
print([dest_sub '_ampdist_spontVevoked.eps'], '-depsc');
print([dest_sub '_ampdist_spontVevoked.pdf'], '-dpdf');

%collect events across cells- 7 press/release events per cell
release_peak_norm = [];
press_peak_norm = [];
spont_peak_norm = [];
cell_list = [];
for ic = 1:nCells 
    if pressN(ic)>6
        if releaseN(ic)>6
            samp_r = randsample(size(release(ic).peak_norm,1),7);
            samp_p = randsample(size(press(ic).peak_norm,1),7);
            release_peak_norm = [release_peak_norm; release(ic).peak_norm(samp_r,:)];
            press_peak_norm = [press_peak_norm; press(ic).peak_norm(samp_p,:)];
            spont_peak_norm = [spont_peak_norm; events(ic).peak_norm];
            cell_list = [cell_list ic];
        end
    end
end
fig=figure;
subplot(1,2,1)
[h_rel] = cdfplot(release_peak_norm);
set(h_rel,'Color','k')
hold on;
[h_press] = cdfplot(press_peak_norm);
set(h_press,'Color','c')
hold on
[h_spont] = cdfplot(spont_peak_norm);
set(h_spont,'Color','b')
[kh_relVpr kp_relVpr] = kstest2(release_peak_norm,press_peak_norm);
[kh_relVspont kp_relVspont] = kstest2(release_peak_norm,spont_peak_norm);
[kh_prVspont kp_prVspont] = kstest2(spont_peak_norm,press_peak_norm);
title(['All cells- n = ' num2str(length(cell_list)) ' cells' ])%rVp = ' num2str(chop(kp_relVpr,2)) '; rVs = ' num2str(chop(kp_relVspont,2)) '; pVs = ' num2str(chop(kp_prVspont,2))])

release_peak_norm = [];
press_peak_norm = [];
spont_peak_norm = [];
cell_list = [];
for ic = find(resp_h) 
    if pressN(ic)>6
        if releaseN(ic)>6
            samp_r = randsample(size(release(ic).peak_norm,1),7);
            samp_p = randsample(size(press(ic).peak_norm,1),7);
            release_peak_norm = [release_peak_norm; release(ic).peak_norm(samp_r,:)];
            press_peak_norm = [press_peak_norm; press(ic).peak_norm(samp_p,:)];
            spont_peak_norm = [spont_peak_norm; events(ic).peak_norm];
            cell_list = [cell_list ic];
        end
    end
end
subplot(1,2,2)
[h_rel] = cdfplot(release_peak_norm);
set(h_rel,'Color','k')
hold on;
[h_press] = cdfplot(press_peak_norm);
set(h_press,'Color','c')
hold on
[h_spont] = cdfplot(spont_peak_norm);
set(h_spont,'Color','b')
[kh_relVpr kp_relVpr] = kstest2(release_peak_norm,press_peak_norm);
[kh_relVspont kp_relVspont] = kstest2(release_peak_norm,spont_peak_norm);
[kh_prVspont kp_prVspont] = kstest2(spont_peak_norm,press_peak_norm);
title(['Resp cells- n = ' num2str(length(cell_list)) ' cells' ])%rVp = ' num2str(chop(kp_relVpr,2)) '; rVs = ' num2str(chop(kp_relVspont,2)) '; pVs = ' num2str(chop(kp_prVspont,2))])
supertitle([date ' ' mouse 'peak event distribution- at least 7 events per cell'])
saveas(fig, [dest_sub '_ampdist_spontVevoked_min7events.fig']);
print([dest_sub '_ampdist_spontVevoked_min7events.eps'], '-depsc');
print([dest_sub '_ampdist_spontVevoked_min7events.pdf'], '-dpdf');

%compare evoked and spontaneous waveforms
fig=figure;
spont_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
release_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
press_norm_avg = zeros(size(press(1).dfoverf_chunk,2),nCells);
spont_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
press_norm_sem = zeros(size(press(1).dfoverf_chunk,2),nCells);
release_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
success_norm_avg = zeros(size(success(1).dfoverf_chunk,2),nCells);
fail_norm_avg = zeros(size(fail(1).dfoverf_chunk,2),nCells);

for ic = 1:nCells
    subplot(n,n2,ic)
    max_amp = nanmean(events(ic).good_event_dfoverf_peaks,1); %normalizing to the average spont, not max
    max_amp = 1;
    spont_norm{ic} = events(ic).dfoverf_chunk./max_amp;
    release_norm{ic} = release(ic).good_event_dfoverf./max_amp;
    press_norm{ic} = press(ic).good_event_dfoverf./max_amp;
    fail_norm{ic} = fail(ic).good_event_dfoverf./max_amp;
    success_norm{ic} = success(ic).good_event_dfoverf./max_amp;
    
    spont_norm_avg(:,ic) = squeeze(mean(spont_norm{ic},1));
    release_norm_avg(:,ic) = squeeze(mean(release_norm{ic},1));
    press_norm_avg(:,ic) = squeeze(mean(press_norm{ic},1));
    success_norm_avg(:,ic) = squeeze(mean(success_norm{ic},1));
    fail_norm_avg(:,ic) = squeeze(mean(fail_norm{ic},1));
    
    spont_norm_sem(:,ic) = squeeze(std(spont_norm{ic},[],1))./sqrt(size(spont_norm{ic},1));
    release_norm_sem(:,ic) = squeeze(std(release_norm{ic},[],1))./sqrt(size(release_norm{ic},1));
    press_norm_sem(:,ic) = squeeze(std(press_norm{ic},[],1))./sqrt(size(press_norm{ic},1));
    errorbar(spont_norm_avg(:,ic), spont_norm_sem(:,ic), 'b')
    hold on
    errorbar(release_norm_avg(:,ic), release_norm_sem(:,ic), 'k')
    hold on
    errorbar(press_norm_avg(:,ic), press_norm_sem(:,ic), 'c')
    title([num2str(sum(~isnan(spont_norm_avg(1,:)),2)) ' spont; ' num2str(sum(~isnan(release_norm_avg(1,:)),2)) ' release; ' num2str(sum(~isnan(press_norm_avg(1,:)),2)) ' press'])
end
supertitle(['Average good events- norm to average Spont'])
saveas(fig, [dest_sub '_good_event_releaseVspontVpress.fig']);
print([dest_sub '_good_event_releaseVspontVpress.eps'], '-depsc');
print([dest_sub '_good_event_releaseVspontVpress.pdf'], '-dpdf');

%average all cells
release_renorm_avg = bsxfun(@rdivide, release_norm_avg, max(spont_norm_avg,[],1));
spont_renorm_avg = bsxfun(@rdivide, spont_norm_avg, max(spont_norm_avg,[],1));
press_renorm_avg = bsxfun(@rdivide, press_norm_avg, max(spont_norm_avg,[],1));
fig=figure;
subplot(1,2,1)
errorbar(nanmean(release_renorm_avg,2), nanstd(release_renorm_avg,[],2)./sqrt(sum(~isnan(release_renorm_avg(1,:)),2)),'k')
hold on
errorbar(nanmean(spont_renorm_avg,2), nanstd(spont_renorm_avg,[],2)./sqrt(sum(~isnan(spont_renorm_avg(1,:)),2)),'b')
hold on
errorbar(nanmean(press_renorm_avg,2), nanstd(press_renorm_avg,[],2)./sqrt(sum(~isnan(press_renorm_avg(1,:)),2)),'c')
title(['All cells- ' num2str(sum(~isnan(spont_norm_avg(1,:)),2)) ' spont; ' num2str(sum(~isnan(release_norm_avg(1,:)),2)) ' release; ' num2str(sum(~isnan(press_norm_avg(1,:)),2)) ' press'])

subplot(1,2,2)
errorbar(nanmean(release_renorm_avg(:,find(resp_h)),2), nanstd(release_renorm_avg(:,find(resp_h)),[],2)./sqrt(sum(~isnan(release_renorm_avg(1,find(resp_h))),2)),'k')
hold on
errorbar(nanmean(spont_renorm_avg(:,find(resp_h)),2), nanstd(spont_renorm_avg(:,find(resp_h)),[],2)./sqrt(sum(~isnan(spont_renorm_avg(1,find(resp_h))),2)),'b')
hold on
errorbar(nanmean(press_renorm_avg(:,find(resp_h)),2), nanstd(press_renorm_avg(:,find(resp_h)),[],2)./sqrt(sum(~isnan(press_renorm_avg(1,find(resp_h))),2)),'c')
title(['Responsive cells- ' num2str(sum(~isnan(spont_norm_avg(1,find(resp_h))),2)) ' spont; ' num2str(sum(~isnan(release_norm_avg(1,find(resp_h))),2)) ' release; ' num2str(sum(~isnan(press_norm_avg(1,find(resp_h))),2)) ' press'])

supertitle(['Average good events- Spont and Evoked- norm to average Spont'])
saveas(fig, [dest_sub '_avg_releaseVspontVpress.fig']);
print([dest_sub '_avg_releaseVspontVpress.eps'], '-depsc');
print([dest_sub '_avg_releaseVspontVpress.pdf'], '-dpdf');

save([dest_sub '_norm2spont.mat'], 'release_norm_avg', 'spont_norm_avg', 'press_norm_avg', 'success_norm_avg', 'fail_norm_avg');


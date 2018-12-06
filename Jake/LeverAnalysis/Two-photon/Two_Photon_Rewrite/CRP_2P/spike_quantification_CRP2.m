% quantifies amplitude and timing of events for each condition
% 1. calculate average timecourses for NR/OR/UR events
% 2. compare timecourse of NR/OR/UR events
% 3. compare PSTHs for NR/OR/UR 
% 4. compare amplitudes of evoked and spontaneous events

%load data
load([dest_sub '_spont_events.mat']);
load([dest_sub '_evoked_events.mat']);
load([dest 'parse_behavior.mat']);
load([dest '_pvals.mat']);
load([dest '_cell_categories.mat']);
load([dest '_cue_movies_lick.mat']); 
load([dest '_cue_movies.mat']); 
%load([dest 'ROI_TCs.mat']);
pre_rew_frames = pre_cue_frames + round((input.RewardDelayDurationMs + mean(cell2mat(input.reactTimesMs)))/ifi);

%define useful variables
nCells = size(data_tc_spont,2);
RS_cells = allresp_cells;
n = ceil(sqrt(nCells));
if nCells <((n.^2)-n)
    n2= n-1;
else
    n2 = n;
end

%% 1- calculate average events for different conditions 
for ic = 1:nCells
    %get mean and sem across trials for each cell
    events(ic).dfoverf_avg = nanmean(events(ic).dfoverf_chunk,1);
    events(ic).dfoverf_sem = std(events(ic).dfoverf_chunk,[],1)./sqrt(size(events(ic).dfoverf_chunk,1));
    
    %Reward driven events 
    normalR(ic).good_event_dfoverf = normalR(ic).dfoverf_chunk(find(normalR(ic).good_event==1),:);
    if ~isempty(omitR)
        omitR(ic).good_event_dfoverf = omitR(ic).dfoverf_chunk(find(omitR(ic).good_event==1),:);
    end
    
    %Cue driven events 
    normalCue(ic).good_event_dfoverf = normalCue(ic).dfoverf_chunk(find(normalCue(ic).good_event==1),:);
    if ~isempty(omitCue)
        omitCue(ic).good_event_dfoverf = omitCue(ic).dfoverf_chunk(find(omitCue(ic).good_event==1),:);
    end
    
    %mean reward driven events
    if isempty( normalR(ic).dfoverf_chunk(find(normalR(ic).good_event==1),:) )
        normalR(ic).good_event_dfoverf_avg = zeros(size(nanmean(normalR(ic).good_event_dfoverf,1)));
    else
        normalR(ic).good_event_dfoverf_avg = nanmean(normalR(ic).good_event_dfoverf,1);
    end
    
    %mean cue driven events
    if isempty( normalCue(ic).dfoverf_chunk(find(normalCue(ic).good_event==1),:) )
        normalCue(ic).good_event_dfoverf_avg = zeros(size(nanmean(normalCue(ic).good_event_dfoverf,1)));
    else
        normalCue(ic).good_event_dfoverf_avg = nanmean(normalCue(ic).good_event_dfoverf,1);
    end
    
    %calculate sem for reward and cue driven events
    normalR(ic).good_event_dfoverf_sem = std(normalR(ic).good_event_dfoverf,[],1)./sqrt(size(normalR(ic).good_event_dfoverf,1));
    normalCue(ic).good_event_dfoverf_sem = std(normalCue(ic).good_event_dfoverf,[],1)./sqrt(size(normalCue(ic).good_event_dfoverf,1));
    %define event peak as the frame data_start+100ms into the chunk????????????????????????????????????????????????????
    events(ic).good_event_dfoverf_peaks = events(ic).dfoverf_chunk(:,data_start+ceil(100/double(ifi))); 
    %store the number of each event type
    normalCue_nevents = size(normalCue(ic).good_event_dfoverf,1);
    normalR_nevents = size(normalR(ic).good_event_dfoverf,1);
    spont_nevents = size(events(ic).dfoverf_chunk,1);
end
%  save([dest_sub '_event_summary.mat'], 'normalR', 'normalCue', 'events', 'omitR', 'omitCue', 'normalCue_nevents', 'normalR_nevents', 'spont_nevents');

%% 2) plots - compare timecourse of NR/OR/UR events?
ts = [-pre_buffer_cue:post_buffer_cue].*double(ifi);

% %% compare timecourse of aligned good events- evoked
% %overlay individual events

%by cells
% fig=figure; 
% for ic = 1:nCells
%     subplot(n,n2,ic)
%     errorbar(normalR(ic).good_event_dfoverf_avg, normalR(ic).good_event_dfoverf_sem, 'k')
%     hold on
%     errorbar( normalCue(ic).good_event_dfoverf_avg, normalCue(ic).good_event_dfoverf_sem, 'g')
%     hold on
%     errorbar( events(ic).dfoverf_avg, events(ic).dfoverf_sem, 'm')
%     xlabel('Time (ms)')
%     ylabel('dF/F')
% %     if release_h(ic)
% %         rel_col = '\color{red}';
% %     else
% %         rel_col = '\color{black}';
% %     end
% %     if press_h(ic)
% %         pr_col = '\color{red}';
% %     else
% %         pr_col = '\color{black}';
% %     end
% %     title([rel_col num2str(size(success(ic).good_event_dfoverf,1)) ' S; ' rel_col num2str(size(fail(ic).good_event_dfoverf,1)) ' F; ' pr_col num2str(size(press(ic).good_event_dfoverf,1)) ' P'])
% end
% supertitle([mouse ' ' date '- Good events- Black: Reward; Green: Cue; Purple: Spont events'])
% saveas(fig, [dest_sub '_good_event_dfoverf.fig']);
% print([dest_sub '_good_event_dfoverf.eps'], '-depsc');
% print([dest_sub '_good_event_dfoverf.pdf'], '-dpdf');

resp_h = find(allresp_cells);

%% 3) compare PSTHs for NR/OR/UR 
% fig=figure;
% plot(repmat(spike_frames,1,2)', bsxfun(@times,[ie:ie+1],ones(length(spike_frames),2))', 'k')

sz1 = size(NRCue_tc, 1); % #frames
sz1_half = floor(sz1/2);
sz_n = size(normalCue(1).f_chunk,1);% #events
if ~isempty(omitCue)
    sz_o = size(omitCue(1).f_chunk,1);
end
if ~isempty(unexpCue)
    sz_u = size(unexpCue(1).f_chunk,1);
end

all_NR_hist = zeros(sz1,nCells);
all_OR_hist = zeros(sz1,nCells);
all_UR_hist = zeros(sz1,nCells);
all_NR_nolick_hist = zeros(sz1,nCells);
all_OR_nolick_hist = zeros(sz1,nCells);
all_UR_nolick_hist = zeros(sz1,nCells);
NegRate_cell = [];

for ic = 1:nCells
    %for each cell convert spike indeces to histograms 
    normalCue(ic).hist = zeros(sz1,sz_n);
    for is = 1:sz_n
        spike_frames = normalCue(ic).ind{is};
        normalCue(ic).hist(spike_frames,is) = 1;
    end
    if ~isempty(omitCue)
        omitCue(ic).hist = zeros(sz1,sz_o);
        for ie = 1:sz_o
            spike_frames = omitCue(ic).ind{ie};
            omitCue(ic).hist(spike_frames,ie) = 1;
        end
    end
    if ~isempty(unexpCue)
        unexpCue(ic).hist = zeros(sz1,sz_u);
        for ie = 1:sz_u
            spike_frames = unexpCue(ic).ind{ie};
            unexpCue(ic).hist(spike_frames,ie) = 1;
        end
    end
    
    %exclude frames from visual stimulus
    if strcmp(session, '180502_img084') | strcmp(session, '180403_img077') | strcmp(session, '180404_img077')
        if strcmp(session, '180502_img084')
            artifact_ind = [2:6];
        elseif strcmp(session, '180403_img077')
            artifact_ind = [2:5];
        elseif strcmp(session, '180404_img077')
            artifact_ind = [3:6];
        end
        normalCue(ic).hist([pre_cue_frames+1+artifact_ind],:) = NaN;
        omitCue(ic).hist([pre_cue_frames+1+artifact_ind],:) = NaN;
    end
    
    %calculate mean histogram matrices for each condition
    %normal reward
    all_NR_hist(:,ic) = mean(normalCue(ic).hist,2);
    %unexpected reward
    if ~isempty(unexpCue)
        all_UR_hist(:,ic) = mean(unexpCue(ic).hist,2);
    end
    %omitted reward (also looks for negatively modulated neurons)
    if ~isempty(omitCue)
        all_OR_hist(:,ic) = mean(omitCue(ic).hist,2);
        all_OR_hist_baseRate(:,ic) = all_OR_hist(pre_cue_frames - 6 : pre_cue_frames, ic);
        all_OR_hist_NegRate(:,ic) = all_OR_hist(pre_rew_frames : pre_rew_frames + 12, ic);
        
        diff_NegRate = diff(all_OR_hist_NegRate(:,ic));
        decrease_idx = find(diff_NegRate < 0, 1, 'first');
        
        if ( diff_NegRate(2) < 0) && ( mean(all_OR_hist_NegRate(1:2,ic),1) < mean(all_OR_hist_baseRate(:,ic),1) )
            NegRate_cell = [NegRate_cell ic];
        end
%         v = all_OR_hist_NegRate(:,ic);
%         [maxtab, mintab] = peakdet(v, 0.001);
%         if ~isempty(mintab) && ~isempty(maxtab)
%             max_peak = maxtab(maxtab(:,2) == max(maxtab(:,2)),:);
%             max_peak_idx = max_peak(1,1);
%             max_peak = max_peak(1,2);
%
%             min_peak = mintab(mintab(:,2) == min(mintab(:,2)),:);
%             min_peak_idx = min_peak(1,1);
%             min_peak = min_peak(1,2);
%             
%             if min_peak < mean(all_OR_hist_baseRate(:,ic)) && max_peak_idx > min_peak_idx
%                 NegRate_cell = [NegRate_cell ic];
%             end
%         end
    end

    %calculate nolick condition histograms
    all_NR_nolick_hist(:,ic) = mean(normalCue(ic).hist(:,logical(NR_lick_info.no_lick_cue_to_500)),2);  %logical(NR_lick_info.no_lick_cue_to_500)   %~NR_lick_info.lickTrial
    if ~isempty(omitCue)
        all_OR_nolick_hist(:,ic) = mean(omitCue(ic).hist(:,logical(OR_lick_info.no_lick_cue_to_500)),2);
    end
    if ~isempty(unexpCue)
        all_UR_nolick_hist(:,ic) = mean(unexpCue(ic).hist(:,logical(UR_lick_info.no_lick_cue_to_500)),2);
    end
%     all_NR_syn(:,:,ic) = normalR(ic).hist;
%     all_OR_syn(:,:,ic) = omitR(ic).hist;
end

% save([dest_sub '_evoked_events.mat'], 'NR_tc', 'OR_tc', 'UR_tc', 'NRCue_tc', 'ORCue_tc', 'URCue_tc', 'normalR', 'omitR', 'unexpR', 'normalCue', 'omitCue', 'unexpCue', 'min_iei', 'opt', 'pre_buffer_cue', 'post_buffer_cue',...
%     'NRdata_start_rew', 'NRdata_end_rew', 'Cuedata_start', 'Cuedata_end', 'post_buffer_rew', 'pre_buffer_rew')

fig=figure;
% errorbar(ts,nanmean(all_NR_hist.*(1000/double(ifi)),2),nanstd(all_NR_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_NR_hist,2)),'k')
plot(ts,nanmean(all_NR_hist.*(1000/double(ifi)),2), 'k');
if ~isempty(omitCue)
    hold on
%     errorbar(ts,nanmean(all_OR_hist.*(1000/double(ifi)),2),nanstd(all_OR_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_OR_hist,2)),'r')
    plot(ts, nanmean(all_OR_hist.*(1000/double(ifi)),2), 'r');
end

if ~isempty(unexpCue)
    hold on
%     errorbar(ts,nanmean(all_UR_hist.*(1000/double(ifi)),2),nanstd(all_UR_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_UR_hist,2)),'g')
    plot(ts, nanmean(all_UR_hist.*(1000/double(ifi)),2), 'g');
end

title('All NR (black) and all OR (red) PSTHs. green=UR');
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');

% saveas(fig, [dest_sub '_all_event_avgPSTH.fig']);
% print([dest_sub '_all_event_avgPSTH.eps'], '-depsc');
% print([dest_sub '_all_event_avgPSTH.pdf'], '-dpdf');

% errorbar(ts,nanmean(all_NR_hist.*(1000/double(ifi)),2),nanstd(all_NR_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_NR_hist,2)),'k')

if ~isempty(omitCue) && ~isempty(NegRate_cell)
    fig = figure;
    
    nn = ceil(sqrt(length(NegRate_cell)));
    if length(NegRate_cell) <((nn.^2)-nn)
        nn2= nn-1;
    else
        nn2 = nn;
    end
    for ic = 1:length(NegRate_cell)
        subplot(nn,nn2,ic)
        %     errorbar(tt,avg_NR(ic,:), sem_NR(ic,:),'k');hold on;
        %     errorbar(tt,avg_OR(ic,:), sem_OR(ic,:),'r');
        %     errorbar(tt,avg_UR(ic,:), sem_UR(ic,:),'b');
        plot(ts,all_OR_hist(:,NegRate_cell(ic)).*(1000/double(ifi)),'r');hold on;
        hold on; vline((pre_rew_frames - pre_buffer_cue)*33, '--k');
        ylim([0 4])
        xlim([500 1000])
        %     set(gca,'XTick',[)
    end
    suptitle([session_date ' ' mouse_ID ' Average cue resp: Normal- black (n = ' num2str(size(NR_movie,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie,1)) ' trials)'])
    orient landscape
    
    fig=figure;
    hold on
    %     errorbar(ts,nanmean(all_OR_hist.*(1000/double(ifi)),2),nanstd(all_OR_hist.*(1000/double(ifi)),[],2)./sqrt(size(all_OR_hist,2)),'r')
    plot(ts, nanmean(all_OR_hist(:, NegRate_cell).*(1000/double(ifi)),2), 'r');
    vline((pre_rew_frames - pre_buffer_cue)*33, '--k');
    title(['nCells with decrease firing rate = ', num2str(length(NegRate_cell))]);
    xlabel('Time (ms)')
    ylabel('Firing rate (Hz)')
    saveas(fig, [dest_sub '_neg_event_avgPSTH.fig']);
    print([dest_sub '_neg_event_avgPSTH.eps'], '-depsc');
    print([dest_sub '_neg_event_avgPSTH.pdf'], '-dpdf');
    all_OR_hist_negR = all_OR_hist(:, NegRate_cell);
else
    all_OR_hist_negR = [];
end

% find peak and latency of PSTH
% fig=figure;
NR_ind_all = [];
NR_rate_all = [];
OR_ind_all = [];
OR_rate_all = [];
for ic = 1:nCells
    [NR_rate, NR_ind_all(ic)] = max(mean(normalCue(ic).hist((sz1_half):(sz1_half)+ceil(500/double(ifi)),:),2),[],1);
%     basef = nanmean(mean(fail(ic).hist(1:sz1_half-ceil(65/double(ifi)),:),1),2);
    if ~isempty(omitR)
        [OR_rate, OR_ind_all(ic)] = max(mean(omitCue(ic).hist((sz1_half):(sz1_half)+ceil(500/double(ifi)),:),2),[],1);
    end
end

NR_ind_RS = [];
OR_ind_RS = [];
start = 1;
OR_rate = [];
for ic = find(resp_h)
    [NR_rate, NR_ind_RS(start)] = max(mean(normalCue(ic).hist((sz1_half):(sz1_half)+ceil(500/double(ifi)),:),2),[],1);
    if ~isempty(omitR)
        [OR_rate, OR_ind_RS(start)] = max(mean(omitCue(ic).hist((sz1_half):(sz1_half)+ceil(500/double(ifi)),:),2),[],1);
    end
    start = start+1;
end

% save([dest_sub '_event_hist.mat'], 'NR_ind_all', 'OR_ind_all', 'all_NR_hist', 'all_OR_hist', 'all_UR_hist', 'all_NR_nolick_hist', ...
%     'all_OR_nolick_hist', 'all_UR_nolick_hist', 'NR_ind_RS', 'OR_ind_RS', 'all_OR_hist_negR');

%% 4) Identify responsive cells based on PSTHs

%POSITIVE CUE responsive neurons 
base_cue_window = [pre_cue_frames-round(500/double(ifi)):pre_cue_frames];
resp_cue_window = [pre_cue_frames+1+round(100/double(ifi)):pre_cue_frames+1+round(600/double(ifi))];
effect_sign = 'pos';
[NR_Cue_h_pos, NR_Cue_p_pos, NR_Cue_resp_cells_pos, NR_Cue_resp_avg_pos, NR_Cue_resp_sem_pos, NR_Cue_base_pos, NR_Cue_resp_pos] = findRespCell_spike(normalCue, base_cue_window, resp_cue_window, effect_sign);
[OR_Cue_h_pos, OR_Cue_p_pos, OR_Cue_resp_cells_pos, OR_Cue_resp_avg_pos, OR_Cue_resp_sem_pos, OR_Cue_base_pos, OR_Cue_resp_pos] = findRespCell_spike(omitCue, base_cue_window, resp_cue_window, effect_sign);
[UR_Cue_h_pos, UR_Cue_p_pos, UR_Cue_resp_cells_pos, UR_Cue_resp_avg_pos, UR_Cue_resp_sem_pos, UR_Cue_base_pos, UR_Cue_resp_pos] = findRespCell_spike(unexpCue, base_cue_window, resp_cue_window, effect_sign);

%NEGATIVE CUE responsive neurons
effect_sign = 'neg';
[NR_Cue_h_neg, NR_Cue_p_neg, NR_Cue_resp_cells_neg, NR_Cue_resp_avg_neg, NR_Cue_resp_sem_neg, NR_Cue_base_neg, NR_Cue_resp_neg] = findRespCell_spike(normalCue, base_cue_window, resp_cue_window, effect_sign);
[OR_Cue_h_neg, OR_Cue_p_neg, OR_Cue_resp_cells_neg, OR_Cue_resp_avg_neg, OR_Cue_resp_sem_neg, OR_Cue_base_neg, OR_Cue_resp_neg] = findRespCell_spike(omitCue, base_cue_window, resp_cue_window, effect_sign);
[UR_Cue_h_neg, UR_Cue_p_neg, UR_Cue_resp_cells_neg, UR_Cue_resp_avg_neg, UR_Cue_resp_sem_neg, UR_Cue_base_neg, UR_Cue_resp_neg] = findRespCell_spike(unexpCue, base_cue_window, resp_cue_window, effect_sign);

%POSITIVE REWARD responsive neurons
pre_rew_frames = pre_cue_frames + round((trial_outcome.normalReward(1) - trial_outcome.normalRewardCue(1))/ifi);
base_reward_window = base_cue_window;
resp_reward_window=[ pre_rew_frames +round(100/double(ifi)) : pre_rew_frames +round(600/double(ifi))];
effect_sign = 'pos';
[NR_Rew_h_pos, NR_Rew_p_pos, NR_Rew_resp_cells_pos, NR_Rew_resp_avg_pos, NR_Rew_resp_sem_pos, NR_Rew_base_pos, NR_Rew_resp_pos] = findRespCell_spike(normalCue, base_reward_window, resp_reward_window, effect_sign);
[OR_Rew_h_pos, OR_Rew_p_pos, OR_Rew_resp_cells_pos, OR_Rew_resp_avg_pos, OR_Rew_resp_sem_pos, OR_Rew_base_pos, OR_Rew_resp_pos] = findRespCell_spike(omitCue, base_reward_window, resp_reward_window, effect_sign);
[UR_Rew_h_pos, UR_Rew_p_pos, UR_Rew_resp_cells_pos, UR_Rew_resp_avg_pos, UR_Rew_resp_sem_pos, UR_Rew_base_pos, UR_Rew_resp_pos] = findRespCell_spike(unexpCue, base_reward_window, resp_reward_window, effect_sign);

%NEGATIVE REWARD responsive neurons
effect_sign = 'neg';
[NR_Rew_h_neg, NR_Rew_p_neg, NR_Rew_resp_cells_neg, NR_Rew_resp_avg_neg, NR_Rew_resp_sem_neg, NR_Rew_base_neg, NR_Rew_resp_neg] = findRespCell_spike(normalCue, base_reward_window, resp_reward_window, effect_sign);
[OR_Rew_h_neg, OR_Rew_p_neg, OR_Rew_resp_cells_neg, OR_Rew_resp_avg_neg, OR_Rew_resp_sem_neg, OR_Rew_base_neg, OR_Rew_resp_neg] = findRespCell_spike(omitCue, base_reward_window, resp_reward_window, effect_sign);
[UR_Rew_h_neg, UR_Rew_p_neg, UR_Rew_resp_cells_neg, UR_Rew_resp_avg_neg, UR_Rew_resp_sem_neg, UR_Rew_base_neg, UR_Rew_resp_neg] = findRespCell_spike(unexpCue, base_reward_window, resp_reward_window, effect_sign);

% save([dest_sub '_cell_resp_spike.mat'], 'NR_Cue_base_pos', 'NR_Cue_resp_pos', 'OR_Cue_base_pos',  'OR_Cue_resp_pos', 'UR_Cue_base_pos', 'UR_Cue_resp_pos',...
%     'NR_Rew_base_pos', 'NR_Rew_resp_pos', 'OR_Rew_base_pos', 'OR_Rew_resp_pos', 'UR_Rew_base_pos', 'UR_Rew_resp_pos', ...
% 	'NR_Cue_base_neg', 'NR_Cue_resp_neg', 'OR_Cue_base_neg', 'OR_Cue_resp_neg', 'UR_Cue_base_neg', 'UR_Cue_resp_neg', ...
%     'NR_Rew_base_neg', 'NR_Rew_resp_neg', 'OR_Rew_base_neg', 'OR_Rew_resp_neg', 'UR_Rew_base_neg', 'UR_Rew_resp_neg');

%Make combinations of responsive subtypes for plotting later
if ~isempty(OR_movie)
    allresp_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos | NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos | NR_Rew_h_pos | OR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg | NR_Rew_h_neg | OR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | OR_Cue_h_pos | NR_Cue_h_neg | OR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos | OR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg | OR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | OR_Rew_h_pos | NR_Rew_h_neg | OR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos | OR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg | OR_Rew_h_neg;
else
    allresp_cells = NR_Cue_h_pos | NR_Rew_h_pos | NR_Cue_h_neg | NR_Rew_h_neg;
    allresp_cells_pos = NR_Cue_h_pos | NR_Rew_h_pos;
    allresp_cells_neg = NR_Cue_h_neg | NR_Rew_h_neg;
    cue_cells = NR_Cue_h_pos | NR_Cue_h_neg;
    cue_cells_pos = NR_Cue_h_pos;
    cue_cells_neg = NR_Cue_h_neg;
    rew_cells = NR_Rew_h_pos | NR_Rew_h_neg;
    rew_cells_pos = NR_Rew_h_pos;
    rew_cells_neg = NR_Rew_h_neg;
end

if ~isempty(UR_movie)
    allresp_cells = allresp_cells | UR_Cue_h_pos | UR_Rew_h_pos | UR_Cue_h_neg | UR_Rew_h_neg;
    allresp_cells_pos = allresp_cells_pos | UR_Cue_h_pos | UR_Rew_h_pos;
    allresp_cells_neg = allresp_cells_neg | UR_Cue_h_neg | UR_Rew_h_neg;
    cue_cells = cue_cells | UR_Cue_h_pos | UR_Cue_h_neg;
    cue_cells_pos = cue_cells_pos | UR_Cue_h_pos;
    cue_cells_neg = cue_cells_neg | UR_Cue_h_neg;
    rew_cells = rew_cells | UR_Rew_h_pos | UR_Rew_h_neg;
    rew_cells_pos = rew_cells_pos | UR_Rew_h_pos;
    rew_cells_neg = rew_cells_neg | UR_Rew_h_neg;
end

cue_cells_perc = sum(cue_cells) / nCells;
cue_cells_perc_pos = sum(cue_cells_pos) / nCells;
cue_cells_perc_neg = sum(cue_cells_neg) / nCells;
rew_cells_perc = sum(rew_cells) / nCells;
rew_cells_perc_pos = sum(rew_cells_pos) / nCells;
rew_cells_perc_neg = sum(rew_cells_neg) / nCells;

% save([dest_sub '_cell_categories_spike.mat'], 'NR_Cue_resp_cells_pos', 'NR_Cue_resp_cells_neg', 'OR_Cue_resp_cells_pos', 'OR_Cue_resp_cells_neg', ...  %'OR_Cue_resp_cells2_pos', 'OR_Cue_resp_cells2_neg',
%     'UR_Cue_resp_cells_pos', 'UR_Cue_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg', 'OR_Rew_resp_cells_pos', 'OR_Rew_resp_cells_neg', 'NR_Rew_resp_cells_pos', 'NR_Rew_resp_cells_neg',...
%     'allresp_cells', 'allresp_cells_pos', 'allresp_cells_neg', 'nCells', 'cue_cells', 'cue_cells_pos', 'cue_cells_neg',  'rew_cells', 'rew_cells_pos', 'rew_cells_neg', ...
%     'cue_cells_perc', 'cue_cells_perc_pos', 'cue_cells_perc_neg', 'rew_cells_perc', 'rew_cells_perc_pos', 'rew_cells_perc_neg');
% save([dest_sub '_pvals_spike.mat'], 'NR_Cue_p_pos', 'NR_Cue_p_neg', 'NR_Cue_h_pos', 'NR_Cue_h_neg', 'NR_Rew_p_pos', 'NR_Rew_p_neg', 'NR_Rew_h_pos', 'NR_Rew_h_neg', ...
%     'OR_Cue_p_pos', 'OR_Cue_p_neg', 'OR_Cue_h_pos', 'OR_Cue_h_neg', 'OR_Rew_p_pos', 'OR_Rew_p_pos', 'OR_Rew_h_pos', 'OR_Rew_h_neg', ... %'OR_Cue_h2_pos', 'OR_Cue_h2_neg', 'OR_Cue_p2_pos', 'OR_Cue_p2_neg', 
%     'UR_Cue_p_pos', 'UR_Cue_p_neg', 'UR_Cue_h_pos', 'UR_Cue_h_neg', 'UR_Rew_p_pos', 'UR_Rew_p_neg', 'UR_Rew_h_pos', 'UR_Rew_h_neg');

%% 5) ID which trials have a Ca event after the rew
if session_subset == 2
    num_NR = length(normalCue(1).ind);
    num_OR = length(omitCue(1).ind);
    rew_event_win = pre_rew_frames+1+2 : pre_rew_frames+1+2+(400/double(ifi));
    NR_rew_event_trials = zeros(length(normalCue), num_NR);
    OR_rew_event_trials = zeros(length(normalCue), num_OR);
    for cell_num = 1:length(normalCue)
        for trial_num = 1:num_NR
            this_trial_ind = normalCue(cell_num).ind{trial_num}';
            NR_rew_event_trials(cell_num, trial_num) = length(find(this_trial_ind >= rew_event_win(1) & this_trial_ind <= rew_event_win(end)));
            if trial_num <= num_OR
                this_trial_ind = omitCue(cell_num).ind{trial_num}';
                OR_rew_event_trials(cell_num, trial_num) = length(find(this_trial_ind >= rew_event_win(1) & this_trial_ind <= rew_event_win(end)));
            end
        end
    end
    save([dest_sub '_rew_event_trials.mat'], 'NR_rew_event_trials', 'OR_rew_event_trials', 'rew_event_win');
    
    %look for trials with a Ca event after the cue
    cue_event_win =  pre_cue_frames+5:pre_rew_frames-1;
    NR_cue_event_trials = zeros(length(normalCue), num_NR);
    OR_cue_event_trials = zeros(length(normalCue), num_OR);
    for cell_num = 1:length(normalCue)
        for trial_num = 1:num_NR
            this_trial_ind = normalCue(cell_num).ind{trial_num}';
            NR_cue_event_trials(cell_num, trial_num) = length(find(this_trial_ind >= cue_event_win(1) & this_trial_ind <= cue_event_win(end)));
            if trial_num <= num_OR
                this_trial_ind = omitCue(cell_num).ind{trial_num}';
                OR_cue_event_trials(cell_num, trial_num) = length(find(this_trial_ind >= cue_event_win(1) & this_trial_ind <= cue_event_win(end)));
            end
        end
    end
    %save([dest_sub '_cue_event_trials.mat'], 'NR_cue_event_trials', 'OR_cue_event_trials', 'cue_event_win');
end
% 
% %% 5) compare spontaneous and evoked event amplitudes and waveforms
% 
% %compare evoked and spontaneous waveforms
% fig=figure('rend', 'painters', 'pos', [50 100 1200 700]);
% spont_norm_avg = zeros(size(events(1).dfoverf_chunk,2),nCells);
% spont_norm_sem = zeros(size(events(1).dfoverf_chunk,2),nCells);
% normalR_norm_avg = zeros(size(normalR(1).dfoverf_chunk,2),nCells);
% normalCue_norm_avg = zeros(size(normalCue(1).dfoverf_chunk,2),nCells);
% normalR_norm_sem = zeros(size(normalR(1).dfoverf_chunk,2),nCells);
% normalCue_norm_sem = zeros(size(normalCue(1).dfoverf_chunk,2),nCells);
% if ~isempty(omitR)
%     omitR_norm_avg = zeros(size(omitR(1).dfoverf_chunk,2),nCells);
%     omitCue_norm_avg = zeros(size(omitCue(1).dfoverf_chunk,2),nCells);
%     omitR_norm_sem = zeros(size(omitR(1).dfoverf_chunk,2),nCells);
%     omitCue_norm_sem = zeros(size(omitCue(1).dfoverf_chunk,2),nCells);
% end
% 
% if ~isempty(omitR)
%     for ic = 1:nCells
%         subplot(n,n2,ic)
%         %      max_amp = nanmean(events(ic).good_event_dfoverf_peaks,1); %normalizing to the average spont, not max
%         max_amp = 1;
%         spont_norm1{ic} = events(ic).dfoverf_chunk./max_amp;
%         normalR_norm1{ic} = normalR(ic).good_event_dfoverf./max_amp;
%         omitR_norm1{ic} = omitR(ic).good_event_dfoverf./max_amp;
%         normalCue_norm1{ic} = normalCue(ic).good_event_dfoverf./max_amp;
%         omitCue_norm1{ic} = omitCue(ic).good_event_dfoverf./max_amp;
%         
%         spont_norm1_avg(:,ic) = squeeze(mean(spont_norm1{ic},1));
%         normalR_norm1_avg(:,ic) = squeeze(mean(normalR_norm1{ic},1));
%         omitR_norm1_avg(:,ic) = squeeze(mean(omitR_norm1{ic},1));
%         normalCue_norm1_avg(:,ic) = squeeze(mean(normalCue_norm1{ic},1));
%         omitCue_norm1_avg(:,ic) = squeeze(mean(omitCue_norm1{ic},1));
%         
%         spont_norm1_sem(:,ic) = squeeze(std(spont_norm1{ic},[],1))./sqrt(size(spont_norm1{ic},1));
%         normalR_norm1_sem(:,ic) = squeeze(std(normalR_norm1{ic},[],1))./sqrt(size(normalR_norm1{ic},1));
%         omitR_norm1_sem(:,ic) = squeeze(std(omitR_norm1{ic},[],1))./sqrt(size(omitR_norm1{ic},1));
%         normalCue_norm1_sem(:,ic) = squeeze(std(normalCue_norm1{ic},[],1))./sqrt(size(normalCue_norm1{ic},1));
%         omitCue_norm1_sem(:,ic) = squeeze(std(omitCue_norm1{ic},[],1))./sqrt(size(omitCue_norm1{ic},1));
%         
%         %normalize to average spont
%         max_amp = nanmean(events(ic).good_event_dfoverf_peaks,1); %normalizing to the average spont, not max
%         spont_norm{ic} = events(ic).dfoverf_chunk./max_amp;
%         normalR_norm{ic} = normalR(ic).good_event_dfoverf./max_amp;
%         omitR_norm{ic} = omitR(ic).good_event_dfoverf./max_amp;
%         normalCue_norm{ic} = normalCue(ic).good_event_dfoverf./max_amp;
%         omitCue_norm{ic} = omitCue(ic).good_event_dfoverf./max_amp;
%         
%         spont_norm_avg(:,ic) = squeeze(mean(spont_norm{ic},1));
%         normalR_norm_avg(:,ic) = squeeze(mean(normalR_norm{ic},1));
%         omitR_norm_avg(:,ic) = squeeze(mean(omitR_norm{ic},1));
%         normalCue_norm_avg(:,ic) = squeeze(mean(normalCue_norm{ic},1));
%         omitCue_norm_avg(:,ic) = squeeze(mean(omitCue_norm{ic},1));
%         
%         spont_norm_sem(:,ic) = squeeze(std(spont_norm{ic},[],1))./sqrt(size(spont_norm{ic},1));
%         normalR_norm_sem(:,ic) = squeeze(std(normalR_norm{ic},[],1))./sqrt(size(normalR_norm{ic},1));
%         omitR_norm_sem(:,ic) = squeeze(std(omitR_norm{ic},[],1))./sqrt(size(omitR_norm{ic},1));
%         normalCue_norm_sem(:,ic) = squeeze(std(normalCue_norm{ic},[],1))./sqrt(size(normalCue_norm{ic},1));
%         omitCue_norm_sem(:,ic) = squeeze(std(omitCue_norm{ic},[],1))./sqrt(size(omitCue_norm{ic},1));
%         errorbar(spont_norm_avg(:,ic), spont_norm_sem(:,ic), 'b')
%         hold on
%         errorbar(normalR_norm_avg(:,ic), normalR_norm_sem(:,ic), 'k')
%         hold on
%         errorbar(normalCue_norm_avg(:,ic), normalCue_norm_sem(:,ic), 'c')
%     end
%     suptitle(['Blue: ', num2str(sum(~isnan(spont_norm_avg(:,1)))) ' spont;   Black: ' num2str(sum(~isnan(normalR_norm_avg(:,1)))) ' normalR;   Cyan:' num2str(sum(~isnan(normalCue_norm_avg(:,1)))) ' normalCue' ]);
%     savefig([dest_sub, 'cmp_event_TC_each_cell_spont_Rew_Cue']);
%     save([dest_sub '_norm2spont.mat'], 'spont_norm_avg', 'normalR_norm_avg', 'omitR_norm_avg', 'normalCue_norm_avg', 'omitCue_norm_avg',...
%         'spont_norm1_avg', 'normalR_norm1_avg', 'omitR_norm1_avg', 'normalCue_norm1_avg', 'omitCue_norm1_avg');
% end



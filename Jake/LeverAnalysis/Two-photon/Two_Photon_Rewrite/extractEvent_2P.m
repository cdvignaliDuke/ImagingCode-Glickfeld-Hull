function [success_tc, success] = extractEvent_2P(data_tc, trial_outcome, frame_info, pre_buffer, post_buffer, opt, ifi, data_start, data_end)
nIC = size(data_tc,2);
success_fr = zeros(length(trial_outcome),pre_buffer+post_buffer+1);

success_tc = [];

neg_ind = find(trial_outcome<0);
if isempty(neg_ind)
    start_is = 1;
else
    start_is = neg_ind(end) + 1;
end
for is = start_is:length(trial_outcome)    %this could be problematic due to unremoved NaNs
    
    if trial_outcome(is) < length(frame_info.counter) && frame_info.counter(trial_outcome(is))-pre_buffer > 0
        success_fr(is,:) = [frame_info.counter(trial_outcome(is))-pre_buffer:frame_info.counter(trial_outcome(is))+post_buffer];
        if and(min(success_fr(is,:),[],2)>0, max(success_fr(is,:),[],2)<size(data_tc,1))
            success_tc = cat(3, success_tc, data_tc(success_fr(is,:),:));
        end
    elseif trial_outcome(is) >= length(frame_info.counter)
        success_fr(is:end,:) = [];
        break
    end
end
if ~isempty(neg_ind)
    success_fr(1:start_is-1 ,:)=[];
end
[xx, yy] = find(success_fr == 0);
success_fr(unique(xx),:) = [];

success = [];
win_release_start = round(300./double(ifi));
win_release_end = round(100./double(ifi));
min_iei = 650;
for ic = 1:nIC
    success_diff = diff(squeeze(success_tc(:,ic,:)),[],1);
    success(ic).f_chunk = [];
    success(ic).ind = {};
    success(ic).good_ind = {};
    for is = 1:size(success_diff,2)
        [~,ind,~] = CellsortFindspikes(success_diff(:,is), opt.thresh, opt.dt, opt.deconvtau, opt.normalization);
        %         ind = find(success_diff(:,is)>thresh(ic));
        ind(find(diff(ind)==1)+1)=[];
        double_events = find(diff(ind)<min_iei);
        ind_good = ind;
        ind_good([double_events (double_events+1)])=[];
        success(ic).ind{is} = ind;
        success(ic).good_ind{is} = ind_good;
        if intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))
            wind = success_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
            if intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end))
                success(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-win_release_start), find(ind_good<pre_buffer+win_release_end)));
                first_event = success_fr(is,success(ic).event_ind(:,is));
                success(ic).good_event(:,is) = 1;
                success(ic).event(:,is) = 1;
            else
                success(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer-win_release_start), find(ind<pre_buffer+win_release_end))),[],1);
                first_event = success_fr(is,success(ic).event_ind(:,is));
                success(ic).good_event(:,is) = 0;
                success(ic).event(:,is) = 1;
            end
            success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
        else
            success(ic).event_ind(:,is) = NaN;
            wind = success_fr(is,pre_buffer-win_release_start:pre_buffer+win_release_end);
            first_event = wind(1);
            success(ic).good_event(:,is) = NaN;
            success(ic).event(:,is) = 0;
            success(ic).f_chunk = [success(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; %data_tc(success_fr(is,:),ic)']; %
        end
    end
    success(ic).f_base = mean(success(ic).f_chunk(:,data_start-round(50./double(ifi)):data_start+1),2); %13:15),2); %
    success(ic).df_chunk = bsxfun(@minus, success(ic).f_chunk, success(ic).f_base);
    success(ic).dfoverf_chunk = bsxfun(@rdivide, success(ic).df_chunk, success(ic).f_base);
end

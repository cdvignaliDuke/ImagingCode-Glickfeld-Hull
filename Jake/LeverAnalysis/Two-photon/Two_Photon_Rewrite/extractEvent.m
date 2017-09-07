function [NR_tc, normalR] = extractEvent(data_tc, frame_info, trial_outcome, lick, opt, pre_buffer, post_buffer, data_start, data_end, min_iei)
NR_fr = zeros(length(trial_outcome),pre_buffer+post_buffer+1);
nIC = size(data_tc,2);
NR_tc = [];

neg_ind = find(trial_outcome<0);

if isempty(neg_ind)
    start_is = 1;
else
    start_is = neg_ind(end) + 1;
end

for is = start_is:length(trial_outcome)    
    if trial_outcome(is) < length(frame_info.counter) && is <= length(lick.lickTrial)
        NR_fr(is,:) = [frame_info.counter(trial_outcome(is))-pre_buffer:frame_info.counter(trial_outcome(is))+post_buffer];
        if and(min(NR_fr(is,:),[],2)>0, max(NR_fr(is,:),[],2)<size(data_tc,1))
            NR_tc = cat(3, NR_tc, data_tc(NR_fr(is,:),:));
        end
    else
        NR_fr(is:end,:) = [];
        break
    end
end

if ~isempty(neg_ind)
    NR_fr(1:start_is-1,:)=[];
end

%find all events during trials
normalR = [];
win_release_end = round(500./double(frame_info.ifi));
for ic = 1:nIC
    
    normalR_diff = diff(squeeze(NR_tc(:,ic,:)),[],1);
    normalR(ic).f_chunk = [];
    normalR(ic).ind = {};
    normalR(ic).good_ind = {};
    for is = 1:size(normalR_diff,2)
        
        [~,ind,~] = CellsortFindspikes(normalR_diff(:,is), opt.thresh, opt.dt, opt.deconvtau, opt.normalization);
        %         ind = find(success_diff(:,is)>thresh(ic));
        ind(find(diff(ind)==1)+1)=[];
        double_events = find(diff(ind)<min_iei);
        if ~isempty(double_events)
            ind_good = [];
        else
            ind_good = ind;
        end
%         ind_good([double_events (double_events+1)])=[];
        normalR(ic).ind{is} = ind;
        normalR(ic).good_ind{is} = ind_good;
        % exclude trials with licking bout
        if lick.lickTrial(is) ~= 1 && ~isempty(intersect(find(ind>pre_buffer), find(ind<pre_buffer+win_release_end))) % find events between onset and 500ms
            wind = NR_fr(is,pre_buffer:pre_buffer+win_release_end);
            if intersect(find(ind_good>pre_buffer), find(ind_good<pre_buffer+win_release_end))
                normalR(ic).event_ind(:,is) =ind_good(intersect(find(ind_good>pre_buffer-1), find(ind_good<pre_buffer+win_release_end)));
                first_event = NR_fr(is,normalR(ic).event_ind(:,is));
                normalR(ic).good_event(:,is) = 1;
                normalR(ic).event(:,is) = 1;
            else
                normalR(ic).event_ind(:,is) = min(ind(intersect(find(ind>pre_buffer), find(ind<pre_buffer+win_release_end))),[],1);
                first_event = NR_fr(is,normalR(ic).event_ind(:,is));
                normalR(ic).good_event(:,is) = 0;
                normalR(ic).event(:,is) = 1;
            end
            normalR(ic).f_chunk = [normalR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; 
        else
            normalR(ic).event_ind(:,is) = NaN;
            wind = NR_fr(is,pre_buffer:pre_buffer+1);
            first_event = wind(1);
            normalR(ic).good_event(:,is) = NaN;
            normalR(ic).event(:,is) = 0;
            normalR(ic).f_chunk = [normalR(ic).f_chunk; data_tc(first_event-data_start:first_event+data_end,ic)']; 
        end
    end
    normalR(ic).f_base = mean(normalR(ic).f_chunk(:,data_start-round(50./double(frame_info.ifi)):data_start+1),2); %13:15),2); %
    normalR(ic).df_chunk = bsxfun(@minus, normalR(ic).f_chunk, normalR(ic).f_base);
    normalR(ic).dfoverf_chunk = bsxfun(@rdivide, normalR(ic).df_chunk, normalR(ic).f_base);
end
end
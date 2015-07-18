nIC = size(data_tc,1);
n = ceil(sqrt(nIC));
if nIC <((n.^2)-n)
    n2= n-1;
else
    n2 = n;
end

events = [];
data_temp = zeros(floor(size(data_tcd,2)/1000),1000,nIC);
for ic = 1:nIC
    start = 1;
    for i = 1:floor(size(data_tc,2)/1000)
        data_temp(i,:,ic) = data_tc(ic,start:start+1000-1);
        start = start+1000;
    end
end
            
figure;
thresh = 1000;
for ic = 1:nIC
    diff_x = diff(data_temp(:,:,ic),[],2);
    ind = find(diff_x(1,:)>thresh);
    ind(find(diff(ind)==1)+1)=[];
    ind_low = find(diff_x(1,:)<1500);
    ind_high = find(diff_x(1,:)>1500);
    ind_close = intersect(ind,ind_low);
    ind_above = intersect(ind,ind_high);
    figure;
    subplot(2,2,1)
    for ii = 1:length(ind_close)
        if ind_close(ii)-5>0
            if ind_close(ii)+10<1000
                plot(data_temp(1,ind_close(ii)-5:ind_close(ii)+10,ic),'b')
                hold on
            end
        end
    end
    ylim([10000 25000])
    subplot(2,2,3)
    for ii = 1:length(ind_above)
        if ind_above(ii)-5>0
            if ind_above(ii)+10<1000
                plot(data_temp(1,ind_above(ii)-5:ind_above(ii)+10,ic),'k')
                hold on
            end
        end
    end
    ylim([10000 25000])
    ind_low = find(diff_x(1,:)<2000);
    ind_high = find(diff_x(1,:)>2000);
    ind_close = intersect(ind,ind_low);
    ind_above = intersect(ind,ind_high);
    subplot(2,2,2)
    for ii = 1:length(ind_close)
        if ind_close(ii)-5>0
            if ind_close(ii)+10<1000
                plot(data_temp(1,ind_close(ii)-5:ind_close(ii)+10,ic),'b')
                hold on
            end
        end
    end
    ylim([10000 25000])
    subplot(2,2,4)
    for ii = 1:length(ind_above)
        if ind_above(ii)-5>0
            if ind_above(ii)+10<1000
                plot(data_temp(1,ind_above(ii)-5:ind_above(ii)+10,ic),'k')
                hold on
            end
        end
    end
    ylim([10000 25000])
    suptitle(['Cell #' num2str(ic) '; ' num2str(length(ind)) ' events'])
end
%custom threshold;
thresh = [ones(1,nIC).*1000];
thresh(1) = 1500;

%find all events to get rate
for ic = 1:nIC
    events(ic).diff_x = diff(data_tc(ic,:),[],2);
    events(ic).ind = find(events(ic).diff_x>thresh(ic));
    events(ic).ind(find(diff(events(ic).ind)==1)+1)=[];
    events(ic).rate = length(events(ic).ind)./((double(ifi)./1000)*size(data_tc,2));
end

min_iei = 10;
figure;
for ic = 1:nIC
    double_events = find(diff(events(ic).ind)<min_iei);
    events(ic).ind_good = events(ic).ind;
    events(ic).ind_good([double_events (double_events+1)])=[];
    ind_out = [];
    events(ic).f_chunk = [];
    for ii = 1: length(events(ic).ind_good)
        if and(events(ic).ind_good(ii)-5>0, events(ic).ind_good(ii)+10<size(data_tc,2))
            events(ic).f_chunk = [events(ic).f_chunk; data_tc(ic,events(ic).ind_good(ii)-5:events(ic).ind_good(ii)+10)];
        else
            ind_out = [ind_out ii];
        end
    end
    events(ic).ind_good(ind_out)=[];
    events(ic).f_base = mean(events(ic).f_chunk(:,4:6),2);
    events(ic).df_chunk = bsxfun(@minus, events(ic).f_chunk, events(ic).f_base);
    events(ic).dfoverf_chunk = bsxfun(@rdivide, events(ic).df_chunk, events(ic).f_base);
    events(ic).dfoverf_chunk_avg = mean(events(ic).dfoverf_chunk,1);
    events(ic).dfoverf_chunk_sem = std(events(ic).dfoverf_chunk,1)./(sqrt(size(events(ic).dfoverf_chunk,1)));
    subplot(4,4,ic)
    ts = [-5:10].*double(ifi);
    %plot(ts, events(ic).dfoverf_chunk)
    errorbar(ts,events(ic).dfoverf_chunk_avg,events(ic).dfoverf_chunk_sem)
end

    
pre_buffer = ceil(300/double(ifi));
post_buffer = ceil(1000/double(ifi));
fr_lever = [];
for ipress = 1:length(lever.press)    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.press(ipress)))-pre_buffer:frame_info.counter(round(lever.press(ipress)))+post_buffer];
end 
for irelease = 1:length(lever.release)    %this could be problematic due to unremoved NaNs
    fr_lever = [fr_lever frame_info.counter(round(lever.release(irelease)))-pre_buffer:frame_info.counter(round(lever.release(irelease)))+post_buffer];
end 

fr_lever = unique(fr_lever);
fr_lever(find(fr_lever<1)) = [];
fr_lever(find(fr_lever>size(data_tc,2))) = [];
fr_lever(find(fr_lever>length(ftimes.frame_times))) = [];

data_tc_spont = data_tc; 
data_tc_spont(:,fr_lever) = [];

data_tc_spont_z = bsxfun(@minus, data_tc_spont, (median(data_tc_spont,2)));
tc_diff = diff(data_tc_spont_z')';
thresh = 1000;
min_iei = 4;

for ic = 1:size(data_tc,1)
    event(ic).frames = find(tc_diff(ic,:)>thresh);
    iei = diff(event(ic).frames);
    event(ic).frames(find(iei==1)) = [];
    iei = diff(event(ic).frames);
%     double_events = find(iei<min_iei);
%     event(ic).frames([double_events double_events+1]) = [];
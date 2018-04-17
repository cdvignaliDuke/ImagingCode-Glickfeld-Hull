function event_std = calcEventStd(spike_per_event_f, spike_per_event_l, zero_event)
 
spike_per_event = fliplr(spike_per_event_f);
[peak_x, peak_y] = find(spike_per_event);
peak_xy = [peak_x peak_y];
[~,idx] = unique(peak_xy(:,2));
peak_ind_f = [peak_x(idx,:) peak_y(idx,:)];

spike_per_event = spike_per_event_l;
[peak_x, peak_y] = find(spike_per_event);
peak_xy = [peak_x peak_y];
[~,idx] = unique(peak_xy(:,2));
peak_ind_l = [peak_x(idx,:) peak_y(idx,:)];

c_ind = find(ismember(peak_ind_l(:,2),peak_ind_f(:,2)));

add_ind = find(~ismember(peak_ind_l(:,2),peak_ind_f(:,2)));

peak_ind = [peak_ind_f; peak_ind_l(add_ind,:)];

if ~isempty(c_ind)
    for ii = c_ind'
        temp = peak_ind_f(peak_ind_f(:,2) == peak_ind_l(ii,2),1);
        if (peak_ind_l(ii,1) < temp)
            peak_ind(peak_ind_f(:,2) == peak_ind_l(ii,2),1) =  peak_ind_l(ii,1);
        end
    end
else
    peak_ind = [peak_ind_f;peak_ind_l(:,1) peak_ind_l(:,2)];
end

peak_ind_t = find(zero_event);

if ~isempty(peak_ind_t)
    t_ind = find(ismember(peak_ind(:,2), peak_ind_t));
    if ~isempty(t_ind)
        peak_ind(t_ind,1) = 0;
    else
        peak_ind = [peak_ind;zeros(length(t_ind),1), t_ind];
    end
end
event_std = std(peak_ind(:,1));
end
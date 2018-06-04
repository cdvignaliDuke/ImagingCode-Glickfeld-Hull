function peak_ind = calcEventStd2(spike_per_event_f, spike_per_event_l, zero_event, cue_offset)
%editted version of calcEventStd which finds the std of the nearest spike relative to lever release. Can be used to find the std across trials for each neuron or across neurons for each trial 
%peak_ind should be an N by 2 matrix. N=the number of cells (or trials). The right column values are the cell (or trial) number. The left column
%values are the indeces of the spike closest to the lever release. 

%flip the pre-lever window so the distance from release is consiststent pre vs post release
spike_per_event = fliplr(spike_per_event_f);
%generate indeces of all spike events. Y axis is cell #
[peak_x, peak_y] = find(spike_per_event);
peak_xy = [peak_x peak_y];
[~,idx] = unique(peak_xy(:,2));
%make the pre-lever window frame #s negative 
peak_ind_f = [peak_x(idx,:)*(-1) peak_y(idx,:)];

%collect frame nums for the post-lever window
spike_per_event = spike_per_event_l;
[peak_x, peak_y] = find(spike_per_event);
peak_xy = [peak_x peak_y];
[~,idx] = unique(peak_xy(:,2));
peak_ind_l = [peak_x(idx,:) peak_y(idx,:)];

%find the index elements of post-lever window which are also in the prelever window (cells which have a spike in both windows)
c_ind = find(ismember(peak_ind_l(:,2),peak_ind_f(:,2)));

%find the indeces of the elements in the post-lever window which are not also in the pre-lever window
add_ind = find(~ismember(peak_ind_l(:,2),peak_ind_f(:,2)));
%add the unique elements to the pre-lever window. So peak_ind has all unique cells from each category 
peak_ind = [peak_ind_f; peak_ind_l(add_ind,:)];

%make sure that this can distinguish times that are closer to the release
%despite negative values.
if ~isempty(c_ind)  %if there are cells with spikes in both windows...
    for ii = c_ind'
        temp = peak_ind_f(peak_ind_f(:,2) == peak_ind_l(ii,2),1); %find the indeces of the cells in post-L window with spikes in both windows
        if (peak_ind_l(ii,1) < abs(temp))  %if the spike in the pre-L window is equidistant or close then use the pre-L window spike. Otherwise use the post-L window spike
            peak_ind(peak_ind_f(:,2) == peak_ind_l(ii,2),1) =  peak_ind_l(ii,1);
        end
    end
end

%find the indeces of cells with spikes in the L window
if size(zero_event,1) > 1
    zero_event = zero_event';
end
peak_ind_t = find(zero_event);
%find cells which are unique to the L window
add_t_ind = find(~ismember(peak_ind_t, peak_ind(:,2)));

if ~isempty(peak_ind_t)
    t_ind = find(ismember(peak_ind(:,2), peak_ind_t)); %find cells which have a spike at L and in either the pre or post L window
    if ~isempty(t_ind)
        peak_ind(t_ind,1) = 0;
        if ~isempty(add_t_ind)
            peak_ind = [peak_ind;        zeros(length(add_t_ind),1), peak_ind_t(add_t_ind)'  ];
        end
    else %if there are no cells wit spikes in L window and one other window then just add the L window cells to the list
        peak_ind = [peak_ind;   zeros(length(peak_ind_t),1), peak_ind_t'];
    end
end

%make sure that no cells are counted twice
assert(length(peak_ind(:,2)) == length(unique(peak_ind(:,2))));

%if you are aligning to cue then shift the spike rasters for each trial 
if ~isempty(cue_offset)
    peak_ind(:,1) = peak_ind(:,1) + double(cue_offset([peak_ind(:,2)]))';
end

%find the standard deviation and return
%event_std = std(peak_ind(:,1));
end








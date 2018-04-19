[x, y, z] = size(downsampled_movie);
events = reshape(downsampled_movie, x*y, z);
events_diff_x = bsxfun(@diff, events, events);
nIC = x*y;
normalization = 1;
deconvtau = 0;
dt = 1;
for ic = 1:nIC
    %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
    [~, events_ind{ic}, ~] = CellsortFindspikes(events_diff_x(:,ic), thresh, dt, deconvtau, normalization);
    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];
    events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
end
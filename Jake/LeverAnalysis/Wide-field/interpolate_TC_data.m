function inerpolated_TCs = interpolate_TC_data(initial_TCs, step_size);
%function for interpolating df/f timecourses: used in peak_f_vs_bx_events
%and struct_maker
XVector = [1:.01:size(initial_TCs,3)];
inerpolated_TCs = nan(length(XVector), size(initial_TCs,1), size(initial_TCs,2)); %dims: 1=T2 2=trial# 3=ROI#
for ROI_num = 1:size(initial_TCs,2); %for each roi...
    TC_temp = squeeze(initial_TCs(:,ROI_num,:))';  %dims 1=Time 2=trial number
    TC_temp = interp1(TC_temp, XVector);
    inerpolated_TCs(:,:,ROI_num) = TC_temp;  %dims 1=Time 2=trialNumber  3=roi
end

return
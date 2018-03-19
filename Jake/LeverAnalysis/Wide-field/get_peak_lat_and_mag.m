function [TC_tbyt_lat_peak, TC_tbyt_peak_mag] = get_peak_lat_and_mag(TC_roi, analysis_window, release_frame, peak_mag_window)
   
%Determine time of the peak df/f on a trial by trial basis for each ROI
%TC dims must be as follows: dim1=trial_num dim2=ROI_num dim3=frame_num
%or dim1= trial_num dim2=frame_num

%if two dims then convert to three
if length(size(TC_roi))<3
    TC_roi = reshape(TC_roi, [size(TC_roi,1), 1, size(TC_roi,2)]);
end

%for each ROI go through each trial and find the max value within the window
TC_tbyt_lat_peak = [];
TC_tbyt_peak_mag = [];
for ROI_num = 1:size(TC_roi,2)
    TC_tbyt_lat_peak_roi = [];
    TC_tbyt_peak_mag_roi = [];
    for trial_num = 1:size(TC_roi,1);
        temp_val = find(TC_roi(trial_num,ROI_num,:)==max(TC_roi(trial_num, ROI_num, [analysis_window])));  %finds the peak time of each trials TC
        if size(temp_val,1) > 1 %in case of finding multiple values it will select the first one
            temp_val = temp_val(find(temp_val >= analysis_window(1),1, 'first'));
        end
        TC_tbyt_lat_peak_roi = [TC_tbyt_lat_peak_roi, temp_val]; %give peak time for each trial. Will occur in increments of 10 because the peak will always be on one of the actual frames. Not interpolated data points
        TC_tbyt_peak_mag_roi = [TC_tbyt_peak_mag_roi, mean(TC_roi(trial_num, ROI_num, [temp_val-peak_mag_window:temp_val+peak_mag_window]))]; %collects all the trial by trial peak magnitudes
    end
    if ROI_num == 1
        TC_tbyt_lat_peak = TC_tbyt_lat_peak_roi;
        TC_tbyt_peak_mag = TC_tbyt_peak_mag_roi;
    elseif ROI_num > 1
        TC_tbyt_lat_peak = [TC_tbyt_lat_peak; TC_tbyt_lat_peak_roi]; %dim1=ROI    dim2=trial_num
        TC_tbyt_peak_mag = [TC_tbyt_peak_mag; TC_tbyt_peak_mag_roi];
    end
end
TC_tbyt_lat_peak = TC_tbyt_lat_peak-release_frame;

return

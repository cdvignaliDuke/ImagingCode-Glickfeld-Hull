function licking_data = bin_licking_by_frame(lickTimes, frame_info);
%bins # of licks by frame times. Used in WF_lever_plotting_TCs
frameTimes = frame_info.times; %frame times is monotonic by the time it gets here.
licksByFrame = [];
for i = 1:length(frameTimes);
    if i == 1; %if this is the first frameTime
        licksThisFrame = sum(lickTimes>=1 & lickTimes<=frameTimes(i));
        licksByFrame = [licksByFrame licksThisFrame];
    else
        licksThisFrame = sum(lickTimes>frameTimes(i-1) & lickTimes<=frameTimes(i));
        licksByFrame = [licksByFrame licksThisFrame];
    end
end
licking_data = [];
licking_data.lickTimes = lickTimes;
licking_data.licksByFrame = licksByFrame;
return
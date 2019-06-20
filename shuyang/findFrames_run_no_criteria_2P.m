% NOTE:similar as findFrame_runWindows_2P, but this one includes all running
%      windows regardless if the animal is stationary all the time during that
%      1s before/after running 
%generate cell for all windows right before and after running (300ms)
% matrix for run triggered average (500ms before run and 1s from running onset)
% matrix for running windows w/ 300ms before and after
function[frames_befo_run_cell,frames_aft_run_cell,frames_runTrigger_NC_mat,frms_runoff_NC_mat,...
    frames_run_mat]= findFrames_run_no_criteria_2P (speed,frames_run_cell,period,befoRun,totalT,aftRunOff)

frames = 1: length(speed);

%% generate frames_befo/aft run, and runTriggers

%period = 9; % # of frames before and after running
%befoRun = 30; %1s, # of frames that speed =0 before running
%totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
%aftRunOff = 30; %1s,# of frames that speed = 0 after running
frames_befo_run_cell = {};
%was trying to generate a matrix in the for loop but then when m=1, if it doesn't fullfill the requirement and just continues, 
%the first line is going to be zeros. and matlab doesn't do (end,:) if the variable is initialized to []. so cell is easier
frames_aft_run_cell = {}; 
frames_runTrigger_NC_mat = []; % frames_runTrigger is for run triggered average analysis. 
frms_runoff_NC_mat = [];
for m = 1: size(frames_run_cell,2)
        if (frames_run_cell{m}(1)-period < 1) || (frames_run_cell{end}(end)+ period > frames(end))
            continue
        elseif (frames_run_cell{m}(1)-befoRun <1) || (frames_run_cell{end}(end)+ aftRunOff > frames(end))
            continue
        else
            frames_befo_run_cell = cat(2, frames_befo_run_cell, frames_run_cell{m}(1)-period:frames_run_cell{m}(1)-1);
            frames_aft_run_cell =  cat(2, frames_aft_run_cell, frames_run_cell{m}(end)+1:frames_run_cell{m}(end)+ period);
            frames_runTrigger_temp = frames_run_cell{m}(1)-befoRun : frames_run_cell{m}(end);
            if length(frames_runTrigger_temp) >= totalT %&& sum(speed(frames_runTrigger_temp(1:befoRun)) == 0) == befoRun 
               frames_runTrigger_NC_mat = cat(2, frames_runTrigger_NC_mat, frames_runTrigger_temp(1:totalT));
            end
            frms_runoff_temp = frames_run_cell{m}(1):frames_run_cell{m}(end)+aftRunOff;
            if length(frms_runoff_temp) >= totalT %&& sum(speed(frms_runoff_temp(end-aftRunOff+1:end))==0)== aftRunOff
                frms_runoff_NC_mat = cat(2,frms_runoff_NC_mat, frms_runoff_temp(end-totalT+1:end));
            end
         
            % frames_runTrigger is the vectors contain 500ms still before run and the first 1s of each running window
        end
end

%% create matrix for frames_runTrigger, this can be used for triggered_averaging plot
frames_runTrigger_NC_mat = reshape(frames_runTrigger_NC_mat, totalT, length(frames_runTrigger_NC_mat)/totalT);
frms_runoff_NC_mat = reshape(frms_runoff_NC_mat, totalT, length(frms_runoff_NC_mat)/totalT);

%% create matrix for frames of running windows with buffer of 3 frames before and after all windows, this will be useful when drawing the figures. 
%each line in frames_run_buffer_mat is a single running window.
run_length_temp = cell2mat(cellfun(@size,frames_run_cell, 'UniformOutput',0));
run_length = max(run_length_temp);
frames_run_mat = nan(size(frames_run_cell,2), run_length);
for n = 1: size(frames_run_mat,1)
    temp = frames_run_cell{n};
    frames_run_mat(n,1:size(temp,2)) = temp;
end

end






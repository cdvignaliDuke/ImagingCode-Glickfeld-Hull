%used in twoP_MDspikes
%generate cell for all windows right before and after running (300ms)
% matrix for run triggered average (500ms before run and 1s from running onset)
% matrix for running windows w/ 300ms before and after
function[frames_befo_run_cell,frames_aft_run_cell,frames_runTrigger_mat,frms_runoff_mat, frames_run_mat]= findFrames_runWindows_2P (speed,frames_run_cell)

frames = 1: length(speed);

%% generate frames_befo/aft run, and runTriggers

period = 9; % 300ms % numbers in this script should change b/c the sampling rate is now 3 times of wide field
befoRunStay = 15; %500ms
runTriggerDura = 45; %1.5s
aftRunOff = 15; %500ms
frames_befo_run_cell = {};
%was trying to generate a matrix in the for loop but then when m=1, if it doesn't fullfill the requirement and just continues, 
%the first line is going to be zeros. and matlab doesn't do (end,:) if the variable is initialized to []. so cell is easier
frames_aft_run_cell = {}; 
frames_runTrigger_mat = []; % frames_runTrigger is for run triggered average analysis. 
frms_runoff_mat = [];
for m = 1: size(frames_run_cell,2)
        if (frames_run_cell{m}(1)-period < 1) || (frames_run_cell{end}(end)+ period > frames(end))
            continue
        elseif (frames_run_cell{m}(1)-befoRunStay <1) || (frames_run_cell{end}(end)+ aftRunOff > frames(end))
            continue
        else
            frames_befo_run_cell = cat(2, frames_befo_run_cell, frames_run_cell{m}(1)-period:frames_run_cell{m}(1)-1);
            frames_aft_run_cell =  cat(2, frames_aft_run_cell, frames_run_cell{m}(end)+1:frames_run_cell{m}(end)+ period);
            frames_runTrigger_temp = frames_run_cell{m}(1)-befoRunStay : frames_run_cell{m}(end);
            if length(frames_runTrigger_temp) >= runTriggerDura && sum(speed(frames_runTrigger_temp(1:befoRunStay)) == 0) == befoRunStay 
               frames_runTrigger_mat = cat(2, frames_runTrigger_mat, frames_runTrigger_temp(1:runTriggerDura));
            end
            frms_runoff_temp = frames_run_cell{m}(1):frames_run_cell{m}(end)+aftRunOff;
            if length(frms_runoff_temp) >= runTriggerDura && sum(speed(frms_runoff_temp(end-aftRunOff+1:end))==0)== aftRunOff
                frms_runoff_mat = cat(2,frms_runoff_mat, frms_runoff_temp(end-runTriggerDura+1:end));
            end
         
            % frames_runTrigger is the vectors contain 500ms still before run and the first 1s of each running window
        end
end

%% create matrix for frames_runTrigger, this can be used for triggered_averaging plot
runTriggerDura = 45;
frames_runTrigger_mat = reshape(frames_runTrigger_mat, runTriggerDura, length(frames_runTrigger_mat)/runTriggerDura);
frms_runoff_mat = reshape(frms_runoff_mat, runTriggerDura, length(frms_runoff_mat)/runTriggerDura);

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






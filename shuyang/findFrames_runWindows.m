function[frames_befo_run,frames_aft_run,frames_runTrigger,frames_run_buffer_mat]= findFrames_runWindows (speed,frames_run_cell)

frames = 1: length(speed);

%% generate frames_befo/aft run,frames_run_buffer(run windows+a buffer of 3 frames before and after each window), and runTriggers

frames_befo_run = {};
frames_aft_run = {};
frames_run_buffer ={}; %frames_run_buffer is for analysis of df/f right before and right after every running window.
frames_runTrigger = []; % frames_runTrigger is for run triggered average analysis.
period = 3;
befoRunStay = 5;
runTriggerDura = 15;
for m = 1: size(frames_run_cell,2)
        if (frames_run_cell{m}(1)-period < 1) || (frames_run_cell{end}(end)+ period > frames(end))
            continue
        elseif (frames_run_cell{m}(1)-befoRunStay <1) || (frames_run_cell{end}(end)+ befoRunStay > frames(end))
            continue
        else
            frames_befo_run = cat(2, frames_befo_run,frames_run_cell{m}(1)-period:frames_run_cell{m}(1)-1);
            frames_aft_run = cat(2, frames_aft_run, frames_run_cell{m}(end)+1:frames_run_cell{m}(end)+ period);
            frames_run_buffer = cat(2, frames_run_buffer, frames_run_cell{m}(1)-period:frames_run_cell{m}(end)+ period);
            frames_runTrigger_temp = frames_run_cell{m}(1)-befoRunStay : frames_run_cell{m}(end);
            if length(frames_runTrigger_temp) >= runTriggerDura && sum(speed(frames_runTrigger_temp(1:5)) == 0) == 5
               frames_runTrigger = cat(2, frames_runTrigger, frames_runTrigger_temp(1:15));
            end
            % frames_runTrigger is the vectors contain 500ms still
            % before run and the first 1s of each running window
        end
end

%% create matrix for frames of running windows with buffer of 3 frames before and after all windows, this will be useful when drawing the figures. 
%each line in frames_run_buffer_mat is a single running window.
run_length_temp = cell2mat(cellfun(@size,frames_run_buffer, 'UniformOutput',0));
run_length = max(run_length_temp);
frames_run_buffer_mat = nan(size(frames_run_buffer,2), run_length);
for n = 1: size(frames_run_buffer_mat,1)
    temp = frames_run_buffer{n};
    frames_run_buffer_mat(n,1:size(temp,2)) = temp;
end

end






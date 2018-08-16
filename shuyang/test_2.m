frames_bf_cell = {};
frames_run_cell = {};
frames_befo_run = {};
frames_aft_run = {};
frames_run_buffer ={}; %frames_run_buffer is for analysis of df/f right before and right after every running window.
frames_runTrigger = []; % frames_runTrigger is for run triggered average analysis.
period = 3;
befoRunStay = 5;
runTriggerDura = 15;
for m = 1: size(frames_move_cell,2)
    bf = find(max(speed(frames_move_cell{m})) <= 10);
    if ~isempty(bf) &&  min(speed(frames_move_cell{m})) < 0
        frames_bf_cell = cat(2, frames_bf_cell, frames_move_cell{m});
    else
        frames_run_cell = cat(2, frames_run_cell, frames_move_cell{m});
        if (frames_move_cell{m}(1)-period < 1) || (frames_move_cell{end}(end)+ period > frames(end))
            continue
        elseif (frames_move_cell{m}(1)-befoRunStay <1) || (frames_move_cell{end}(end)+ befoRunStay > frames(end))
            continue
        else
            frames_befo_run = cat(2, frames_befo_run,frames_move_cell{m}(1)-period:frames_move_cell{m}(1)-1);
            frames_aft_run = cat(2, frames_aft_run, frames_move_cell{m}(end)+1:frames_move_cell{m}(end)+ period);
            frames_run_buffer = cat(2, frames_run_buffer, frames_move_cell{m}(1)-period:frames_move_cell{m}(end)+ period);
            frames_runTrigger_temp = frames_move_cell{m}(1)-befoRunStay : frames_move_cell{m}(end);
            if length(frames_runTrigger_temp) >= runTriggerDura && sum(speed(frames_runTrigger_temp(1:5)) == 0) == 5
               frames_runTrigger = cat(2, frames_runTrigger, frames_runTrigger_temp(1:15));
            end
            % frames_runTrigger is the vectors contain 500ms still
            % before run and the first 1s of each running window
        end
    end
end
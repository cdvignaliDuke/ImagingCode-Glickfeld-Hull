dfOvF_runtrig_mat = zeros(size(frames_runTrigger_mat,1),size(frames_runTrigger_mat,2),size(TCave,2));
for i = 1: size(frames_runTrigger_mat,2)                                    % for each running trig window
    dfOvF_runtrig_mat(:,i,:) = dfOvF(frames_runTrigger_mat(:,i),:);                 % df/f for all cells during that window, frame*cells
    
end
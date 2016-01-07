%run initial analysis before combining runs
clear all
Jake_2P_exptlist
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    nrun = nrun_mat(id,:);
    for irun = 1:nrun
        run = run_mat(id,:,irun);
        disp([date ' ' mouse])
        run_name = [date '_' mouse '_run' run(length(run)-2:end)];
        out_path = fullfile(out_base,run_name);
        dest =  fullfile(out_path);   %altered 11/8/15 removed "run_name" duplicate
        dest_sub = fullfile([dest '_nosub']);  %altered 11/9/15 jh
        prepare_movie_for_analysis_2P   %added by JH 11/9/15
        HAD_cmp_success_fail_2P_ROIs
        HAD_2P_event_detection 
    end
end

%run summary analysis after combining runs
clear all
Jake_2P_exptlist
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    nrun = nrun_mat(id,:);
    run = run_mat(id,:,1);
    disp([date ' ' mouse])
    if nrun == 1
        run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    else
        run_name = [date '_' mouse '_run' run(length(run)-2:end) '-00' num2str(nrun-1)];
    end
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);
    dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
    %HAD_2P_TC_quantification  
    HAD_2P_event_quantification
    close all
end
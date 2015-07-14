clear all
base_dir =  '\\CRASH.dhe.duke.edu\data\home\jake\Data\2P_imaging\';
date_mat = strvcat('150703', '150703','150704','150704');
run_mat = strvcat('_000_000', '_000_001', '_000_001', '_000_000');
mouse_mat = strvcat('img24','img25','img24','img25');
subNum_mat = strvcat('924','925','924','925');
time_mat = strvcat('1821', '2005','1845','1937');
out_base = 'Z:\home\lindsey\Analysis\2P\Jake\';
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    run = run_mat(id,:);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);

    run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);
    dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
    load([dest '_ROI_TCs.mat']);
%     load([dest '_npSubTCs.mat']);
%     data_tc = npSubTC;
    
    
    load([dest_sub '_press_resp_by_hold.mat']);
    load([dest_sub '_release_resp_by_outcome.mat']);
    load([dest '_ICs.mat']);
    %prepare_movie_for_analysis_2P
    %HAD_cmp_success_fail_2P_ROIs
    HAD_2P_quantification
end
suptitle([date ' ' mouse ' Avg resp:' tc_type ' Success (black), Fail (red), Press (cyan)']);
print([out_base 'Summary_allcell_responses_' tc_type '.eps'], '-depsc');
print([out_base 'Summary_allcell_responses_' tc_type '.pdf'], '-dpdf');

clear all
date_mat = strvcat('150703', '150703','150704','150704');
run_mat = strvcat('_000_000', '_000_001', '_000_001', '_000_000');
mouse_mat = strvcat('img24','img25','img24','img25');
out_base = 'Z:\home\lindsey\Analysis\2P\Jake';
figure;
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    run = run_mat(id,:);
    mouse = mouse_mat(id,:);
    
    run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);
    load([dest '_press_resp_by_hold.mat']);
    load([dest '_release_resp_by_outcome.mat']);
    
    nCells = size(success_movie,2);
    success_base = squeeze(mean(success_movie(:,:,1:pre_release_frames),3));
    success_resp = squeeze(mean(success_movie(:,:,pre_release_frames+1:pre_release_frames+2),3));
    [success_h, success_p] = ttest(success_base, success_resp, 'dim', 1, 'tail', 'both');
    
    fail_base = squeeze(mean(fail_movie(:,:,1:pre_release_frames),3));
    fail_resp = squeeze(mean(fail_movie(:,:,pre_release_frames+1:pre_release_frames+2),3));
    [fail_h, fail_p] = ttest(fail_base, fail_resp, 'dim', 1, 'tail', 'both');
    
    success_resp_cells = find(success_h);
    success_resp_avg = mean((success_resp-success_base),1);
    success_resp_sem = std((success_resp-success_base),[],1)./sqrt(size(success_resp,1));

    fail_resp_cells = find(fail_h);
    fail_resp_avg = mean((fail_resp-fail_base),1);
    fail_resp_sem = std((fail_resp-fail_base),[],1)./sqrt(size(fail_resp,1));
    
    fail_only_cells = fail_resp_cells(find(ismember(fail_resp_cells, success_resp_cells)==0));
    success_only_cells = success_resp_cells(find(ismember(success_resp_cells, fail_resp_cells)==0));
    fail_and_success_cells = find(and(fail_h, success_h));
    
    [fail_v_success_h, fail_v_success_p] = ttest(fail_resp_avg(1,fail_and_success_cells),success_resp_avg(1,fail_and_success_cells), 'tail', 'left');
    
    press_base = squeeze(mean(press_long_movie(:,:,1:pre_press_frames-5),3));
    press_resp = squeeze(mean(press_long_movie(:,:,pre_press_frames-1:pre_press_frames),3));
    [press_h, press_p] = ttest(press_base, press_resp, 'dim', 1, 'tail', 'both');
    
    press_resp_cells = find(press_h);
    press_resp_avg = mean((press_resp-press_base),1);
    press_resp_sem = std((press_resp-press_base),[],1)./sqrt(size(press_resp,1));
    
    press_and_success_cells = find(and(success_h,press_h));
    press_notsuccess_cells = press_resp_cells(find(ismember(press_resp_cells, success_resp_cells)==0));
    success_notpress_cells = success_resp_cells(find(ismember(success_resp_cells, press_resp_cells)==0));
    
    noresponse_cells = find(sum(success_h,fail_h,press_h)==0);
end
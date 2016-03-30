clear all
Jake_2P_exptlist
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    run = run_mat(id,:,1);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    nrun = nrun_mat(id,:);
    if nrun == 1
        run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    else
        run_name = [date '_' mouse '_run' run(length(run)-2:end) '-00' num2str(nrun-1)];
    end
    
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);
    
    dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
    load([dest_sub '_cell_resp.mat']);
    load([dest_sub '_cell_categories.mat']);
    
    ncells(id) = size(press_resp,2);
    tot_resp(id) = length(unique([success_resp_cells fail_resp_cells press_resp_cells]));
    release_resp = unique([success_resp_cells fail_resp_cells]);
    tot_release(id) = length(release_resp);
    pct_resp(id) = tot_resp(id)./ncells(id);
    pct_release(id) = tot_release(id)./tot_resp(id);
    pct_success(id) = length(success_resp_cells)./tot_resp(id);
    pct_fail(id) = length(fail_resp_cells)./tot_resp(id);
    pct_press(id) = length(press_resp_cells)./tot_resp(id);
    overlap_release(id) = length(intersect(success_resp_cells,fail_resp_cells))./tot_release(id);
    overlap_press(id) = length(intersect(press_resp_cells,release_resp))./tot_resp(id);
    success_only(id) = length(success_only_cells)./tot_resp(id);
    fail_only(id) = length(fail_only_cells)./tot_resp(id);
    press_only(id) = length(press_resp_cells(find(ismember(press_resp_cells,release_resp)==0)))./tot_resp(id);
    release_only(id) = length(release_resp(find(ismember(release_resp,press_resp_cells)==0)))./tot_resp(id);
end
figure;
col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c');
%number of cells
subplot(3,3,1), scatter(1:size(date_mat,1), ncells, 'ok')
xlabel('Expt #')
ylabel('#')
ylim([0 25])
xlim([0 size(date_mat,1)+1])
title('# of cells')
%number driven cells
subplot(3,3,2), scatter(1:size(date_mat,1), tot_resp, 'ok')
xlabel('Expt #')
ylabel('#')
ylim([0 25])
xlim([0 size(date_mat,1)+1])
title('# of responsive cells')
%number driven cells
subplot(3,3,3), scatter(1:size(date_mat,1), pct_resp, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% of responsive cells')

%percent cells driven by release/press
subplot(3,3,4), scatter(1:size(date_mat,1), pct_release, 'ok'); hold on; scatter(1:size(date_mat,1), pct_press, 'oc')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by release/press')
%percent driven by success/fail
subplot(3,3,5), scatter(1:size(date_mat,1), pct_success, 'ok'); hold on; scatter(1:size(date_mat,1), pct_fail, 'or');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by success/fail')

%percent overlap of success/fail
subplot(3,3,6), scatter(1:size(date_mat,1), overlap_release, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by success & fail')
%percent overlap of press/release
subplot(3,3,7), scatter(1:size(date_mat,1), overlap_press, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by press & release')
%percent success only and fail only
subplot(3,3,8), scatter(1:size(date_mat,1), success_only, 'ok'); hold on; scatter(1:size(date_mat,1), fail_only, 'or');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by success/fail only')
%percent release only and press only
subplot(3,3,9), scatter(1:size(date_mat,1), release_only, 'ok'); hold on; scatter(1:size(date_mat,1), press_only, 'oc');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(date_mat,1)+1])
title('% driven by release/press only')


suptitle(['Summary of % cell response stats'])
    print([out_base 'Summary_pct_cell_response_stats.eps'], '-depsc');
    print([out_base 'Summary_pct_cell_response_stats.pdf'], '-dpdf');
    

figure;
subplot(2,2,1)
for id = 1:size(date_mat,1)
    scatter(pct_release(id), pct_press(id), col_mat(id,:))
    hold on
end
x = 0:.1:1;
y = x;
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release responsive (%)')
ylabel('Press responsive (%)')
subplot(2,2,2)
for id = 1:size(date_mat,1)
    scatter(release_only(id), press_only(id), col_mat(id,:))
    hold on
end
x = 0:.1:1;
y = x;
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Release only (%)')
ylabel('Press only (%)')
subplot(2,2,3)
for id = 1:size(date_mat,1)
    scatter(pct_success(id), pct_fail(id), col_mat(id,:))
    hold on
end
x = 0:.1:1;
y = x;
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success responsive (%)')
ylabel('Fail responsive (%)')
subplot(2,2,4)
for id = 1:size(date_mat,1)
    scatter(success_only(id), fail_only(id), col_mat(id,:))
    hold on
end
x = 0:.1:1;
y = x;
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success only (%)')
ylabel('Fail only (%)')
suptitle('% of all responsive cells- Red: img24; Blue: img25; Magenta: img27 Purple: img28 Cyan: img32;')
    print([out_base 'Summary_pct_cell_response_scatter.eps'], '-depsc');
    print([out_base 'Summary_pct_cell_response_scatter.pdf'], '-dpdf');


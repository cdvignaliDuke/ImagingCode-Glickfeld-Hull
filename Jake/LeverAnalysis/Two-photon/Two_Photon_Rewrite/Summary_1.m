clear 
file_info
% out_base = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis\');
out_base = fullfile('Z:\home\ziye\2P_Analysis\2P_Analysis\');
mouseID = mouseID(1:30);
tot_cells = 0; tot_sus = 0; tot_trans = 0; press_trans_tot = 0; press_sus_tot = 0;
tot_success = 0; tot_fail = 0; tot_press = 0;
for id = 1:size(mouseID,2)
    for rID  = 1:2
%         dest_sub  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        dest_sub = ['Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure\', date{id}, '_', runID{rID}, '_', mouseID{id}, '\'];
        if exist(dest_sub)
            load([dest_sub '_cell_resp.mat']);
            load([dest_sub '_cell_categories.mat']);
            load([dest_sub 'ROI_TCs.mat']);
            load([out_base, 'cell_count.mat']);
%             load([dest_sub 'Reg_out.mat']);

            tot_cells = size(tc_avg,2) + tot_cells;
            tot_sus = size(sustain_cell_ind,1) + tot_sus;
            tot_trans = size(trans_cell_ind,1) + tot_trans;
            press_trans = intersect(press_resp_cells,trans_cell_ind');
            press_sus = intersect(press_resp_cells,sustain_cell_ind');
            
            press_trans_tot = size(press_trans,2) + press_trans_tot;
            press_sus_tot = size(press_sus,2) + press_sus_tot;
            
            ncells(id) = size(press_resp,2);
            tot_resp(id) = length(unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]));
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
            
            tot_success = tot_success + length(success_resp_cells);
            tot_fail = tot_fail + length(fail_resp_cells);
            tot_press = tot_press + length(press_resp_cells);
        end
    end
end
tot_resp_cells = sum(tot_resp);
pct_success_tot = tot_success/tot_resp_cells*100;
pct_fail_tot = tot_fail/tot_resp_cells*100;
pct_press_tot = tot_press/tot_resp_cells*100;
% col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c');
col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0; 0.5 0.5 0.5; 0.3 0.5 1; 0.1 0.5 0.7; 0 0.6 0.2;0.8 0.8 0.4;0.1 0.1 0.1;0.3 0.7 0; 0 0 0; 0.1 0.5 0.5];
% col_mat = repmat([0 0 0], id, 1);

%number of cells
fig = figure;
subplot(3,3,1), scatter(1:size(mouseID,2), ncells, 'ok')
xlabel('Expt #')
ylabel('#')
ylim([0 25])
xlim([0 size(mouseID,2)+1])
title('# of cells')
%number driven cells
subplot(3,3,2), scatter(1:size(mouseID,2), tot_resp, 'ok')
xlabel('Expt #')
ylabel('#')
ylim([0 25])
xlim([0 size(mouseID,2)+1])
title('# of responsive cells')
%number driven cells
subplot(3,3,3), scatter(1:size(mouseID,2), pct_resp, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% of responsive cells')

%percent cells driven by release/press
subplot(3,3,4), scatter(1:size(mouseID,2), pct_release, 'ok'); hold on; scatter(1:size(mouseID,2), pct_press, 'oc')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by release/press')
%percent driven by success/fail
subplot(3,3,5), scatter(1:size(mouseID,2), pct_success, 'ok'); hold on; scatter(1:size(mouseID,2), pct_fail, 'or');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by success/fail')

%percent overlap of success/fail
subplot(3,3,6), scatter(1:size(mouseID,2), overlap_release, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by success & fail')
%percent overlap of press/release
subplot(3,3,7), scatter(1:size(mouseID,2), overlap_press, 'ok')
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by press & release')
%percent success only and fail only
subplot(3,3,8), scatter(1:size(mouseID,2), success_only, 'ok'); hold on; scatter(1:size(mouseID,2), fail_only, 'or');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by success/fail only')
%percent release only and press only
subplot(3,3,9), scatter(1:size(mouseID,2), release_only, 'ok'); hold on; scatter(1:size(mouseID,2), press_only, 'oc');
xlabel('Expt #')
ylabel('%')
ylim([0 1])
xlim([0 size(mouseID,2)+1])
title('% driven by release/press only')


supertitle(['Summary of % cell response stats'])
saveas(fig, [out_base 'Summary_pct_cell_response_stats.fig']);
print([out_base 'Summary_pct_cell_response_stats.eps'], '-depsc');
print([out_base 'Summary_pct_cell_response_stats.pdf'], '-dpdf');


fig=figure;
subplot(2,2,1)
cid = 1;
for id = 1:size(mouseID,2)
    id
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(pct_release(id), pct_press(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid + 1; 
            scatter(pct_release(id), pct_press(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(pct_release(id), pct_press(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
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
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(release_only(id), press_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid +1;
            scatter(release_only(id), press_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(release_only(id), press_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
    
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
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(pct_success(id), pct_fail(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid +1;
            scatter(pct_success(id), pct_fail(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
           
        end
    else
        scatter(pct_success(id), pct_fail(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
   
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
cid = 1;
for id = 1:size(mouseID,2)
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(success_only(id), fail_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid +1;
            scatter(success_only(id), fail_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(success_only(id), fail_only(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
   
    hold on
end
x = 0:.1:1;
y = x;
plot(x,y,'-k')
xlim([0 1])
ylim([0 1])
xlabel('Success only (%)')
ylabel('Fail only (%)')
supertitle('% of all responsive cells')
legend(mouseID)
saveas(fig, [out_base 'Summary_pct_cell_response_scatter.fig']);
print([out_base 'Summary_pct_cell_response_scatter.eps'], '-depsc');
print([out_base 'Summary_pct_cell_response_scatter.pdf'], '-dpdf');


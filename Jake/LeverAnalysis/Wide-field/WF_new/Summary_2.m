clear
file_info_WF
out_base = 'Z:\public\Motor Timing Paper\Ziye\WF_PCA_ICA\';

for id = 1:size(mouseID,2)
        data_dir = fullfile('Z:\home\jake\Analysis\WF Lever Analysis\PCA_ICA_output\');
        dest_sub = fullfile(data_dir,[date{id}, '\']);
        
%         dest_sub  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        if exist(dest_sub)
            load([dest_sub '_cell_TCs.mat']);
            load([dest_sub '_cell_resp.mat']);
%             load([dest_sub '_cell_resp_amp.mat']);
            load([dest_sub '_cell_categories.mat']);
            load([dest_sub '_release_movies.mat'])
            
            ncells{id} = size(press_resp,2);
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells tooFast_resp_cells]);
            tot_resp{id} = length(RS_cells{id});
            
            TC_length{id} = size(avg_release,2);
            TC_ifi{id} = ifi;
            pre_frames{id} = pre_release_frames;
            post_frames{id} = post_release_frames;
            
            RL_cells{id} = release_resp_cells;
            tot_RL_resp{id} = length(RL_cells{id});
            RL_trans_cells{id} = trans_cell_ind;
            RL_sustain_cells{id} = sustain_cell_ind;
            RL_trans_cellNum(id) = length(trans_cell_ind);
            RL_sustain_cellNum(id) = length(sustain_cell_ind);
            
            press_trans_cells(id) = length(intersect(press_resp_cells, trans_cell_ind'));
            press_sus_cells(id) = length(intersect(press_resp_cells, sustain_cell_ind'));
            
            release_resp_all{id} = mean((release_resp-release_base),1);
            release_resp_mean{id} = mean(release_resp_all{id},2);
            release_resp_sem{id} = std(release_resp_all{id},[],2)./sqrt(ncells{id});
            release_resp_RL{id} = mean((release_resp(:,release_resp_cells)-release_base(:,release_resp_cells)),1);
            release_resp_RL_mean{id} = mean(release_resp_RL{id},2);
            release_resp_RL_sem{id} = std(release_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            release_resp_RS{id} = mean((release_resp(:,RS_cells{id})-release_base(:,RS_cells{id})),1);
            release_resp_RS_mean{id} = mean(release_resp_RS{id},2);
            release_resp_RS_sem{id} = std(release_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
%             release_resp_peak_mean{id} = mean(mean(release_peak));
            
            press_TC{id} = avg_press;
            press_TC_RS{id} = avg_press(RS_cells{id},:);
            press_TC_RL{id} = avg_press(RL_cells{id},:);
            press_TC_mean{id} = mean(avg_press,1);
            press_TC_sem{id} = std(avg_press,1)./sqrt(size(avg_press,1));
            press_TC_RS_mean{id} = mean(avg_press(RS_cells{id},:),1);
            press_TC_RS_sem{id} = std(avg_press(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            press_TC_RL_mean{id} = mean(avg_press(release_resp_cells,:),1);
            press_TC_RL_sem{id} = std(avg_press(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            press_resp_all{id} = mean((press_resp-press_base),1);
            press_resp_mean{id} = mean(press_resp_all{id},2);
            press_resp_sem{id} = std(press_resp_all{id},[],2)./sqrt(ncells{id});
            press_resp_RL{id} = mean((press_resp(:,release_resp_cells)-press_base(:,release_resp_cells)),1);
            press_resp_RL_mean{id} = mean(press_resp_RL{id},2);
            press_resp_RL_sem{id} = std(press_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            press_resp_RS{id} = mean((press_resp(:,RS_cells{id})-press_base(:,RS_cells{id})),1);
            press_resp_RS_mean{id} = mean(press_resp_RS{id},2);
            press_resp_RS_sem{id} = std(press_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
%             press_resp_peak_mean{id} = mean(mean(press_peak));
            
            success_TC{id} = avg_success;
            success_TC_RS{id} = avg_success(RS_cells{id},:);
            success_TC_RL{id} = avg_success(RL_cells{id},:);
            success_TC_trans{id} = avg_success(trans_cell_ind,:);
            success_TC_sus{id} = avg_success(sustain_cell_ind,:);
            success_TC_mean{id} = mean(avg_success,1);
            success_TC_sem{id} = std(avg_success,1)./sqrt(size(avg_success,1));
            success_TC_RS_mean{id} = mean(avg_success(RS_cells{id},:),1);
            success_TC_RS_sem{id} = std(avg_success(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            success_TC_RL_mean{id} = mean(avg_success(release_resp_cells,:),1);
            success_TC_RL_sem{id} = std(avg_success(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            success_resp_all{id} = mean((success_resp-success_base),1);
            success_resp_mean{id} = mean(success_resp_all{id},2);
            success_resp_sem{id} = std(success_resp_all{id},[],2)./sqrt(ncells{id});
            success_resp_RS{id} = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
            success_resp_RS_mean{id} = mean(success_resp_RS{id},2);
            success_resp_RS_sem{id} = std(success_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            success_resp_RL{id} = mean((success_resp(:,release_resp_cells)-success_base(:,release_resp_cells)),1);
            success_resp_RL_mean{id} = mean(success_resp_RL{id},2);
            success_resp_RL_sem{id} = std(success_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            success_resp_SC{id} = mean((success_resp(:,success_resp_cells)-success_base(:,success_resp_cells)),1);
            success_resp_SC_mean{id} = mean(success_resp_SC{id},2);
            success_resp_SC_sem{id} = std(success_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            success_resp_RL_trans{id} = mean((success_resp(:,trans_cell_ind)-success_base(:,trans_cell_ind)),1);
            success_resp_RL_sus{id} = mean((success_resp(:,sustain_cell_ind)-success_base(:,sustain_cell_ind)),1);
%             success_resp_peak_mean{id} = mean(mean(success_peak));
            
            tooFast_TC{id} = avg_tooFast;
            tooFast_TC_RS{id} = avg_tooFast(RS_cells{id},:);
            tooFast_TC_RL{id} = avg_tooFast(RL_cells{id},:);
            tooFast_resp_all{id} = mean((tooFast_resp-tooFast_base),1);
            tooFast_resp_RS{id} = mean((tooFast_resp(:,RS_cells{id})-tooFast_base(:,RS_cells{id})),1);
            tooFast_resp_RL{id} = mean((tooFast_resp(:,release_resp_cells)-tooFast_base(:,release_resp_cells)),1);
            
            fail_TC{id} = avg_fail;
            fail_TC_RS{id} = avg_fail(RS_cells{id},:);
            fail_TC_RL{id} = avg_fail(RL_cells{id},:);
            fail_TC_trans{id} = avg_fail(trans_cell_ind,:);
            fail_TC_sus{id} = avg_fail(sustain_cell_ind,:);
            fail_TC_mean{id} = mean(avg_fail,1);
            fail_TC_sem{id} = std(avg_fail,1)./sqrt(size(avg_fail,1));
            fail_TC_RS_mean{id} = mean(avg_fail(RS_cells{id},:),1);
            fail_TC_RS_sem{id} = std(avg_fail(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
            fail_TC_RL_mean{id} = mean(avg_fail(release_resp_cells,:),1);
            fail_TC_RL_sem{id} = std(avg_fail(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
            fail_resp_all{id} = mean((fail_resp-fail_base),1);
            fail_resp_mean{id} = mean(fail_resp_all{id},2);
            fail_resp_sem{id} = std(fail_resp_all{id},[],2)./sqrt(ncells{id});
            fail_resp_RS{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),1);
            fail_resp_RS_mean{id} = mean(fail_resp_RS{id},2);
            fail_resp_RS_sem{id} = std(fail_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
            fail_resp_RL{id} = mean((fail_resp(:,release_resp_cells)-fail_base(:,release_resp_cells)),1);
            fail_resp_RL_mean{id} = mean(fail_resp_RL{id},2);
            fail_resp_RL_sem{id} = std(fail_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
            fail_resp_SC{id} = mean((fail_resp(:,success_resp_cells)-fail_base(:,success_resp_cells)),1);
            fail_resp_SC_mean{id} = mean(fail_resp_SC{id},2);
            fail_resp_SC_sem{id} = std(fail_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
            fail_resp_RL_trans{id} = mean((fail_resp(:,trans_cell_ind)-fail_base(:,trans_cell_ind)),1);
            fail_resp_RL_sus{id} = mean((fail_resp(:,sustain_cell_ind)-fail_base(:,sustain_cell_ind)),1);
        end
    
end

total_cells = sum(cell2mat(ncells));
total_resp = sum(cell2mat(tot_resp));
total_RL_resp = sum(cell2mat(tot_RL_resp));

save([out_base, 'cell_count.mat'], 'total_cells', 'total_resp', 'total_RL_resp');

%summary of average response amplitudes
fig=figure;
%resp amplitudes- press/release
subplot(3,2,1)
errorbar(1:size(mouseID,2), cell2mat(release_resp_mean), cell2mat(release_resp_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(press_resp_mean), cell2mat(press_resp_sem), 'oc')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Press/release amplitude')
%resp amplitudes- success/fail
subplot(3,2,2)
errorbar(1:size(mouseID,2), cell2mat(success_resp_mean), cell2mat(success_resp_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(fail_resp_mean), cell2mat(fail_resp_sem), 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Success/fail amplitude')

%resp amplitudes- press/release - driven by press or release
subplot(3,2,3)
errorbar(1:size(mouseID,2), cell2mat(release_resp_RS_mean), cell2mat(release_resp_RS_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(press_resp_RS_mean), cell2mat(press_resp_RS_sem), 'oc')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Press/release: driven')
%resp amplitudes- success/fail - driven press or release
subplot(3,2,4)
errorbar(1:size(mouseID,2), cell2mat(success_resp_RS_mean), cell2mat(success_resp_RS_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(fail_resp_RS_mean), cell2mat(fail_resp_RS_sem), 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Success/fail: driven ')

%resp amplitudes- success/fail - driven by release
subplot(3,2,5)
errorbar(1:size(mouseID,2), cell2mat(success_resp_RL_mean), cell2mat(success_resp_RL_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(fail_resp_RL_mean), cell2mat(fail_resp_RL_sem), 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Success/fail driven by release')
%resp amplitudes- success/fail - driven by fail
subplot(3,2,6)
errorbar(1:size(mouseID,2), cell2mat(success_resp_SC_mean), cell2mat(success_resp_SC_sem), 'ok')
hold on
errorbar(1:size(mouseID,2), cell2mat(fail_resp_SC_mean), cell2mat(fail_resp_SC_sem), 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(mouseID,2)+1])
title('Success/fail driven by success')

supertitle(['Summary of cell response amplitudes'])
saveas(fig, [out_base 'Summary_cell_response_amp.fig']);
print([out_base 'Summary_cell_response_amp.eps'], '-depsc');
print([out_base 'Summary_cell_response_amp.pdf'], '-dpdf');

%scatter of response amplitudes
fig=figure;
x = [-.05:.01:.2];
y = x;
r = [];
p = [];
% col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c'); %hardcoded
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
    0.5 0.4 0];

col_mat_s = repmat([0.5 0.5 0.5], 11, 1);

subplot(3,3,1)
[h_prall, p_prall] = scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(p_prall)])

subplot(3,3,2)
[h_prRS, p_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Release dF/F')
ylabel('Press dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(p_prRS)])
subplot(3,3,3)
[h_prRL, p_prRL] = scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(p_prRL)])

subplot(3,3,4)
[h_sfall, p_sfall] = scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat_s);
hold on 
xlim([-.025 .3]);
ylim([-.005 .25]);
plot(x,y,'-k')
xlabel('Success dF/F')
ylabel('Fail dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(chop(p_sfall,2))])
subplot(3,3,5)
[h_sfRS, p_sfRS] = scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('Fail dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])

subplot(3,3,6)
[h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('Fail dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])

subplot(3,3,7)
[h_stall, p_sfall] = scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('TooFast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(chop(p_sfall,2))])

subplot(3,3,8)
[h_sfRS, p_sfRS] = scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('TooFast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])

subplot(3,3,9)
[h_stRL, p_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('TooFast dF/F')

hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])

supertitle(['Summary of cell response amplitudes'])

saveas(fig, [out_base 'Summary_cell_response_amp_scatter.fig']);
print([out_base 'Summary_cell_response_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_cell_response_amp_scatter.pdf'], '-dpdf');

%new figure for comparing transient and sustain response
fig = figure;
x = [-.05:.01:.3];
y = x;
subplot(1,2,1)
scatter_plot(mouseID, success_resp_RL_trans, fail_resp_RL_trans, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
axis square
xlabel('Success dF/F')
ylabel('Fail dF/F')
title(['Transient response cells n=',num2str(sum(RL_trans_cellNum))])
subplot(1,2,2)
scatter_plot(mouseID, success_resp_RL_sus, fail_resp_RL_sus, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
axis square
xlabel('Success dF/F')
ylabel('Fail dF/F')
title(['Sustain response cells n=',num2str(sum(RL_sustain_cellNum))])
saveas(fig, [out_base 'Summary_transient_sustain_resp_amp.fig']);
print([out_base 'Summary_transient_sustain_resp_amp.eps'], '-depsc');
print([out_base 'Summary_transient_sustain_resp_amp.pdf'], '-dpdf');

%avg scatter of amplitudes by mouse
fig=figure;
x = [-.05:.01:.25];
y = x;
% col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm', 'c');   %HARDCODED
subplot(3,3,1)
scatter_plot(mouseID, release_resp_all, press_resp_all, col_mat);

hold on
plot(x,y,'-k')
xlim([-.025 .3]);
ylim([-.025 .2]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells'])
subplot(3,3,2)
scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Release dF/F')
ylabel('Press dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells'])
subplot(3,3,3)
[h_prRS, p_prRS]=scatter_plot(mouseID, release_resp_RS, press_resp_RS, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Release dF/F')
ylabel('Press dF/F')
title(['Resp cells'])
% subplot(3,3,3)
% scatter_plot(mouseID, release_resp_RL, press_resp_RL, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Release dF/F')
% ylabel('Press dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Release resp cells'])
subplot(3,3,4)
scatter_plot(mouseID, success_resp_all, fail_resp_all, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('Early dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells'])
% subplot(3,3,5)
% scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Fail dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])
subplot(3,3,5)
scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('Early dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells'])

subplot(3,3,6)
[h_sfRL, p_sfRL] = scatter_plot(mouseID, success_resp_RL, fail_resp_RL, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('Early dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells'])
% title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])

subplot(3,3,7)
scatter_plot(mouseID, success_resp_all, tooFast_resp_all, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('TooFast Correct dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells'])

% subplot(3,3,8)
% scatter_plot(mouseID, success_resp_RS, tooFast_resp_RS, col_mat);
% hold on
% plot(x,y,'-k')
% xlim([-.05 .3]);
% ylim([-.05 .3]);
% xlabel('Success dF/F')
% ylabel('Toofast Correct dF/F')
% hold on
% vline(0,'--k')
% hline(0,'--k')
% title(['Resp cells'])

subplot(3,3,8)
scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('Toofast Correct dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells'])

subplot(3,3,9)
[h_stRL, p_stRL] = scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat_s);
hold on 
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Correct dF/F')
ylabel('TooFast Correct dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells'])
% title(['Release resp cells- p = ' num2str(chop(p_stRL,2))])

supertitle(['Summary of cell response amplitudes'])
saveas(fig, [out_base 'Summary_avg_response_amp_scatter.fig']);
print([out_base 'Summary_avg_response_amp_scatter.eps'], '-depsc');
print([out_base 'Summary_avg_response_amp_scatter.pdf'], '-dpdf');

% plot avg response amp scatter for only fail/success and tooFast/success
fig = figure;
subplot(1,2,1)
scatter_plot(mouseID, success_resp_RS, fail_resp_RS, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('Fail dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells'])
axis square;
subplot(1,2,2)
scatter_plot(mouseID, success_resp_RL, tooFast_resp_RL, col_mat);
hold on
plot(x,y,'-k')
xlim([-.05 .3]);
ylim([-.05 .3]);
xlabel('Success dF/F')
ylabel('Toofast Correct dF/F')
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells'])
axis square;

supertitle(['Summary of cell response amplitudes'])
saveas(fig, [out_base 'Summary_avg_response_amp_scatter_2.fig']);
print([out_base 'Summary_avg_response_amp_scatter_2.eps'], '-depsc');
print([out_base 'Summary_avg_response_amp_scatter_2.pdf'], '-dpdf');
%average timecourse for expts

fig=figure;
for id = 1:size(mouseID,2)
    subplot(5,5,id)
    tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
    shadedErrorBar(tt, success_TC_mean{id},success_TC_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_mean{id},fail_TC_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_mean{id},press_TC_sem{id}, 'c');
    xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date{id} ' ' mouseID{id}])
end
supertitle(['Summary of all cell timecourses'])
saveas(fig, [out_base 'Summary_allexptTCs_allcells.fig']);
print([out_base 'Summary_allexptTCs_allcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_allcells.pdf'], '-dpdf');

fig=figure;
for id = 1:size(mouseID,2)
    subplot(5,5,id)
    tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
    shadedErrorBar(tt, success_TC_RS_mean{id},success_TC_RS_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_RS_mean{id},fail_TC_RS_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_RS_mean{id},press_TC_RS_sem{id}, 'c');
    xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date{id} ' ' mouseID{id}])
end
supertitle(['Summary of responsive cell timecourses'])
saveas(fig, [out_base 'Summary_allexptTCs_respcells.fig']);
print([out_base 'Summary_allexptTCs_respcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_respcells.pdf'], '-dpdf');

fig=figure;
for id = 1:size(mouseID,2)
    subplot(5,5,id)
    tt =((-pre_frames{id}:post_frames{id}).*double(TC_ifi{id}))./1000;
    shadedErrorBar(tt, success_TC_RL_mean{id},success_TC_RL_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_RL_mean{id},fail_TC_RL_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_RL_mean{id},press_TC_RL_sem{id}, 'c');
    xlim([-pre_frames{id}.*double(TC_ifi{id})./1000 post_frames{id}.*double(TC_ifi{id})./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date{id} ' ' mouseID{id}])
end
supertitle(['Summary of release responsive cell timecourses'])
saveas(fig, [out_base 'Summary_allexptTCs_relcells.fig']);
print([out_base 'Summary_allexptTCs_relcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_relcells.pdf'], '-dpdf');

%% commented for now until a decision is made on how to average across experiments with different acquisition rates
%averaging across all cells- specific to different acquisition rates

frame_size = cell2mat(cellfun(@size, success_TC, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

success_TC_all = interp_frame(success_TC, max_frame);
fail_TC_all    = interp_frame(fail_TC, max_frame);
press_TC_all   = interp_frame(press_TC, max_frame);
tooFast_TC_all = interp_frame(tooFast_TC, max_frame);

success_TC_all_RS = interp_frame(success_TC_RS, max_frame);
fail_TC_all_RS    = interp_frame(fail_TC_RS, max_frame);
press_TC_all_RS   = interp_frame(press_TC_RS, max_frame);
tooFast_TC_all_RS = interp_frame(tooFast_TC_RS, max_frame);

success_TC_all_RL = interp_frame(success_TC_RL, max_frame);
fail_TC_all_RL    = interp_frame(fail_TC_RL, max_frame);
press_TC_all_RL   = interp_frame(press_TC_RL, max_frame);
tooFast_TC_all_RL = interp_frame(tooFast_TC_RL, max_frame);

tt =(-max(cell2mat(pre_frames)):max(cell2mat(post_frames))).*double(min(cell2mat(TC_ifi)))./1000;
fig=figure;
subplot(3,1,1)
shadedErrorBar(tt, mean(success_TC_all,1), std(success_TC_all,[],1)./sqrt(size(success_TC_all,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all,1), std(fail_TC_all,[],1)./sqrt(size(fail_TC_all,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all,1), std(press_TC_all,[],1)./sqrt(size(press_TC_all,1)), 'c');
shadedErrorBar(tt, mean(tooFast_TC_all,1), std(tooFast_TC_all,[],1)./sqrt(size(tooFast_TC_all,1)), 'g');
title(['All cells- n = ' num2str(size(press_TC_all,1))])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,2)
shadedErrorBar(tt, mean(success_TC_all_RS,1), std(success_TC_all_RS,[],1)./sqrt(size(success_TC_all_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RS,1), std(fail_TC_all_RS,[],1)./sqrt(size(fail_TC_all_RS,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RS,1), std(press_TC_all_RS,[],1)./sqrt(size(press_TC_all_RS,1)), 'c');
shadedErrorBar(tt, mean(tooFast_TC_all_RS,1), std(tooFast_TC_all_RS,[],1)./sqrt(size(tooFast_TC_all_RS,1)), 'g');
title(['Responsive cells- n = ' num2str(size(press_TC_all_RS,1))])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,3)
shadedErrorBar(tt, mean(success_TC_all_RL,1), std(success_TC_all_RL,[],1)./sqrt(size(success_TC_all_RL,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RL,1), std(fail_TC_all_RL,[],1)./sqrt(size(fail_TC_all_RL,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RL,1), std(press_TC_all_RL,[],1)./sqrt(size(press_TC_all_RL,1)), 'c');
shadedErrorBar(tt, mean(tooFast_TC_all_RL,1), std(tooFast_TC_all_RL,[],1)./sqrt(size(tooFast_TC_all_RL,1)), 'g');
title(['Release responsive cells- n = ' num2str(size(press_TC_all_RL,1))])
xlabel('Time (ms)')
ylabel('dF/F')
% supertitle('Average all cells collected at 30 Hz')
supertitle(['Timecourses of average all cells across experiments Black- Correct; Red- Early; Cyan- Press; Green- TooFast']);
saveas(fig, [out_base 'Summary_avgexpt_TCs.fig']);
print([out_base 'Summary_avgexpt_TCs.eps'], '-depsc');
print([out_base 'Summary_avgexpt_TCs.pdf'], '-dpdf');

%% plot TC for sustain and transient cells
success_TC_sus = success_TC_sus(~cellfun('isempty', success_TC_sus));
fail_TC_sus = fail_TC_sus(~cellfun('isempty', fail_TC_sus));

frame_size = cell2mat(cellfun(@size, success_TC_trans, 'UniformOutput', 0));
max_frame = max(frame_size(2:2:end));

success_TC_all_sus = interp_frame(success_TC_sus, max_frame);
fail_TC_all_sus    = interp_frame(fail_TC_sus, max_frame);

success_TC_all_trans = interp_frame(success_TC_trans, max_frame);
fail_TC_all_trans    = interp_frame(fail_TC_trans, max_frame);

fig = figure;
tt =(-max(cell2mat(pre_frames)):max(cell2mat(post_frames))).*double(min(cell2mat(TC_ifi)))./1000;
subplot(1,2,1)
shadedErrorBar(tt, mean(success_TC_all_trans,1), std(success_TC_all_trans,[],1)./sqrt(size(success_TC_all_trans,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_trans,1), std(fail_TC_all_trans,[],1)./sqrt(size(fail_TC_all_trans,1)), 'r');
ylim([-0.03 0.15])
xlim([-0.5,2.5])
xlabel('time (ms)')
ylabel('dF/F')
title(['transient cells n=',num2str(sum(RL_trans_cellNum))])

subplot(1,2,2)
shadedErrorBar(tt, mean(success_TC_all_sus,1), std(success_TC_all_sus,[],1)./sqrt(size(success_TC_all_sus,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_sus,1), std(fail_TC_all_sus,[],1)./sqrt(size(fail_TC_all_sus,1)), 'r');
ylim([-0.03 0.15])
xlim([-0.5,2.5])
xlabel('time (ms)')
ylabel('dF/F')
title(['sustain cells n=',num2str(sum(RL_sustain_cellNum))])


supertitle(['Timecourses of sustain and transient cells'])
saveas(fig, [out_base 'Summary_avgexpt_transient_sus_TCs.fig']);
print([out_base 'Summary_avgexpt_transient_sus_TCs.eps'], '-depsc');
print([out_base 'Summary_avgexpt_transient_sus_TCs.pdf'], '-dpdf');
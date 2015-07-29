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
    load([dest_sub '_cell_TCs.mat']);
    load([dest_sub '_cell_resp.mat']);
    load([dest_sub '_cell_categories.mat']);
    load([dest_sub '_release_movies.mat'])
    
    ncells(id) = size(press_resp,2);
    RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells]);
    tot_resp(id) = length(RS_cells{id});
    
    TC_length(id) = size(avg_release,2);
    TC_ifi(id) = ifi;
    pre_frames(id) = pre_release_frames;
    post_frames(id) = post_release_frames;
    
    RL_cells{id} = release_resp_cells;
    
    release_resp_all{id} = mean((release_resp-release_base),1);
    release_resp_mean(id) = mean(release_resp_all{id},2);
    release_resp_sem(id) = std(release_resp_all{id},[],2)./sqrt(ncells(id));
    release_resp_RL{id} = mean((release_resp(:,release_resp_cells)-release_base(:,release_resp_cells)),1);
    release_resp_RL_mean(id) = mean(release_resp_RL{id},2);
    release_resp_RL_sem(id) = std(release_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
    release_resp_RS{id} = mean((release_resp(:,RS_cells{id})-release_base(:,RS_cells{id})),1);
    release_resp_RS_mean(id) = mean(release_resp_RS{id},2);
    release_resp_RS_sem(id) = std(release_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    
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
    press_resp_mean(id) = mean(press_resp_all{id},2);
    press_resp_sem(id) = std(press_resp_all{id},[],2)./sqrt(ncells(id));
    press_resp_RL{id} = mean((press_resp(:,release_resp_cells)-press_base(:,release_resp_cells)),1);
    press_resp_RL_mean(id) = mean(press_resp_RL{id},2);
    press_resp_RL_sem(id) = std(press_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
    press_resp_RS{id} = mean((press_resp(:,RS_cells{id})-press_base(:,RS_cells{id})),1);
    press_resp_RS_mean(id) = mean(press_resp_RS{id},2);
    press_resp_RS_sem(id) = std(press_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    
    success_TC{id} = avg_success;
    success_TC_RS{id} = avg_success(RS_cells{id},:);
    success_TC_RL{id} = avg_success(RL_cells{id},:);
    success_TC_mean{id} = mean(avg_success,1);
    success_TC_sem{id} = std(avg_success,1)./sqrt(size(avg_success,1));
    success_TC_RS_mean{id} = mean(avg_success(RS_cells{id},:),1);
    success_TC_RS_sem{id} = std(avg_success(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
    success_TC_RL_mean{id} = mean(avg_success(release_resp_cells,:),1);
    success_TC_RL_sem{id} = std(avg_success(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
    success_resp_all{id} = mean((success_resp-success_base),1);
    success_resp_mean(id) = mean(success_resp_all{id},2);
    success_resp_sem(id) = std(success_resp_all{id},[],2)./sqrt(ncells(id));
    success_resp_RS{id} = mean((success_resp(:,RS_cells{id})-success_base(:,RS_cells{id})),1);
    success_resp_RS_mean(id) = mean(success_resp_RS{id},2);
    success_resp_RS_sem(id) = std(success_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    success_resp_RL{id} = mean((success_resp(:,release_resp_cells)-success_base(:,release_resp_cells)),1);
    success_resp_RL_mean(id) = mean(success_resp_RL{id},2);
    success_resp_RL_sem(id) = std(success_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
    success_resp_SC{id} = mean((success_resp(:,success_resp_cells)-success_base(:,success_resp_cells)),1);
    success_resp_SC_mean(id) = mean(success_resp_SC{id},2);
    success_resp_SC_sem(id) = std(success_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
    
    fail_TC{id} = avg_fail;
    fail_TC_RS{id} = avg_fail(RS_cells{id},:);
    fail_TC_RL{id} = avg_fail(RL_cells{id},:);
    fail_TC_mean{id} = mean(avg_fail,1);
    fail_TC_sem{id} = std(avg_fail,1)./sqrt(size(avg_fail,1));
    fail_TC_RS_mean{id} = mean(avg_fail(RS_cells{id},:),1);
    fail_TC_RS_sem{id} = std(avg_fail(RS_cells{id},:),1)./sqrt(size(RS_cells{id},2));
    fail_TC_RL_mean{id} = mean(avg_fail(release_resp_cells,:),1);
    fail_TC_RL_sem{id} = std(avg_fail(release_resp_cells,:),1)./sqrt(size(release_resp_cells,2));
    fail_resp_all{id} = mean((fail_resp-fail_base),1);
    fail_resp_mean(id) = mean(fail_resp_all{id},2);
    fail_resp_sem(id) = std(fail_resp_all{id},[],2)./sqrt(ncells(id));
    fail_resp_RS{id} = mean((fail_resp(:,RS_cells{id})-fail_base(:,RS_cells{id})),1);
    fail_resp_RS_mean(id) = mean(fail_resp_RS{id},2);
    fail_resp_RS_sem(id) = std(fail_resp_RS{id},[],2)./sqrt(length(RS_cells{id}));
    fail_resp_RL{id} = mean((fail_resp(:,release_resp_cells)-fail_base(:,release_resp_cells)),1);
    fail_resp_RL_mean(id) = mean(fail_resp_RL{id},2);
    fail_resp_RL_sem(id) = std(fail_resp_RL{id},[],2)./sqrt(length(release_resp_cells));
    fail_resp_SC{id} = mean((fail_resp(:,success_resp_cells)-fail_base(:,success_resp_cells)),1);
    fail_resp_SC_mean(id) = mean(fail_resp_SC{id},2);
    fail_resp_SC_sem(id) = std(fail_resp_SC{id},[],2)./sqrt(length(success_resp_cells));
end
%summary of average response amplitudes
figure; 
%resp amplitudes- press/release
subplot(3,2,1)
errorbar(1:size(date_mat,1), release_resp_mean, release_resp_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), press_resp_mean, press_resp_sem, 'oc')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Press/release amplitude')
%resp amplitudes- success/fail
subplot(3,2,2)
errorbar(1:size(date_mat,1), success_resp_mean, success_resp_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), fail_resp_mean, fail_resp_sem, 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Success/fail amplitude')

%resp amplitudes- press/release - driven by press or release
subplot(3,2,3)
errorbar(1:size(date_mat,1), release_resp_RS_mean, release_resp_RS_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), press_resp_RS_mean, press_resp_RS_sem, 'oc')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Press/release: driven')
%resp amplitudes- success/fail - driven press or release
subplot(3,2,4)
errorbar(1:size(date_mat,1), success_resp_RS_mean, success_resp_RS_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), fail_resp_RS_mean, fail_resp_RS_sem, 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Success/fail: driven ')

%resp amplitudes- success/fail - driven by release
subplot(3,2,5)
errorbar(1:size(date_mat,1), success_resp_RL_mean, success_resp_RL_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), fail_resp_RL_mean, fail_resp_RL_sem, 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Success/fail driven by release')
%resp amplitudes- success/fail - driven by fail
subplot(3,2,6)
errorbar(1:size(date_mat,1), success_resp_SC_mean, success_resp_SC_sem, 'ok')
hold on
errorbar(1:size(date_mat,1), fail_resp_SC_mean, fail_resp_SC_sem, 'or')
xlabel('Expt #')
ylabel('dF/F')
xlim([0 size(date_mat,1)+1])
title('Success/fail driven by success')


suptitle(['Summary of cell response amplitudes'])
    print([out_base 'Summary_cell_response_amp.eps'], '-depsc');
    print([out_base 'Summary_cell_response_amp.pdf'], '-dpdf');
    
%scatter of response amplitudes
figure;
x = [-.05:.01:.2];
y = x;
r = [];
p = [];
col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm');
subplot(2,3,1)
for id = 1:length(date_mat)
    scatter(release_resp_all{id}, press_resp_all{id}, col_mat(id,:))
    r = [r release_resp_all{id}];
    p = [p press_resp_all{id}];
    hold on
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Release dF/F')
ylabel('Press dF/F')
[h_prall p_prall] = ttest(r,p);
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(p_prall)])
subplot(2,3,2)
r = [];
p = [];
for id = 1:length(date_mat)
    scatter(release_resp_RS{id}, press_resp_RS{id}, col_mat(id,:))
    hold on
    r = [r release_resp_RS{id}];
    p = [p press_resp_RS{id}];
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Release dF/F')
ylabel('Press dF/F')
[h_prRS p_prRS] = ttest(r,p);
hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(p_prRS)])
subplot(2,3,3)
r = [];
p = [];
for id = 1:length(date_mat)
    scatter(release_resp_RL{id}, press_resp_RL{id}, col_mat(id,:))
    hold on
    r = [r release_resp_RL{id}];
    p = [p press_resp_RL{id}];
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Release dF/F')
ylabel('Press dF/F')
[h_prRL p_prRL] = ttest(r,p);
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(p_prRL)])
subplot(2,3,4)
s = [];
f = [];
for id = 1:length(date_mat)
    scatter(success_resp_all{id}, fail_resp_all{id}, col_mat(id,:))
    hold on
    s = [s success_resp_all{id}];
    f = [f fail_resp_all{id}];
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Success dF/F')
ylabel('Fail dF/F')
[h_sfall p_sfall] = ttest(s,f);
hold on
vline(0,'--k')
hline(0,'--k')
title(['All cells- p = ' num2str(chop(p_sfall,2))])
subplot(2,3,5)
s = [];
f = [];
for id = 1:length(date_mat)
    scatter(success_resp_RS{id}, fail_resp_RS{id}, col_mat(id,:))
    hold on
    s = [s success_resp_RS{id}];
    f = [f fail_resp_RS{id}];
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Success dF/F')
ylabel('Fail dF/F')
[h_sfRS p_sfRS] = ttest(s,f);
hold on
vline(0,'--k')
hline(0,'--k')
title(['Resp cells- p = ' num2str(chop(p_sfRS,2))])
subplot(2,3,6)
s = [];
f = [];
for id = 1:length(date_mat)
    scatter(success_resp_RL{id}, fail_resp_RL{id}, col_mat(id,:))
    hold on
    s = [s success_resp_RL{id}];
    f = [f fail_resp_RL{id}];
end
plot(x,y,'-k')
xlim([-.05 .1]);
ylim([-.05 .1]);
xlabel('Success dF/F')
ylabel('Fail dF/F')
[h_sfRL p_sfRL] = ttest(s,f);
hold on
vline(0,'--k')
hline(0,'--k')
title(['Release resp cells- p = ' num2str(chop(p_sfRL,2))])
suptitle(['Summary of cell response amplitudes- Red: img24; Blue: img25; Green: img28; Magenta: img27'])
    print([out_base 'Summary_cell_response_amp_scatter.eps'], '-depsc');
    print([out_base 'Summary_cell_response_amp_scatter.pdf'], '-dpdf');
    
%average timecourse for expts

figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    tt =((-pre_frames(id):post_frames(id)).*double(TC_ifi(id)))./1000;
    shadedErrorBar(tt, success_TC_mean{id},success_TC_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_mean{id},fail_TC_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_mean{id},press_TC_sem{id}, 'c');
    xlim([-pre_frames(id).*double(TC_ifi(id))./1000 post_frames(id).*double(TC_ifi(id))./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date_mat(id,:) ' ' mouse_mat(id,:)])
end
suptitle(['Summary of all cell timecourses'])
print([out_base 'Summary_allexptTCs_allcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_allcells.pdf'], '-dpdf');

figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    tt =((-pre_frames(id):post_frames(id)).*double(TC_ifi(id)))./1000;
    shadedErrorBar(tt, success_TC_RS_mean{id},success_TC_RS_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_RS_mean{id},fail_TC_RS_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_RS_mean{id},press_TC_RS_sem{id}, 'c');
    xlim([-pre_frames(id).*double(TC_ifi(id))./1000 post_frames(id).*double(TC_ifi(id))./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date_mat(id,:) ' ' mouse_mat(id,:)])
end
suptitle(['Summary of responsive cell timecourses'])
print([out_base 'Summary_allexptTCs_respcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_respcells.pdf'], '-dpdf');

figure;
for id = 1:length(date_mat)
    subplot(2,3,id)
    tt =((-pre_frames(id):post_frames(id)).*double(TC_ifi(id)))./1000;
    shadedErrorBar(tt, success_TC_RL_mean{id},success_TC_RL_sem{id}, 'k');
    hold on
    shadedErrorBar(tt, fail_TC_RL_mean{id},fail_TC_RL_sem{id}, 'r');
    hold on
    shadedErrorBar(tt, press_TC_RL_mean{id},press_TC_RL_sem{id}, 'c');
    xlim([-pre_frames(id).*double(TC_ifi(id))./1000 post_frames(id).*double(TC_ifi(id))./1000])
    xlabel('Time (ms)')
    ylabel('dF/F')
    title([date_mat(id,:) ' ' mouse_mat(id,:)])
end
suptitle(['Summary of release responsive cell timecourses'])
print([out_base 'Summary_allexptTCs_relcells.eps'], '-depsc');
print([out_base 'Summary_allexptTCs_relcells.pdf'], '-dpdf');

%averaging across all cells- specific to different acquisition rates
success_TC_all = [];
fail_TC_all = [];
press_TC_all = [];
success_TC_all_RS = [];
fail_TC_all_RS = [];
press_TC_all_RS = [];
success_TC_all_RL = [];
fail_TC_all_RL = [];
press_TC_all_RL = [];
for id = 1:4
    success_TC_all = [success_TC_all; success_TC{id}];
    fail_TC_all = [fail_TC_all; fail_TC{id}];
    press_TC_all = [press_TC_all; press_TC{id}];
    success_TC_all_RS = [success_TC_all_RS; success_TC_RS{id}];
    fail_TC_all_RS = [fail_TC_all_RS; fail_TC_RS{id}];
    press_TC_all_RS = [press_TC_all_RS; press_TC_RS{id}];
    success_TC_all_RL = [success_TC_all_RL; success_TC_RL{id}];
    fail_TC_all_RL = [fail_TC_all_RL; fail_TC_RL{id}];
    press_TC_all_RL = [press_TC_all_RL; press_TC_RL{id}];
end
tt =((-pre_frames(1):post_frames(1)).*double(TC_ifi(1)))./1000;
figure;
subplot(3,1,1)
shadedErrorBar(tt, mean(success_TC_all,1), std(success_TC_all,[],1)./sqrt(size(success_TC_all,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all,1), std(fail_TC_all,[],1)./sqrt(size(fail_TC_all,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all,1), std(press_TC_all,[],1)./sqrt(size(press_TC_all,1)), 'c');
title(['All cells- n = ' num2str(size(press_TC_all,1))])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,2)
shadedErrorBar(tt, mean(success_TC_all_RS,1), std(success_TC_all_RS,[],1)./sqrt(size(success_TC_all_RS,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RS,1), std(fail_TC_all_RS,[],1)./sqrt(size(fail_TC_all_RS,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RS,1), std(press_TC_all_RS,[],1)./sqrt(size(press_TC_all_RS,1)), 'c');
title(['Responsive cells- n = ' num2str(size(press_TC_all_RS,1))])
xlabel('Time (ms)')
ylabel('dF/F')

subplot(3,1,3)
shadedErrorBar(tt, mean(success_TC_all_RL,1), std(success_TC_all_RL,[],1)./sqrt(size(success_TC_all_RL,1)), 'k');
hold on;
shadedErrorBar(tt, mean(fail_TC_all_RL,1), std(fail_TC_all_RL,[],1)./sqrt(size(fail_TC_all_RL,1)), 'r');
hold on
shadedErrorBar(tt, mean(press_TC_all_RL,1), std(press_TC_all_RL,[],1)./sqrt(size(press_TC_all_RL,1)), 'c');
title(['Release responsive cells- n = ' num2str(size(press_TC_all_RL,1))])
xlabel('Time (ms)')
ylabel('dF/F')
suptitle('Average all cells collected at 15 Hz')
suptitle(['Summary of release responsive cell timecourses'])
print([out_base 'Summary_15HzTCs.eps'], '-depsc');
print([out_base 'Summary_15HzTCs.pdf'], '-dpdf');

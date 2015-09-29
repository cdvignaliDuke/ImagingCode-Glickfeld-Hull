close all
rc = behavConstsAV;
av = behavParamsAV;
mouse_str = [];
for imouse = 1:size(av,2)
    mouse_str = [mouse_str 'i' num2str(av(imouse).mouse) '_'];  
end
fnout = fullfile(rc.caOutputDir, [date '_' mouse_str]);

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_pre_all = [];
    FAs_pre_all = [];
    hits_post_all = [];
    FAs_post_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_pre = mouse(imouse).expt(iexp).align(2).av(1).outcome(1).pre_resp(1,:);
        FAs_pre = mouse(imouse).expt(iexp).align(2).av(1).outcome(3).pre_resp(1,:);
        scatter(hits_pre, FAs_pre,['.'  col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_pre,2), mean(FAs_pre,2), std(hits_pre,[],2)./sqrt(size(hits_pre,2)), std(FAs_pre,[],2)./sqrt(size(FAs_pre,2)), {['o' col], col, col});
        hold on
        hits_pre_all = [hits_pre_all hits_pre];
        FAs_pre_all = [FAs_pre_all FAs_pre];
        subplot(3,3,4)
        hits_post = mouse(imouse).expt(iexp).align(2).av(1).outcome(1).post_resp(1,:)-mouse(imouse).expt(iexp).align(2).av(1).outcome(1).pre_resp(1,:);
        FAs_post = mouse(imouse).expt(iexp).align(2).av(1).outcome(3).post_resp(1,:)-mouse(imouse).expt(iexp).align(2).av(1).outcome(3).pre_resp(1,:);
        scatter(hits_post, FAs_post,['.'  col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_post,2), mean(FAs_post,2), std(hits_post,[],2)./sqrt(size(hits_post,2)), std(FAs_post,[],2)./sqrt(size(FAs_post,2)), {['o' col], col, col});
        hold on
        hits_post_all = [hits_post_all hits_post];
        FAs_post_all = [FAs_post_all FAs_post];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_pre_all,2), mean(FAs_pre_all,2), std(hits_pre_all,[],2)./sqrt(size(hits_pre_all,2)), std(FAs_pre_all,[],2)./sqrt(size(FAs_pre_all,2)), {['o' col], col, col});
    [h_pre(imouse), p_pre(imouse)] = ttest(hits_pre_all,FAs_pre_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_post_all,2), mean(FAs_post_all,2), std(hits_post_all,[],2)./sqrt(size(hits_post_all,2)), std(FAs_post_all,[],2)./sqrt(size(FAs_post_all,2)), {['o' col], col, col});
    [h_post(imouse), p_post(imouse)] = ttest(hits_post_all,FAs_post_all);
    hold on
    n(imouse) = length(hits_pre_all);
end
for i = 1:6
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.1 0.1])
ylim([-0.1 0.1])
title('Pre-release resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Pre-release resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Pre-release resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_pre(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_pre(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.1 0.1])
ylim([-0.1 0.1])
title('Post-release resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Post-release resp- expt avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Post-release resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_post(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_post(2),3))], 'Color', av(2).col_str);
suptitle('Visual trials- Hits vs FAs- release align- all cells')
print([fnout 'release_align_SE_V.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_pre_all = [];
    FAs_pre_all = [];
    hits_post_all = [];
    FAs_post_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_pre = mouse(imouse).expt(iexp).align(2).av(2).outcome(1).pre_resp(1,:);
        FAs_pre = mouse(imouse).expt(iexp).align(2).av(2).outcome(3).pre_resp(1,:);
        scatter(hits_pre, FAs_pre,['.'  col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_pre,2), mean(FAs_pre,2), std(hits_pre,[],2)./sqrt(size(hits_pre,2)), std(FAs_pre,[],2)./sqrt(size(FAs_pre,2)), {['o' col], col, col});
        hold on
        hits_pre_all = [hits_pre_all hits_pre];
        FAs_pre_all = [FAs_pre_all FAs_pre];
        subplot(3,3,4)
        hits_post = mouse(imouse).expt(iexp).align(2).av(2).outcome(1).post_resp(1,:)-mouse(imouse).expt(iexp).align(2).av(2).outcome(1).pre_resp(1,:);
        FAs_post = mouse(imouse).expt(iexp).align(2).av(2).outcome(3).post_resp(1,:)-mouse(imouse).expt(iexp).align(2).av(2).outcome(3).pre_resp(1,:);
        scatter(hits_post, FAs_post,['.'  col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_post,2), mean(FAs_post,2), std(hits_post,[],2)./sqrt(size(hits_post,2)), std(FAs_post,[],2)./sqrt(size(FAs_post,2)), {['o' col], col, col});
        hold on
        hits_post_all = [hits_post_all hits_post];
        FAs_post_all = [FAs_post_all FAs_post];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_pre_all,2), mean(FAs_pre_all,2), std(hits_pre_all,[],2)./sqrt(size(hits_pre_all,2)), std(FAs_pre_all,[],2)./sqrt(size(FAs_pre_all,2)), {['o' col], col, col});
    [h_pre(imouse), p_pre(imouse)] = ttest(hits_pre_all,FAs_pre_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_post_all,2), mean(FAs_post_all,2), std(hits_post_all,[],2)./sqrt(size(hits_post_all,2)), std(FAs_post_all,[],2)./sqrt(size(FAs_post_all,2)), {['o' col], col, col});
    [h_post(imouse), p_post(imouse)] = ttest(hits_post_all,FAs_post_all);
    hold on
    n(imouse) = length(hits_pre_all);
end
for i = 1:6
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.1 0.1])
ylim([-0.1 0.1])
title('Pre-release resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Pre-release resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Pre-release resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_pre(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_pre(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.1 0.1])
ylim([-0.1 0.1])
title('Post-release resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Post-release resp- expt avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Post-release resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_post(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_post(2),3))], 'Color', av(2).col_str);
suptitle('Auditory trials- Hits vs FAs- release align- all cells')
print([fnout 'release_align_SE_A.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    FAs_trans_all = [];
    FAs_mid_all = [];
    FAs_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(4).av(1).outcome(1).trans_resp(1,:);
        FAs_trans = mouse(imouse).expt(iexp).align(4).av(1).outcome(3).trans_resp(1,:);
        scatter(hits_trans, FAs_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(FAs_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(FAs_trans,[],2)./sqrt(size(FAs_trans,2)), {['o' col], col, col});
        hold on
        subplot(3,3,4)
        hits_mid = mouse(imouse).expt(iexp).align(4).av(1).outcome(1).mid_resp(1,:);
        FAs_mid = mouse(imouse).expt(iexp).align(4).av(1).outcome(3).mid_resp(1,:);
        scatter(hits_mid, FAs_mid, ['.' col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_mid,2), mean(FAs_mid,2), std(hits_mid,[],2)./sqrt(size(hits_mid,2)), std(FAs_mid,[],2)./sqrt(size(FAs_mid,2)), {['o' col], col, col});
        hold on
        subplot(3,3,7)
        hits_late = mouse(imouse).expt(iexp).align(4).av(1).outcome(1).late_resp(1,:);
        FAs_late = mouse(imouse).expt(iexp).align(4).av(1).outcome(3).late_resp(1,:);
        scatter(hits_late, FAs_late, ['.' col]);
        hold on
        subplot(3,3,8)
        errorbarxy(mean(hits_late,2), mean(FAs_late,2), std(hits_late,[],2)./sqrt(size(hits_late,2)), std(FAs_late,[],2)./sqrt(size(FAs_late,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        FAs_trans_all = [FAs_trans_all FAs_trans];
        FAs_mid_all = [FAs_mid_all FAs_mid];
        FAs_late_all = [FAs_late_all FAs_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(FAs_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(FAs_trans_all,[],2)./sqrt(size(FAs_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,FAs_trans_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_mid_all,2), mean(FAs_mid_all,2), std(hits_mid_all,[],2)./sqrt(size(hits_mid_all,2)), std(FAs_mid_all,[],2)./sqrt(size(FAs_mid_all,2)), {['o' col], col, col});
    [h_mid(imouse), p_mid(imouse)] = ttest(hits_mid_all,FAs_mid_all);
    hold on
    subplot(3,3,9)
    errorbarxy(mean(hits_late_all,2), mean(FAs_late_all,2), std(hits_late_all,[],2)./sqrt(size(hits_late_all,2)), std(FAs_late_all,[],2)./sqrt(size(FAs_late_all,2)), {['o' col], col, col});
    [h_late(imouse), p_late(imouse)] = ttest(hits_late_all,FAs_late_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:9
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Mid resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_mid(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_mid(2),3))], 'Color', av(2).col_str);
subplot(3,3,7)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Late resp')
subplot(3,3,8)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- avg')
subplot(3,3,9)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_late(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_late(2),3))], 'Color', av(2).col_str);
suptitle('Visual trials- Hits vs FAs- prev stim align- all cells')
print([fnout 'prevStim_align_SE_V.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    FAs_trans_all = [];
    FAs_mid_all = [];
    FAs_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(4).av(2).outcome(1).trans_resp(1,:);
        FAs_trans = mouse(imouse).expt(iexp).align(4).av(2).outcome(3).trans_resp(1,:);
        scatter(hits_trans, FAs_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(FAs_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(FAs_trans,[],2)./sqrt(size(FAs_trans,2)), {['o' col], col, col});
        hold on
        subplot(3,3,4)
        hits_mid = mouse(imouse).expt(iexp).align(4).av(2).outcome(1).mid_resp(1,:);
        FAs_mid = mouse(imouse).expt(iexp).align(4).av(2).outcome(3).mid_resp(1,:);
        scatter(hits_mid, FAs_mid, ['.' col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_mid,2), mean(FAs_mid,2), std(hits_mid,[],2)./sqrt(size(hits_mid,2)), std(FAs_mid,[],2)./sqrt(size(FAs_mid,2)), {['o' col], col, col});
        hold on
        subplot(3,3,7)
        hits_late = mouse(imouse).expt(iexp).align(4).av(2).outcome(1).late_resp(1,:);
        FAs_late = mouse(imouse).expt(iexp).align(4).av(2).outcome(3).late_resp(1,:);
        scatter(hits_late, FAs_late, ['.' col]);
        hold on
        subplot(3,3,8)
        errorbarxy(mean(hits_late,2), mean(FAs_late,2), std(hits_late,[],2)./sqrt(size(hits_late,2)), std(FAs_late,[],2)./sqrt(size(FAs_late,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        FAs_trans_all = [FAs_trans_all FAs_trans];
        FAs_mid_all = [FAs_mid_all FAs_mid];
        FAs_late_all = [FAs_late_all FAs_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(FAs_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(FAs_trans_all,[],2)./sqrt(size(FAs_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,FAs_trans_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_mid_all,2), mean(FAs_mid_all,2), std(hits_mid_all,[],2)./sqrt(size(hits_mid_all,2)), std(FAs_mid_all,[],2)./sqrt(size(FAs_mid_all,2)), {['o' col], col, col});
    [h_mid(imouse), p_mid(imouse)] = ttest(hits_mid_all,FAs_mid_all);
    hold on
    subplot(3,3,9)
    errorbarxy(mean(hits_late_all,2), mean(FAs_late_all,2), std(hits_late_all,[],2)./sqrt(size(hits_late_all,2)), std(FAs_late_all,[],2)./sqrt(size(FAs_late_all,2)), {['o' col], col, col});
    [h_late(imouse), p_late(imouse)] = ttest(hits_late_all,FAs_late_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:9
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Mid resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_mid(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_mid(2),3))], 'Color', av(2).col_str);
subplot(3,3,7)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Late resp')
subplot(3,3,8)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- avg')
subplot(3,3,9)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_late(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_late(2),3))], 'Color', av(2).col_str);
suptitle('Auditory trials- Hits vs FAs- prev stim align- all cells')
print([fnout 'prevStim_align_SE_A.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    miss_trans_all = [];
    miss_mid_all = [];
    miss_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(4).av(1).outcome(5).trans_resp(1,:);
        miss_trans = mouse(imouse).expt(iexp).align(4).av(1).outcome(2).trans_resp(1,:);
        scatter(hits_trans, miss_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(miss_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(miss_trans,[],2)./sqrt(size(miss_trans,2)), {['o' col], col, col});
        hold on
        subplot(3,3,4)
        hits_mid = mouse(imouse).expt(iexp).align(4).av(1).outcome(5).mid_resp(1,:);
        miss_mid = mouse(imouse).expt(iexp).align(4).av(1).outcome(2).mid_resp(1,:);
        scatter(hits_mid, miss_mid, ['.' col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_mid,2), mean(miss_mid,2), std(hits_mid,[],2)./sqrt(size(hits_mid,2)), std(miss_mid,[],2)./sqrt(size(miss_mid,2)), {['o' col], col, col});
        hold on
        subplot(3,3,7)
        hits_late = mouse(imouse).expt(iexp).align(4).av(1).outcome(5).late_resp(1,:);
        miss_late = mouse(imouse).expt(iexp).align(4).av(1).outcome(2).late_resp(1,:);
        scatter(hits_late, miss_late, ['.' col]);
        hold on
        subplot(3,3,8)
        errorbarxy(mean(hits_late,2), mean(miss_late,2), std(hits_late,[],2)./sqrt(size(hits_late,2)), std(miss_late,[],2)./sqrt(size(miss_late,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        miss_trans_all = [miss_trans_all miss_trans];
        miss_mid_all = [miss_mid_all miss_mid];
        miss_late_all = [miss_late_all miss_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(miss_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(miss_trans_all,[],2)./sqrt(size(miss_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,miss_trans_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_mid_all,2), mean(miss_mid_all,2), std(hits_mid_all,[],2)./sqrt(size(hits_mid_all,2)), std(miss_mid_all,[],2)./sqrt(size(miss_mid_all,2)), {['o' col], col, col});
    [h_mid(imouse), p_mid(imouse)] = ttest(hits_mid_all,miss_mid_all);
    hold on
    subplot(3,3,9)
    errorbarxy(mean(hits_late_all,2), mean(miss_late_all,2), std(hits_late_all,[],2)./sqrt(size(hits_late_all,2)), std(miss_late_all,[],2)./sqrt(size(miss_late_all,2)), {['o' col], col, col});
    [h_late(imouse), p_late(imouse)] = ttest(hits_late_all,miss_late_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:9
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - Misses')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Mid resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_mid(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_mid(2),3))], 'Color', av(2).col_str);
subplot(3,3,7)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Late resp')
subplot(3,3,8)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- avg')
subplot(3,3,9)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_late(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_late(2),3))], 'Color', av(2).col_str);
suptitle('Visual trials- Hits vs Misses- prev stim align- all cells')
print([fnout 'prevStim_align_SM_V.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    miss_trans_all = [];
    miss_mid_all = [];
    miss_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(4).av(2).outcome(5).trans_resp(1,:);
        miss_trans = mouse(imouse).expt(iexp).align(4).av(2).outcome(2).trans_resp(1,:);
        scatter(hits_trans, miss_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(miss_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(miss_trans,[],2)./sqrt(size(miss_trans,2)), {['o' col], col, col});
        hold on
        subplot(3,3,4)
        hits_mid = mouse(imouse).expt(iexp).align(4).av(2).outcome(5).mid_resp(1,:);
        miss_mid = mouse(imouse).expt(iexp).align(4).av(2).outcome(2).mid_resp(1,:);
        scatter(hits_mid, miss_mid, ['.' col]);
        hold on
        subplot(3,3,5)
        errorbarxy(mean(hits_mid,2), mean(miss_mid,2), std(hits_mid,[],2)./sqrt(size(hits_mid,2)), std(miss_mid,[],2)./sqrt(size(miss_mid,2)), {['o' col], col, col});
        hold on
        subplot(3,3,7)
        hits_late = mouse(imouse).expt(iexp).align(4).av(2).outcome(5).late_resp(1,:);
        miss_late = mouse(imouse).expt(iexp).align(4).av(2).outcome(2).late_resp(1,:);
        scatter(hits_late, miss_late, ['.' col]);
        hold on
        subplot(3,3,8)
        errorbarxy(mean(hits_late,2), mean(miss_late,2), std(hits_late,[],2)./sqrt(size(hits_late,2)), std(miss_late,[],2)./sqrt(size(miss_late,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        miss_trans_all = [miss_trans_all miss_trans];
        miss_mid_all = [miss_mid_all miss_mid];
        miss_late_all = [miss_late_all miss_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(miss_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(miss_trans_all,[],2)./sqrt(size(miss_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,miss_trans_all);
    hold on
    subplot(3,3,6)
    errorbarxy(mean(hits_mid_all,2), mean(miss_mid_all,2), std(hits_mid_all,[],2)./sqrt(size(hits_mid_all,2)), std(miss_mid_all,[],2)./sqrt(size(miss_mid_all,2)), {['o' col], col, col});
    [h_mid(imouse), p_mid(imouse)] = ttest(hits_mid_all,miss_mid_all);
    hold on
    subplot(3,3,9)
    errorbarxy(mean(hits_late_all,2), mean(miss_late_all,2), std(hits_late_all,[],2)./sqrt(size(hits_late_all,2)), std(miss_late_all,[],2)./sqrt(size(miss_late_all,2)), {['o' col], col, col});
    [h_late(imouse), p_late(imouse)] = ttest(hits_late_all,miss_late_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:9
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - Misses')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
legend(['n = ' num2str(n(1))], ['n = ' num2str(n(2))]);
subplot(3,3,4)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Mid resp')
subplot(3,3,5)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- avg')
subplot(3,3,6)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Mid resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_mid(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_mid(2),3))], 'Color', av(2).col_str);
subplot(3,3,7)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Late resp')
subplot(3,3,8)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- avg')
subplot(3,3,9)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Late resp- expt avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_late(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_late(2),3))], 'Color', av(2).col_str);
suptitle('Auditory trials- Hits vs Misses- prev stim align- all cells')
print([fnout 'prevStim_align_SM_A.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    FAs_trans_all = [];
    FAs_mid_all = [];
    FAs_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(5).av(1).outcome(1).trans_resp(1,:);
        FAs_trans = mouse(imouse).expt(iexp).align(5).av(1).outcome(3).trans_resp(1,:);
        scatter(hits_trans, FAs_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(FAs_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(FAs_trans,[],2)./sqrt(size(FAs_trans,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        FAs_trans_all = [FAs_trans_all FAs_trans];
        FAs_mid_all = [FAs_mid_all FAs_mid];
        FAs_late_all = [FAs_late_all FAs_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(FAs_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(FAs_trans_all,[],2)./sqrt(size(FAs_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,FAs_trans_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:3
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
suptitle('Visual trials- Hits vs FAs- prev base stim align- all cells')
print([fnout 'prevBaseStim_align_SE_V.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    FAs_trans_all = [];
    FAs_mid_all = [];
    FAs_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(5).av(2).outcome(1).trans_resp(1,:);
        FAs_trans = mouse(imouse).expt(iexp).align(5).av(2).outcome(3).trans_resp(1,:);
        scatter(hits_trans, FAs_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(FAs_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(FAs_trans,[],2)./sqrt(size(FAs_trans,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        FAs_trans_all = [FAs_trans_all FAs_trans];
        FAs_mid_all = [FAs_mid_all FAs_mid];
        FAs_late_all = [FAs_late_all FAs_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(FAs_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(FAs_trans_all,[],2)./sqrt(size(FAs_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,FAs_trans_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:3
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - FAs')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
suptitle('Auditory trials- Hits vs FAs- prev base stim align- all cells')
print([fnout 'prevBaseStim_align_SE_A.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    miss_trans_all = [];
    miss_mid_all = [];
    miss_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(5).av(1).outcome(5).trans_resp(1,:);
        miss_trans = mouse(imouse).expt(iexp).align(5).av(1).outcome(2).trans_resp(1,:);
        scatter(hits_trans, miss_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(miss_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(miss_trans,[],2)./sqrt(size(miss_trans,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        miss_trans_all = [miss_trans_all miss_trans];
        miss_mid_all = [miss_mid_all miss_mid];
        miss_late_all = [miss_late_all miss_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(miss_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(miss_trans_all,[],2)./sqrt(size(miss_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,miss_trans_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:3
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - Misses')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
suptitle('Visual trials- Hits vs Misses- prev base stim align- all cells')
print([fnout 'prevBaseStim_align_SM_V.pdf'], '-dpdf')

figure;
x = [-.5:.01:.5];
y = x;
for imouse = 1:size({av.mouse},2)
    col = av(imouse).col_str;
    hits_trans_all = [];
    hits_mid_all = [];
    hits_late_all = [];
    miss_trans_all = [];
    miss_mid_all = [];
    miss_late_all = [];
    for iexp = 1:size(mouse(imouse).expt,2)
        subplot(3,3,1)
        hits_trans = mouse(imouse).expt(iexp).align(5).av(2).outcome(5).trans_resp(1,:);
        miss_trans = mouse(imouse).expt(iexp).align(5).av(2).outcome(2).trans_resp(1,:);
        scatter(hits_trans, miss_trans, ['.' col]);
        hold on
        subplot(3,3,2)
        errorbarxy(mean(hits_trans,2), mean(miss_trans,2), std(hits_trans,[],2)./sqrt(size(hits_trans,2)), std(miss_trans,[],2)./sqrt(size(miss_trans,2)), {['o' col], col, col});
        hold on
        hits_trans_all = [hits_trans_all hits_trans];
        hits_mid_all = [hits_mid_all hits_mid];
        hits_late_all = [hits_late_all hits_late];
        miss_trans_all = [miss_trans_all miss_trans];
        miss_mid_all = [miss_mid_all miss_mid];
        miss_late_all = [miss_late_all miss_late];
    end
    subplot(3,3,3)
    errorbarxy(mean(hits_trans_all,2), mean(miss_trans_all,2), std(hits_trans_all,[],2)./sqrt(size(hits_trans_all,2)), std(miss_trans_all,[],2)./sqrt(size(miss_trans_all,2)), {['o' col], col, col});
    [h_trans(imouse), p_trans(imouse)] = ttest(hits_trans_all,miss_trans_all);
    hold on
    n(imouse) = length(hits_trans_all);
end
for i = 1:3
subplot(3,3,i)
plot(x,y,'--k')
hold on
hline(0)
vline(0)
xlabel('dF/F - Hits')
ylabel('dF/F - Misses')
end
subplot(3,3,1)
xlim([-0.25 0.25])
ylim([-0.25 0.25])
title('Transient resp')
subplot(3,3,2)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- expt avg')
subplot(3,3,3)
xlim([-.01 0.03])
ylim([-.01 0.03])
title('Transient resp- avg')
text(-0.005, 0.025, ['p = ' num2str(chop(p_trans(1),3))], 'Color', av(1).col_str);
hold on
text(-0.005, 0.02, ['p = ' num2str(chop(p_trans(2),3))], 'Color', av(2).col_str);
suptitle('Auditory trials- Hits vs Misses- prev base stim align- all cells')
print([fnout 'prevBaseStim_align_SM_A.pdf'], '-dpdf')

for imouse = 1:size({av.mouse},2)
    for iOri = 2:size(mouse(imouse).expt(1).cells,2)
        for i = 4
            if size(mouse(imouse).expt,2)>0
                mouse(imouse).align(i).name = mouse(imouse).expt(1).align(i).name;
                mouse(imouse).cells(iOri).name = mouse(imouse).expt(1).cells(iOri).name;
            end
            for ii = 1:2
                for iii = [1 2 3 5]
                    mouse(imouse).align(i).av(ii).outcome(iii).resp = [];
                    mouse(imouse).align(i).av(ii).outcome(iii).n = 0;
                    for iexp = 1:size(mouse(imouse).expt,2)
                        ind = mouse(imouse).expt(iexp).cells(iOri).ind;
                        mouse(imouse).align(i).av(ii).outcome(iii).resp = cat(2, mouse(imouse).align(i).av(ii).outcome(iii).resp, mean(mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).resp(:,ind,:),3));
                        mouse(imouse).align(i).av(ii).outcome(iii).n = mouse(imouse).align(i).av(ii).outcome(iii).n + mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).n;
                    end
                end
            end
            nCells = size(mouse(imouse).align(i).av(1).outcome(1).resp,2);
            figure;
            tt = [-pre_event_frames:post_event_frames-1].*(1000/frame_rate);
            S1 = mean(mouse(imouse).align(i).av(1).outcome(1).resp,2);
            F1 = mean(mouse(imouse).align(i).av(1).outcome(3).resp,2);
            S2 = mean(mouse(imouse).align(i).av(2).outcome(1).resp,2);
            F2 = mean(mouse(imouse).align(i).av(2).outcome(3).resp,2);
            S1M = mean(mouse(imouse).align(i).av(1).outcome(5).resp,2);
            M1M = mean(mouse(imouse).align(i).av(1).outcome(2).resp,2);
            S2M = mean(mouse(imouse).align(i).av(2).outcome(5).resp,2);
            M2M = mean(mouse(imouse).align(i).av(2).outcome(2).resp,2);
            S1s = std(mouse(imouse).align(i).av(1).outcome(1).resp,[],2)./sqrt(nCells);
            F1s = std(mouse(imouse).align(i).av(1).outcome(3).resp,[],2)./sqrt(nCells);
            S2s = std(mouse(imouse).align(i).av(2).outcome(1).resp,[],2)./sqrt(nCells);
            F2s = std(mouse(imouse).align(i).av(2).outcome(3).resp,[],2)./sqrt(nCells);
            S1Ms = std(mouse(imouse).align(i).av(1).outcome(5).resp,[],2)./sqrt(nCells);
            M1Ms = std(mouse(imouse).align(i).av(1).outcome(2).resp,[],2)./sqrt(nCells);
            S2Ms = std(mouse(imouse).align(i).av(2).outcome(5).resp,[],2)./sqrt(nCells);
            M2Ms = std(mouse(imouse).align(i).av(2).outcome(2).resp,[],2)./sqrt(nCells);
            subplot(2,2,1)
            shadedErrorBar(tt,S1,S1s, 'k');
            hold on
            shadedErrorBar(tt,F1,F1s, 'r');
            hold on
            vline(0)
            title(['Visual: Hits (' num2str(mouse(imouse).align(i).av(1).outcome(1).n) ') vs FAs (' num2str(mouse(imouse).align(i).av(1).outcome(3).n) ')'])
            subplot(2,2,2)
            shadedErrorBar(tt,S2,S2s, 'k');
            hold on
            shadedErrorBar(tt,F2,F2s, 'r');
            hold on
            vline(0)
            title(['Auditory: Hits (' num2str(mouse(imouse).align(i).av(2).outcome(1).n) ') vs FAs (' num2str(mouse(imouse).align(i).av(2).outcome(3).n) ')'])
            subplot(2,2,3)
            shadedErrorBar(tt,S1M,S1Ms, 'k');
            hold on
            shadedErrorBar(tt,M1M,M1Ms, 'r');
            hold on
            vline(0)
            title(['Visual: Hits (' num2str(mouse(imouse).align(i).av(1).outcome(5).n) ') vs Misses (' num2str(mouse(imouse).align(i).av(1).outcome(2).n) ')'])
            subplot(2,2,4)
            shadedErrorBar(tt,S2M,S2Ms, 'k');
            hold on
            shadedErrorBar(tt,M2M,M2Ms, 'r');
            hold on
            vline(0)
            title(['Auditory: Hits (' num2str(mouse(imouse).align(i).av(2).outcome(5).n) ') vs Misses (' num2str(mouse(imouse).align(i).av(2).outcome(2).n) ')'])
            suptitle(['i' num2str(av(imouse).mouse) '- Align to ' mouse(imouse).align(i).name '; ' mouse(imouse).cells(iOri).name ' pref cells: n = ' num2str(nCells) ' cells'])
            print([fnout mouse(imouse).align(i).name '_align_TC_' mouse(imouse).cells(iOri).name 'prefcells.pdf'], '-dpdf')
        end
    end
end

imouse = 2;
iexp = 1;
base_dirs = zeros(1,size(mouse(imouse).expt(iexp).target,2));
for iDir = 1:size(mouse(imouse).expt(iexp).target,2)
    base_dirs(1,iDir) = chop(mouse(imouse).expt(iexp).target(iDir).name,2);
end

i = 4;
ii = 1;
iii = 1;
for imouse = 1:size(av,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        figure;
        nDir = size(mouse(imouse).expt(iexp).target,2);
        Dirs = zeros(1,nDir);
        dir_ind = [];
        for iDir = 1:nDir
            Dirs(1,iDir) = chop(mouse(imouse).expt(iexp).target(iDir).name,2);
            if find(base_dirs == Dirs(1,iDir))
                dir_ind = [dir_ind iDir];
            end
        end
        for iOri = 1:4
            ind = mouse(imouse).expt(iexp).cells(iOri+1).ind;
            resp_mat = zeros(nDir,3);
            for iDir = 1:nDir
                resp_avg = squeeze(mean(mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).target(iDir).trans_resp(1,ind),2));
                resp_sem = squeeze(std(mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).target(iDir).trans_resp(1,ind),[],2)./sqrt(length(ind)));
                resp_mat(iDir,:) = [mouse(imouse).expt(iexp).target(iDir).name resp_avg resp_sem];
            end
            for iDir = 1:length(base_dirs)
                targ = dir_ind(iDir);
                if iexp == 1
                    mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp = mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).target(targ).trans_resp(1,ind);
                else
                    mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp = [mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp mouse(imouse).expt(iexp).align(i).av(ii).outcome(iii).target(targ).trans_resp(1,ind)];
                end
            end
            subplot(2,2,iOri)
            errorbar(resp_mat(:,1), resp_mat(:,2), resp_mat(:,3),'-ok');
            title([num2str(mouse(imouse).expt(iexp).cells(iOri+1).name) ' deg pref'])
            xlabel('Target Orientation (deg)')
            ylabel('dF/F')
            ylim([-0.02 0.03])
        end
        suptitle(['i' num2str(av(imouse).mouse)  ' Expt #' mouse(imouse).expt(iexp).date])
        fn = fullfile(rc.caOutputDir, mouse_name, date_name, [date_name '_' mouse_name '_']);
        print([fn 'prevstim_align_OriTuning.pdf'], '-dpdf')
    end
    figure;
    for iOri = 1:4
        subplot(2,2,iOri)
        resp_mat = zeros(size(base_dirs,2),2);
        for iDir = 1:size(base_dirs,2)
            resp_avg = squeeze(mean(mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp,2));
            sz = size(mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp,2);
            resp_sem = squeeze(std(mouse(imouse).align(i).av(ii).outcome(iii).target(iDir).cells(iOri+1).trans_resp,[],2)./sqrt(sz));
            resp_mat(iDir,:) = [resp_avg resp_sem];
        end
        errorbar(base_dirs, resp_mat(:,1), resp_mat(:,2),'-ok');
        title([num2str(mouse(imouse).expt(iexp).cells(iOri+1).name) ' deg pref'])
        xlabel('Target Orientation (deg)')
        ylabel('dF/F')
        ylim([-0.01 0.02])
    end
    alignYaxes
    suptitle(['i' num2str(av(imouse).mouse)])
    fn = fullfile(rc.caOutputDir, [date '_i' num2str(av(imouse).mouse)]);
    print([fn '_prevstim_align_OriTuning.pdf'], '-dpdf')
end
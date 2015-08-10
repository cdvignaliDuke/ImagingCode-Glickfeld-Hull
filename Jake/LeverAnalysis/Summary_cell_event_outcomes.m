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
    load([dest_sub '_event_outcome_summary.mat'])
    load([dest_sub '_nextTrOutcome.mat'])
    load([dest_sub '_cell_categories.mat']);
    
    ncells(id) = size(press,2);
    
    pctCorrectNext_release{id} = binofit(sum(successIx,2),length(successIx));
    pctCorrectNext_success{id} = binofit(sum(nextTrOutcome_success,2),length(nextTrOutcome_success));
    pctCorrectNext_fail{id} = binofit(sum(nextTrOutcome_fail,2),length(nextTrOutcome_fail));
    pctCorrectNext_press{id} = binofit(sum(nextTrOutcome_press,2),length(nextTrOutcome_press));
    
    pctCorrectSecond_success{id} = binofit(sum(secondTrOutcome_success,2),length(secondTrOutcome_success));
    pctCorrectSecond_fail{id} = binofit(sum(secondTrOutcome_fail,2),length(secondTrOutcome_fail));
    pctCorrectSecond_press{id} = binofit(sum(secondTrOutcome_press,2),length(secondTrOutcome_press));
    
    pctCorrectPrev_success{id} = binofit(sum(prevTrOutcome_success,2),length(prevTrOutcome_success));
    pctCorrectPrev_fail{id} = binofit(sum(prevTrOutcome_fail,2),length(prevTrOutcome_fail));
    pctCorrectPrev_press{id} = binofit(sum(prevTrOutcome_press,2),length(prevTrOutcome_press));
    
    pctCorrectNext_release_CS_temp = [];
    pctCorrectNext_release_noCS_temp = [];
    pctCorrectNext_success_CS_temp = [];
    pctCorrectNext_success_noCS_temp = [];
    pctCorrectNext_fail_CS_temp = [];
    pctCorrectNext_fail_noCS_temp = [];
    pctCorrectNext_press_CS_temp = [];
    pctCorrectNext_press_noCS_temp = [];
    for ic = 1:ncells(id)
        pctCorrectNext_release_CS_temp(ic,:) = [release(ic).CS(1).pctCorrectNext release(ic).CS(3).pctCorrectNext];
        pctCorrectNext_release_noCS_temp(ic,:) = [release(ic).CS(2).pctCorrectNext release(ic).CS(4).pctCorrectNext];
        pctCorrectNext_success_CS_temp(ic,:) = [success(ic).CS(1).pctCorrectNext success(ic).CS(3).pctCorrectNext];
        pctCorrectNext_success_noCS_temp(ic,:) = [success(ic).CS(2).pctCorrectNext success(ic).CS(4).pctCorrectNext];
        pctCorrectNext_fail_CS_temp(ic,:) = [fail(ic).CS(1).pctCorrectNext fail(ic).CS(3).pctCorrectNext];
        pctCorrectNext_fail_noCS_temp(ic,:) = [fail(ic).CS(2).pctCorrectNext fail(ic).CS(4).pctCorrectNext];
        pctCorrectNext_press_CS_temp(ic,:) = [press(ic).CS(1).pctCorrectNext press(ic).CS(3).pctCorrectNext];
        pctCorrectNext_press_noCS_temp(ic,:) = [press(ic).CS(2).pctCorrectNext press(ic).CS(4).pctCorrectNext];
    end
    pctCorrectNext_release_CS{id} = pctCorrectNext_release_CS_temp;
    pctCorrectNext_release_noCS{id} = pctCorrectNext_release_noCS_temp;
    pctCorrectNext_success_CS{id} = pctCorrectNext_success_CS_temp;
    pctCorrectNext_success_noCS{id} = pctCorrectNext_success_noCS_temp;
    pctCorrectNext_fail_CS{id} = pctCorrectNext_fail_CS_temp;
    pctCorrectNext_fail_noCS{id} = pctCorrectNext_fail_noCS_temp;
    pctCorrectNext_press_CS{id} = pctCorrectNext_press_CS_temp;
    pctCorrectNext_press_noCS{id} = pctCorrectNext_press_noCS_temp;
    
    pctCorrectThis_press_CS_temp = [];
    pctCorrectThis_press_noCS_temp = [];
    for ic = 1:ncells(id)
        pctCorrectThis_press_CS_temp(ic,:) = [press(ic).CS(1).pctCorrectThis press(ic).CS(3).pctCorrectThis];
        pctCorrectThis_press_noCS_temp(ic,:) = [press(ic).CS(2).pctCorrectThis press(ic).CS(4).pctCorrectThis];
    end
    pctCorrectThis_press_CS{id} = pctCorrectThis_press_CS_temp;
    pctCorrectThis_press_noCS{id} = pctCorrectThis_press_noCS_temp;
    
    pctCorrectRand_release_CS_temp = [];
    pctCorrectRand_release_noCS_temp = [];
    pctCorrectRand_success_CS_temp = [];
    pctCorrectRand_success_noCS_temp = [];
    pctCorrectRand_fail_CS_temp = [];
    pctCorrectRand_fail_noCS_temp = [];
    pctCorrectRand_press_CS_temp = [];
    pctCorrectRand_press_noCS_temp = [];
    for ic = 1:ncells(id)
        pctCorrectRand_release_CS_temp(ic,:) = [release(ic).CS(1).pctCorrectRand release(ic).CS(3).pctCorrectRand];
        pctCorrectRand_release_noCS_temp(ic,:) = [release(ic).CS(2).pctCorrectRand release(ic).CS(4).pctCorrectRand];
        pctCorrectRand_success_CS_temp(ic,:) = [success(ic).CS(1).pctCorrectRand success(ic).CS(3).pctCorrectRand];
        pctCorrectRand_success_noCS_temp(ic,:) = [success(ic).CS(2).pctCorrectRand success(ic).CS(4).pctCorrectRand];
        pctCorrectRand_fail_CS_temp(ic,:) = [fail(ic).CS(1).pctCorrectRand fail(ic).CS(3).pctCorrectRand];
        pctCorrectRand_fail_noCS_temp(ic,:) = [fail(ic).CS(2).pctCorrectRand fail(ic).CS(4).pctCorrectRand];
        pctCorrectRand_press_CS_temp(ic,:) = [press(ic).CS(1).pctCorrectRand press(ic).CS(3).pctCorrectRand];
        pctCorrectRand_press_noCS_temp(ic,:) = [press(ic).CS(2).pctCorrectRand press(ic).CS(4).pctCorrectRand];
    end
    pctCorrectRand_release_CS{id} = pctCorrectRand_release_CS_temp;
    pctCorrectRand_release_noCS{id} = pctCorrectRand_release_noCS_temp;
    pctCorrectRand_success_CS{id} = pctCorrectRand_success_CS_temp;
    pctCorrectRand_success_noCS{id} = pctCorrectRand_success_noCS_temp;
    pctCorrectRand_fail_CS{id} = pctCorrectRand_fail_CS_temp;
    pctCorrectRand_fail_noCS{id} = pctCorrectRand_fail_noCS_temp;
    pctCorrectRand_press_CS{id} = pctCorrectRand_press_CS_temp;
    pctCorrectRand_press_noCS{id} = pctCorrectRand_press_noCS_temp;
    
    pctCorrectSecond_release_CS_temp = [];
    pctCorrectSecond_release_noCS_temp = [];
    pctCorrectSecond_success_CS_temp = [];
    pctCorrectSecond_success_noCS_temp = [];
    pctCorrectSecond_fail_CS_temp = [];
    pctCorrectSecond_fail_noCS_temp = [];
    pctCorrectSecond_press_CS_temp = [];
    pctCorrectSecond_press_noCS_temp = [];
    for ic = 1:ncells(id)
        pctCorrectSecond_release_CS_temp(ic,:) = [release(ic).CS(1).pctCorrectSecond release(ic).CS(3).pctCorrectSecond];
        pctCorrectSecond_release_noCS_temp(ic,:) = [release(ic).CS(2).pctCorrectSecond release(ic).CS(4).pctCorrectSecond];
        pctCorrectSecond_success_CS_temp(ic,:) = [success(ic).CS(1).pctCorrectSecond success(ic).CS(3).pctCorrectSecond];
        pctCorrectSecond_success_noCS_temp(ic,:) = [success(ic).CS(2).pctCorrectSecond success(ic).CS(4).pctCorrectSecond];
        pctCorrectSecond_fail_CS_temp(ic,:) = [fail(ic).CS(1).pctCorrectSecond fail(ic).CS(3).pctCorrectSecond];
        pctCorrectSecond_fail_noCS_temp(ic,:) = [fail(ic).CS(2).pctCorrectSecond fail(ic).CS(4).pctCorrectSecond];
        pctCorrectSecond_press_CS_temp(ic,:) = [press(ic).CS(1).pctCorrectSecond press(ic).CS(3).pctCorrectSecond];
        pctCorrectSecond_press_noCS_temp(ic,:) = [press(ic).CS(2).pctCorrectSecond press(ic).CS(4).pctCorrectSecond];
    end
    pctCorrectSecond_release_CS{id} = pctCorrectSecond_release_CS_temp;
    pctCorrectSecond_release_noCS{id} = pctCorrectSecond_release_noCS_temp;
    pctCorrectSecond_success_CS{id} = pctCorrectSecond_success_CS_temp;
    pctCorrectSecond_success_noCS{id} = pctCorrectSecond_success_noCS_temp;
    pctCorrectSecond_fail_CS{id} = pctCorrectSecond_fail_CS_temp;
    pctCorrectSecond_fail_noCS{id} = pctCorrectSecond_fail_noCS_temp;
    pctCorrectSecond_press_CS{id} = pctCorrectSecond_press_CS_temp;
    pctCorrectSecond_press_noCS{id} = pctCorrectSecond_press_noCS_temp;
    
    pctCorrectPrev_release_CS_temp = [];
    pctCorrectPrev_release_noCS_temp = [];
    pctCorrectPrev_success_CS_temp = [];
    pctCorrectPrev_success_noCS_temp = [];
    pctCorrectPrev_fail_CS_temp = [];
    pctCorrectPrev_fail_noCS_temp = [];
    pctCorrectPrev_press_CS_temp = [];
    pctCorrectPrev_press_noCS_temp = [];
    for ic = 1:ncells(id)
        pctCorrectPrev_release_CS_temp(ic,:) = [release(ic).CS(1).pctCorrectPrev release(ic).CS(3).pctCorrectPrev];
        pctCorrectPrev_release_noCS_temp(ic,:) = [release(ic).CS(2).pctCorrectPrev release(ic).CS(4).pctCorrectPrev];
        pctCorrectPrev_success_CS_temp(ic,:) = [success(ic).CS(1).pctCorrectPrev success(ic).CS(3).pctCorrectPrev];
        pctCorrectPrev_success_noCS_temp(ic,:) = [success(ic).CS(2).pctCorrectPrev success(ic).CS(4).pctCorrectPrev];
        pctCorrectPrev_fail_CS_temp(ic,:) = [fail(ic).CS(1).pctCorrectPrev fail(ic).CS(3).pctCorrectPrev];
        pctCorrectPrev_fail_noCS_temp(ic,:) = [fail(ic).CS(2).pctCorrectPrev fail(ic).CS(4).pctCorrectPrev];
        pctCorrectPrev_press_CS_temp(ic,:) = [press(ic).CS(1).pctCorrectPrev press(ic).CS(3).pctCorrectPrev];
        pctCorrectPrev_press_noCS_temp(ic,:) = [press(ic).CS(2).pctCorrectPrev press(ic).CS(4).pctCorrectPrev];
    end
    pctCorrectPrev_release_CS{id} = pctCorrectPrev_release_CS_temp;
    pctCorrectPrev_release_noCS{id} = pctCorrectPrev_release_noCS_temp;
    pctCorrectPrev_success_CS{id} = pctCorrectPrev_success_CS_temp;
    pctCorrectPrev_success_noCS{id} = pctCorrectPrev_success_noCS_temp;
    pctCorrectPrev_fail_CS{id} = pctCorrectPrev_fail_CS_temp;
    pctCorrectPrev_fail_noCS{id} = pctCorrectPrev_fail_noCS_temp;
    pctCorrectPrev_press_CS{id} = pctCorrectPrev_press_CS_temp;
    pctCorrectPrev_press_noCS{id} = pctCorrectPrev_press_noCS_temp;
end

%% plotting
col_mat = strvcat('r', 'b', 'r', 'b', 'g', 'm');

%Figure 1 -percent correct after a press CS or no CS each trial
figure;
allPrev = [];
allThis = [];
allNext = [];
allSecond = [];
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    plot(pctCorrectPrev_press_CS{id}(:,1), pctCorrectPrev_press_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allPrev = [allPrev; pctCorrectPrev_press_CS{id}(:,1) pctCorrectPrev_press_noCS{id}(:,1)];
    hold on
    subplot(2,2,2)
    plot(pctCorrectThis_press_CS{id}(:,1), pctCorrectThis_press_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allThis = [allThis; pctCorrectThis_press_CS{id}(:,1) pctCorrectThis_press_noCS{id}(:,1)];
    hold on
    subplot(2,2,3)
    plot(pctCorrectNext_press_CS{id}(:,1), pctCorrectNext_press_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allNext = [allNext; pctCorrectNext_press_CS{id}(:,1) pctCorrectNext_press_noCS{id}(:,1)];
    hold on
    subplot(2,2,4)
    plot(pctCorrectSecond_press_CS{id}(:,1), pctCorrectSecond_press_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allSecond = [allSecond; pctCorrectSecond_press_CS{id}(:,1) pctCorrectSecond_press_noCS{id}(:,1)];
    hold on
end
x = 0:.01:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allPrev(:,1),allPrev(:,2));
title(['Prev- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,2)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allThis(:,1),allThis(:,2));
title(['This- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,3)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allNext(:,1),allNext(:,2));
title(['Next- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,4)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allSecond(:,1),allSecond(:,2));
title(['Second- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])


suptitle('Prediction of outcome by trial order by press CS')
print([out_base 'Summary_pctcorr_pressCSnoCS_byoutcome.eps'], '-depsc');
print([out_base 'Summary_pctcorr_pressCSnoCS_byoutcome.pdf'], '-dpdf');


%Figure 1 -percent correct after a press CS or no CS each trial
figure;
allPrev = [];
allThis = [];
allNext = [];
allSecond = [];
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    plot(pctCorrectPrev_success_CS{id}(:,1), pctCorrectPrev_success_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allPrev = [allPrev; pctCorrectPrev_success_CS{id}(:,1) pctCorrectPrev_success_noCS{id}(:,1)];
    hold on
    subplot(2,2,2)
    plot(pctCorrectNext_success_CS{id}(:,1), pctCorrectNext_success_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allNext = [allNext; pctCorrectNext_success_CS{id}(:,1) pctCorrectNext_success_noCS{id}(:,1)];
    hold on
    subplot(2,2,3)
    plot(pctCorrectSecond_success_CS{id}(:,1), pctCorrectSecond_success_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allSecond = [allSecond; pctCorrectSecond_success_CS{id}(:,1) pctCorrectSecond_success_noCS{id}(:,1)];
    hold on
end
x = 0:.01:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allPrev(:,1),allPrev(:,2));
title(['Prev- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,2)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allNext(:,1),allNext(:,2));
title(['Next- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,3)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allSecond(:,1),allSecond(:,2));
title(['Second- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])

suptitle('Prediction of outcome by trial order by success CS')
print([out_base 'Summary_pctcorr_successCSnoCS_byoutcome.eps'], '-depsc');
print([out_base 'Summary_pctcorr_successCSnoCS_byoutcome.pdf'], '-dpdf');

%Figure 1 -percent correct after a press CS or no CS each trial
figure;
allPrev = [];
allThis = [];
allNext = [];
allSecond = [];
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    plot(pctCorrectPrev_fail_CS{id}(:,1), pctCorrectPrev_fail_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allPrev = [allPrev; pctCorrectPrev_fail_CS{id}(:,1) pctCorrectPrev_fail_noCS{id}(:,1)];
    hold on
    subplot(2,2,2)
    plot(pctCorrectNext_fail_CS{id}(:,1), pctCorrectNext_fail_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allNext = [allNext; pctCorrectNext_fail_CS{id}(:,1) pctCorrectNext_fail_noCS{id}(:,1)];
    hold on
    subplot(2,2,3)
    plot(pctCorrectSecond_fail_CS{id}(:,1), pctCorrectSecond_fail_noCS{id}(:,1), ['o' col_mat(id,:)]);
    allSecond = [allSecond; pctCorrectSecond_fail_CS{id}(:,1) pctCorrectSecond_fail_noCS{id}(:,1)];
    hold on
end
x = 0:.01:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allPrev(:,1),allPrev(:,2));
title(['Prev- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,2)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allNext(:,1),allNext(:,2));
title(['Next- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])
subplot(2,2,3)
plot(x,y,'-k')
hold on
[h, p ] = ttest(allSecond(:,1),allSecond(:,2));
title(['Second- p = ' num2str(chop(p,2))])
xlabel('CS')
ylabel('no CS')
xlim([0 1])
ylim([0 1])

suptitle('Prediction of outcome by trial order by fail CS')
print([out_base 'Summary_pctcorr_failCSnoCS_byoutcome.eps'], '-depsc');
print([out_base 'Summary_pctcorr_failCSnoCS_byoutcome.pdf'], '-dpdf');

%Figure 2 -percent correct after a press CS on Prev, this, next, second
%trials
figure;
for id = 1:size(date_mat,1)
    subplot(2,3,1)
    plot(pctCorrectPrev_press_CS{id}(:,1), pctCorrectNext_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,3,2)
    plot(pctCorrectPrev_press_CS{id}(:,1), pctCorrectSecond_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,3,3)
    plot(pctCorrectNext_press_CS{id}(:,1), pctCorrectSecond_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,3,4)
    plot(pctCorrectThis_press_CS{id}(:,1), pctCorrectPrev_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,3,5)
    plot(pctCorrectThis_press_CS{id}(:,1), pctCorrectNext_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,3,6)
    plot(pctCorrectThis_press_CS{id}(:,1), pctCorrectSecond_press_CS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
end
x = 0:.01:1;
y = x;
subplot(2,3,1)
plot(x,y,'-k')
hold on
xlabel('CS- Prev')
ylabel('CS- Next')
xlim([0 1])
ylim([0 1])
subplot(2,3,2)
plot(x,y,'-k')
hold on
xlabel('CS- Prev')
ylabel('CS- Second')
xlim([0 1])
ylim([0 1])
subplot(2,3,3)
plot(x,y,'-k')
hold on
xlabel('CS- Next')
ylabel('CS- Second')
xlim([0 1])
ylim([0 1])
subplot(2,3,4)
plot(x,y,'-k')
hold on
xlabel('CS- This')
ylabel('CS- Prev')
xlim([0 1])
ylim([0 1])
subplot(2,3,5)
plot(x,y,'-k')
hold on
xlabel('CS- This')
ylabel('CS- Next')
xlim([0 1])
ylim([0 1])
subplot(2,3,6)
plot(x,y,'-k')
hold on
xlabel('CS- This')
ylabel('CS- Second')
xlim([0 1])
ylim([0 1])

suptitle('Prediction of outcome on this trial by press CS')
print([out_base 'Summary_probCS_press_byoutcome.eps'], '-depsc');
print([out_base 'Summary_probCS_press_byoutcome.pdf'], '-dpdf');

%Figure 2 -percent correct after a press CS or no CS THIS trial
figure;
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    plot(pctCorrectNext_success_CS{id}(:,1), pctCorrectSecond_success_noCS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,2,2)
    plot(pctCorrectNext_fail_CS{id}(:,1), pctCorrectSecond_fail_noCS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,2,3)
    plot(pctCorrectNext_success_CS{id}(:,1), pctCorrectPrev_success_noCS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
    subplot(2,2,4)
    plot(pctCorrectNext_fail_CS{id}(:,1), pctCorrectPrev_fail_noCS{id}(:,1), ['o' col_mat(id,:)]);
    hold on
end
x = 0:.01:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
title('Success- % correct')
xlabel('CS- Next')
ylabel('CS- Second')
xlim([0 1])
ylim([0 1])
subplot(2,2,2)
plot(x,y,'-k')
hold on
title('Fail- % correct')
xlabel('CS- Next')
ylabel('CS- Second')
xlim([0 1])
ylim([0 1])
subplot(2,2,3)
plot(x,y,'-k')
hold on
title('Success- % correct')
xlabel('CS- Next')
ylabel('CS- Prev')
xlim([0 1])
ylim([0 1])
subplot(2,2,4)
plot(x,y,'-k')
hold on
title('Fail- % correct')
xlabel('CS- Next')
ylabel('CS- Prev')
xlim([0 1])
ylim([0 1])

suptitle('Prediction of outcome on this trial by release CS')
print([out_base 'Summary_probCS_release_byoutcome.eps'], '-depsc');
print([out_base 'Summary_probCS_release_byoutcome.pdf'], '-dpdf');


%difference from average outcome for CS and no CS trials
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectNext_success_CS_sub = pctCorrectNext_success_CS{id}(:,1)-pctCorrectNext_success{id};
    pctCorrectNext_success_noCS_sub = pctCorrectNext_success_noCS{id}(:,1)-pctCorrectNext_success{id};
    scatter(pctCorrectNext_success_CS_sub, pctCorrectNext_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectNext_fail_CS_sub = pctCorrectNext_fail_CS{id}(:,1)-pctCorrectNext_fail{id};
    pctCorrectNext_fail_noCS_sub = pctCorrectNext_fail_noCS{id}(:,1)-pctCorrectNext_fail{id};
    scatter(pctCorrectNext_fail_CS_sub, pctCorrectNext_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectNext_release_CS_sub = pctCorrectNext_release_CS{id}(:,1)-pctCorrectNext_release{id};
    pctCorrectNext_release_noCS_sub = pctCorrectNext_release_noCS{id}(:,1)-pctCorrectNext_release{id};
    scatter(pctCorrectNext_release_CS_sub, pctCorrectNext_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectNext_press_CS_sub = pctCorrectNext_press_CS{id}(:,1)-pctCorrectNext_press{id};
    pctCorrectNext_press_noCS_sub = pctCorrectNext_press_noCS{id}(:,1)-pctCorrectNext_press{id};
    scatter(pctCorrectNext_press_CS_sub, pctCorrectNext_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
x = -1:0.1:1;
y = -x;
subplot(2,2,1)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on next trial'])
print([out_base 'Summary_diff_event_outcome.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome.pdf'], '-dpdf');


%difference from average outcome for CS and no CS trials- short window
figure
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectNext_success_CS_sub = pctCorrectNext_success_CS{id}(:,2)-pctCorrectNext_success{id};
    pctCorrectNext_success_noCS_sub = pctCorrectNext_success_noCS{id}(:,2)-pctCorrectNext_success{id};
    scatter(pctCorrectNext_success_CS_sub, pctCorrectNext_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectNext_fail_CS_sub = pctCorrectNext_fail_CS{id}(:,2)-pctCorrectNext_fail{id};
    pctCorrectNext_fail_noCS_sub = pctCorrectNext_fail_noCS{id}(:,2)-pctCorrectNext_fail{id};
    scatter(pctCorrectNext_fail_CS_sub, pctCorrectNext_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectNext_release_CS_sub = pctCorrectNext_release_CS{id}(:,2)-pctCorrectNext_release{id};
    pctCorrectNext_release_noCS_sub = pctCorrectNext_release_noCS{id}(:,2)-pctCorrectNext_release{id};
    scatter(pctCorrectNext_release_CS_sub, pctCorrectNext_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectNext_press_CS_sub = pctCorrectNext_press_CS{id}(:,2)-pctCorrectNext_press{id};
    pctCorrectNext_press_noCS_sub = pctCorrectNext_press_noCS{id}(:,2)-pctCorrectNext_press{id};
    scatter(pctCorrectNext_press_CS_sub, pctCorrectNext_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
x = -1:0.1:1;
y = -x;
subplot(2,2,1)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
plot(x,y,'-k')
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on next trial- 100 ms'])
print([out_base 'Summary_diff_event_outcome_100ms.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_100ms.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- random
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectRand_success_CS_sub = pctCorrectRand_success_CS{id}(:,1)-pctCorrectNext_success{id};
    pctCorrectRand_success_noCS_sub = pctCorrectRand_success_noCS{id}(:,1)-pctCorrectNext_success{id};
    scatter(pctCorrectRand_success_CS_sub, pctCorrectRand_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectRand_fail_CS_sub = pctCorrectRand_fail_CS{id}(:,1)-pctCorrectNext_fail{id};
    pctCorrectRand_fail_noCS_sub = pctCorrectRand_fail_noCS{id}(:,1)-pctCorrectNext_fail{id};
    scatter(pctCorrectRand_fail_CS_sub, pctCorrectRand_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectRand_release_CS_sub = pctCorrectRand_release_CS{id}(:,1)-pctCorrectNext_release{id};
    pctCorrectRand_release_noCS_sub = pctCorrectRand_release_noCS{id}(:,1)-pctCorrectNext_release{id};
    scatter(pctCorrectRand_release_CS_sub, pctCorrectRand_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectRand_press_CS_sub = pctCorrectRand_press_CS{id}(:,1)-pctCorrectNext_press{id};
    pctCorrectRand_press_noCS_sub = pctCorrectRand_press_noCS{id}(:,1)-pctCorrectNext_press{id};
    scatter(pctCorrectRand_press_CS_sub, pctCorrectRand_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on random trial'])
print([out_base 'Summary_diff_event_outcome_rand.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_rand.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- short window
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectRand_success_CS_sub = pctCorrectRand_success_CS{id}(:,2)-pctCorrectNext_success{id};
    pctCorrectRand_success_noCS_sub = pctCorrectRand_success_noCS{id}(:,2)-pctCorrectNext_success{id};
    scatter(pctCorrectRand_success_CS_sub, pctCorrectRand_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectRand_fail_CS_sub = pctCorrectRand_fail_CS{id}(:,2)-pctCorrectNext_fail{id};
    pctCorrectRand_fail_noCS_sub = pctCorrectRand_fail_noCS{id}(:,2)-pctCorrectNext_fail{id};
    scatter(pctCorrectRand_fail_CS_sub, pctCorrectRand_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectRand_release_CS_sub = pctCorrectRand_release_CS{id}(:,2)-pctCorrectNext_release{id};
    pctCorrectRand_release_noCS_sub = pctCorrectRand_release_noCS{id}(:,2)-pctCorrectNext_release{id};
    scatter(pctCorrectRand_release_CS_sub, pctCorrectRand_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectRand_press_CS_sub = pctCorrectRand_press_CS{id}(:,2)-pctCorrectNext_press{id};
    pctCorrectRand_press_noCS_sub = pctCorrectRand_press_noCS{id}(:,2)-pctCorrectNext_press{id};
    scatter(pctCorrectRand_press_CS_sub, pctCorrectRand_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
x = 0:.1:1;
y = x;
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on random trial- 100 ms'])
print([out_base 'Summary_diff_event_outcome_100ms_rand.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_100ms_rand.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- second
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectSecond_success_CS_sub = pctCorrectSecond_success_CS{id}(:,1)-pctCorrectSecond_success{id};
    pctCorrectSecond_success_noCS_sub = pctCorrectSecond_success_noCS{id}(:,1)-pctCorrectSecond_success{id};
    scatter(pctCorrectSecond_success_CS_sub, pctCorrectSecond_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectSecond_fail_CS_sub = pctCorrectSecond_fail_CS{id}(:,1)-pctCorrectSecond_fail{id};
    pctCorrectSecond_fail_noCS_sub = pctCorrectSecond_fail_noCS{id}(:,1)-pctCorrectSecond_fail{id};
    scatter(pctCorrectSecond_fail_CS_sub, pctCorrectSecond_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectSecond_release_CS_sub = pctCorrectSecond_release_CS{id}(:,1)-pctCorrectNext_release{id};
    pctCorrectSecond_release_noCS_sub = pctCorrectSecond_release_noCS{id}(:,1)-pctCorrectNext_release{id};
    scatter(pctCorrectSecond_release_CS_sub, pctCorrectSecond_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectSecond_press_CS_sub = pctCorrectSecond_press_CS{id}(:,1)-pctCorrectSecond_press{id};
    pctCorrectSecond_press_noCS_sub = pctCorrectSecond_press_noCS{id}(:,1)-pctCorrectSecond_press{id};
    scatter(pctCorrectSecond_press_CS_sub, pctCorrectSecond_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on second trial'])
print([out_base 'Summary_diff_event_outcome_second.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_second.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- second- 100 ms
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectSecond_success_CS_sub = pctCorrectSecond_success_CS{id}(:,2)-pctCorrectSecond_success{id};
    pctCorrectSecond_success_noCS_sub = pctCorrectSecond_success_noCS{id}(:,2)-pctCorrectSecond_success{id};
    scatter(pctCorrectSecond_success_CS_sub, pctCorrectSecond_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectSecond_fail_CS_sub = pctCorrectSecond_fail_CS{id}(:,2)-pctCorrectSecond_fail{id};
    pctCorrectSecond_fail_noCS_sub = pctCorrectSecond_fail_noCS{id}(:,2)-pctCorrectSecond_fail{id};
    scatter(pctCorrectSecond_fail_CS_sub, pctCorrectSecond_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectSecond_release_CS_sub = pctCorrectSecond_release_CS{id}(:,2)-pctCorrectNext_release{id};
    pctCorrectSecond_release_noCS_sub = pctCorrectSecond_release_noCS{id}(:,2)-pctCorrectNext_release{id};
    scatter(pctCorrectSecond_release_CS_sub, pctCorrectSecond_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectSecond_press_CS_sub = pctCorrectSecond_press_CS{id}(:,2)-pctCorrectSecond_press{id};
    pctCorrectSecond_press_noCS_sub = pctCorrectSecond_press_noCS{id}(:,2)-pctCorrectSecond_press{id};
    scatter(pctCorrectSecond_press_CS_sub, pctCorrectSecond_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on second trial- 100 ms'])
print([out_base 'Summary_diff_event_outcome_second_100ms.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_second_100ms.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- prev
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectPrev_success_CS_sub = pctCorrectPrev_success_CS{id}(:,1)-pctCorrectPrev_success{id};
    pctCorrectPrev_success_noCS_sub = pctCorrectPrev_success_noCS{id}(:,1)-pctCorrectPrev_success{id};
    scatter(pctCorrectPrev_success_CS_sub, pctCorrectPrev_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectPrev_fail_CS_sub = pctCorrectPrev_fail_CS{id}(:,1)-pctCorrectPrev_fail{id};
    pctCorrectPrev_fail_noCS_sub = pctCorrectPrev_fail_noCS{id}(:,1)-pctCorrectPrev_fail{id};
    scatter(pctCorrectPrev_fail_CS_sub, pctCorrectPrev_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectPrev_release_CS_sub = pctCorrectPrev_release_CS{id}(:,1)-pctCorrectNext_release{id};
    pctCorrectPrev_release_noCS_sub = pctCorrectPrev_release_noCS{id}(:,1)-pctCorrectNext_release{id};
    scatter(pctCorrectPrev_release_CS_sub, pctCorrectPrev_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectPrev_press_CS_sub = pctCorrectPrev_press_CS{id}(:,1)-pctCorrectPrev_press{id};
    pctCorrectPrev_press_noCS_sub = pctCorrectPrev_press_noCS{id}(:,1)-pctCorrectPrev_press{id};
    scatter(pctCorrectPrev_press_CS_sub, pctCorrectPrev_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on prev trial'])
print([out_base 'Summary_diff_event_outcome_prev.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_prev.pdf'], '-dpdf');

%difference from average outcome for CS and no CS trials- prev- 100 ms
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    pctCorrectPrev_success_CS_sub = pctCorrectPrev_success_CS{id}(:,2)-pctCorrectPrev_success{id};
    pctCorrectPrev_success_noCS_sub = pctCorrectPrev_success_noCS{id}(:,2)-pctCorrectPrev_success{id};
    scatter(pctCorrectPrev_success_CS_sub, pctCorrectPrev_success_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    pctCorrectPrev_fail_CS_sub = pctCorrectPrev_fail_CS{id}(:,2)-pctCorrectPrev_fail{id};
    pctCorrectPrev_fail_noCS_sub = pctCorrectPrev_fail_noCS{id}(:,2)-pctCorrectPrev_fail{id};
    scatter(pctCorrectPrev_fail_CS_sub, pctCorrectPrev_fail_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    pctCorrectPrev_release_CS_sub = pctCorrectPrev_release_CS{id}(:,2)-pctCorrectNext_release{id};
    pctCorrectPrev_release_noCS_sub = pctCorrectPrev_release_noCS{id}(:,2)-pctCorrectNext_release{id};
    scatter(pctCorrectPrev_release_CS_sub, pctCorrectPrev_release_noCS_sub, ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    pctCorrectPrev_press_CS_sub = pctCorrectPrev_press_CS{id}(:,2)-pctCorrectPrev_press{id};
    pctCorrectPrev_press_noCS_sub = pctCorrectPrev_press_noCS{id}(:,2)-pctCorrectPrev_press{id};
    scatter(pctCorrectPrev_press_CS_sub, pctCorrectPrev_press_noCS_sub, ['o' col_mat(id,:)])
    hold on
end
subplot(2,2,1)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Success')
subplot(2,2,2)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Fail')
subplot(2,2,3)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Release')
subplot(2,2,4)
xlim([-1 1])
ylim([-1 1])
xlabel('Diff- CS')
ylabel('Diff- no CS')
title('Press')
suptitle(['Difference from average %correct on prev trial- 100 ms'])
print([out_base 'Summary_diff_event_outcome_prev_100ms.eps'], '-depsc');
print([out_base 'Summary_diff_event_outcome_prev_100ms.pdf'], '-dpdf');



%Probability of a correct afer CS on Success/Fail
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    scatter(pctCorrectNext_success_CS{id}(:,1),pctCorrectNext_fail_CS{id}(:,1), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    scatter(pctCorrectNext_success_CS{id}(:,2),pctCorrectNext_fail_CS{id}(:,2), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    scatter(nanmean(pctCorrectNext_success_CS{id}(:,1),1),nanmean(pctCorrectNext_fail_CS{id}(:,1),1), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    scatter(nanmean(pctCorrectNext_success_CS{id}(:,2),1),nanmean(pctCorrectNext_fail_CS{id}(:,2),1), ['o' col_mat(id,:)])
    hold on
end
x = 0:.1:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS Success')
ylabel('CS Fail')
title('% correct - long')
subplot(2,2,2)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS Success')
ylabel('CS Fail')
title('% correct - short')
subplot(2,2,3)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS Success')
ylabel('CS Fail')
title('avg % correct - long')
subplot(2,2,4)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS Success')
ylabel('CS Fail')
title('avg % correct - short')
suptitle(['%correct on next trial after CS on Success/Fail'])
print([out_base 'Summary_pctcorrect_bySF.eps'], '-depsc');
print([out_base 'Summary_pctcorrect_bySF.pdf'], '-dpdf');

%Probability of a correct afer CS on Press/Release
figure
for id = 1:size(date_mat,1)
    subplot(2,2,1)
    scatter(pctCorrectNext_release_CS{id}(:,1),pctCorrectNext_press_CS{id}(:,1), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,2)
    scatter(pctCorrectNext_release_CS{id}(:,2),pctCorrectNext_press_CS{id}(:,2), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,3)
    scatter(nanmean(pctCorrectNext_release_CS{id}(:,1),1),nanmean(pctCorrectNext_press_CS{id}(:,1),1), ['o' col_mat(id,:)])
    hold on
    subplot(2,2,4)
    scatter(nanmean(pctCorrectNext_release_CS{id}(:,2),1),nanmean(pctCorrectNext_press_CS{id}(:,2),1), ['o' col_mat(id,:)])
    hold on
end
x = 0:.1:1;
y = x;
subplot(2,2,1)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS release')
ylabel('CS press')
title('% correct - long')
subplot(2,2,2)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS release')
ylabel('CS press')
title('% correct - short')
subplot(2,2,3)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS release')
ylabel('CS press')
title('avg % correct - long')
subplot(2,2,4)
plot(x,y,'-k')
hold on
xlim([0 1])
ylim([0 1])
xlabel('CS release')
ylabel('CS press')
title('avg % correct - short')
suptitle(['%correct on next trial after CS on release/press'])
print([out_base 'Summary_pctcorrect_byRP.eps'], '-depsc');
print([out_base 'Summary_pctcorrect_byRP.pdf'], '-dpdf');

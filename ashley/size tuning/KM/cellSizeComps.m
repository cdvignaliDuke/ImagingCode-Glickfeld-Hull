rc = behavConstsAV;
fout = fullfile(rc.ashleyAnalysis, 'SizeTuning');
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\cellBodyData_RF_sizeOnly.mat');
load('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning\cellBodyData.mat');
fout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning';
areas = ["V1","LM","AL","PM"];
nArea = length(areas);

%% RF size vs pref size
for i = 1:nArea
    temp_RF = eval(strcat('RFsize_', areas(i)));
    temp_pref = eval(strcat('prefSize_', areas(i)));
    temp_supp = eval(strcat('suppInd_', areas(i)));
    ism2 = find(temp_supp>0);
    figure(1)
    subplot(2,2,i)
    scatter(temp_RF(ism2,:), temp_pref(:,ism2)./2)
    hold on
    LM = fitlm(temp_RF(ism2,:), temp_pref(:,ism2)./2);
    x = [0:1:40];
    y = LM.Coefficients.Estimate(1)+(LM.Coefficients.Estimate(2).*x);
    plot(x,y);
%     b1 = temp_RF\(temp_pref'./2);
%     YCalc = temp_RF.*b1;
%     plot(temp_RF,YCalc);
%     scatter(temp_RF(ism2,:), temp_pref(:,ism2))
    ylim([0 40])
    ylabel('Pref Size')
    xlim([0 40])
    xlabel('RF Size')
    text(30,30,['p= ' num2str(chop(LM.Coefficients.pValue(2),2))])
    text(30,25,['slope= ' num2str(chop(LM.Coefficients.Estimate(2),2))])
    title(areas(i))
    figure(2)
    subplot(2,2,i)
    [n bin] = histc(temp_RF,[0:5:40]);
    for ii = 1:length(n)
        if n(ii)>5
            ind = intersect(ism2,find(bin==ii));
            errorbarxy(mean(temp_RF(ind,:),1), mean(temp_pref(:,ind)./2,2), std(temp_RF(ind,:),[],1), std(temp_pref(:,ind)./2,[],2));
            hold on
        end
    end
    xlim([0 30])
    xlabel('RF Size')
    ylim([0 30])
    ylabel('Pref Size')
    title(areas(i))
end
figure(1)
print([fout '/RFsizePrefComp_scatter.pdf'],'-dpdf','-bestfit')
figure(2)
print([fout '/RFsizePrefComp_bins.pdf'],'-dpdf','-bestfit')

%% Supp ind for suppressed and all cells

figure;
supp_all = [];
supp_m2 = [];
for i = 1:nArea
    temp_supp = eval(strcat('suppInd_', areas(i)));
    supp_all = [supp_all; temp_supp' i.*ones(size(temp_supp'))];
    ism2 = find(temp_supp>0);
    supp_m2 = [supp_m2; temp_supp(:,ism2)' i.*ones(size(temp_supp(:,ism2)'))];
    subplot(1,2,1)
    errorbar(i,mean(temp_supp,2), std(temp_supp,[],2),'o');
    hold on
    subplot(1,2,2)
    errorbar(i,mean(temp_supp(:,ism2),2), std(temp_supp(:,ism2),[],2),'o');
    hold on
end
subplot(1,2,1)
[p_all tab_all stats_all] = anova1(supp_all(:,1), supp_all(:,2),'off');
title(['All cells- p = ' num2str(chop(p_all,2))])
set(gca, 'XTick',[1:4], 'XTickLabel', areas,'xlim',[0 nArea+1],'ylim', [0 1])
axis square
c_all = multcompare(stats_all,'display','off');
text(0.2, 0.2, num2str(chop(c_all(:,[1 2 6]),2)))
subplot(1,2,2)
[p_m2 tab_m2 stats_m2] = anova1(supp_m2(:,1), supp_m2(:,2),'off');
title(['Suppressed cells- p = ' num2str(chop(p_m2,2))])
set(gca, 'XTick',[1:4], 'XTickLabel', areas,'xlim',[0 nArea+1],'ylim', [0 1])
c_m2 = multcompare(stats_m2,'display','off');
text(0.2, 0.2, num2str(chop(c_m2(:,[1 2 6]),2)))
axis square
print([fout '\SuppInd_CompAllvsM2.pdf'],'-dpdf','-bestfit')

figure;

for i = 1:nArea
    temp_supp = eval(strcat('suppInd_', areas(i)));
    subplot(2,2,i)
    hist(temp_supp)
    title(areas(i))
end

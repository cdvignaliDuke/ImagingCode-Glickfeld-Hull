clc
clear all
close all

areas = ["LM","AL","PM"];
nArea = length(areas);
cell_path = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\SizeTuning';
rc = behavConstsAV;
fout = fullfile(rc.ashleyAnalysis, 'SizeTuning');
%% Load files
fprintf('Loading cell and axon data')
load(fullfile(cell_path, 'cellBodyData.mat'))
load(fullfile(cell_path, 'cellBodyData_RF.mat'))
load(fullfile(fout, 'axonSizeSummary.mat'))

%% Compare V1 neurons vs axons

RFsize = [RFsize_V1; RFsize_axons_LM'; RFsize_axons_AL'; RFsize_axons_PM'];
RFsize_group = [ones(size(RFsize_V1)); 2.*ones(size(RFsize_axons_LM')); 3.*ones(size(RFsize_axons_AL')); 4.*ones(size(RFsize_axons_PM'))];
[p_RFsize_anova, t_RFsize, s_RFsize] = anova1(RFsize,RFsize_group,'off');
c_RFsize = multcompare(s_RFsize,'display','off');

prefSize = [prefSize_V1'; prefSize_axons_LM'; prefSize_axons_AL'; prefSize_axons_PM'];
prefSize_group = [ones(size(prefSize_V1')); 2.*ones(size(prefSize_axons_LM')); 3.*ones(size(prefSize_axons_AL')); 4.*ones(size(prefSize_axons_PM'))];
[p_prefSize_anova, t_prefSize, s_prefSize] = anova1(prefSize,prefSize_group,'off');
c_prefSize = multcompare(s_prefSize,'display','off');

suppInd = [suppInd_V1'; suppInd_axons_LM'; suppInd_axons_AL'; suppInd_axons_PM'];
suppInd_group = [ones(size(suppInd_V1')); 2.*ones(size(suppInd_axons_LM')); 3.*ones(size(suppInd_axons_AL')); 4.*ones(size(suppInd_axons_PM'))];
[p_suppInd_anova, t_suppInd, s_suppInd] = anova1(suppInd,suppInd_group,'off');
c_suppInd = multcompare(s_suppInd,'display','off');

figure;
subplot(3,1,1)
errorbar([1:4], [mean(RFsize_V1',2) mean(RFsize_axons_LM,2) mean(RFsize_axons_AL,2) mean(RFsize_axons_PM,2)],[std(RFsize_V1',[],2) std(RFsize_axons_LM,[],2) std(RFsize_axons_AL,[],2) std(RFsize_axons_PM,[],2)],'o')
ylabel('RF size (deg)')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([0 20])
title(['p = ' num2str(chop(c_RFsize(1:3,6)',2))])

subplot(3,1,2)
errorbar([1:4], [mean(prefSize_V1,2) mean(prefSize_axons_LM,2) mean(prefSize_axons_AL,2) mean(prefSize_axons_PM,2)],[std(prefSize_V1,[],2) std(prefSize_axons_LM,[],2) std(prefSize_axons_AL,[],2) std(prefSize_axons_PM,[],2)],'o')
ylabel('Pref size (deg)')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([0 40])
title(['p = ' num2str(chop(c_prefSize(1:3,6)',2))])

subplot(3,1,3)
errorbar([1:4], [mean(suppInd_V1,2) mean(suppInd_axons_LM,2) mean(suppInd_axons_AL,2) mean(suppInd_axons_PM,2)],[std(suppInd_V1,[],2) std(suppInd_axons_LM,[],2) std(suppInd_axons_AL,[],2) std(suppInd_axons_PM,[],2)],'o')
ylabel('Supp. Index')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([-.5 1.5])
title(['p = ' num2str(chop(c_suppInd(1:3,6)',2))])
suptitle('V1 Cells vs V1 Axons') 

print(fullfile(fout, 'Cell_Axon_V1comp.pdf'),'-dpdf', '-bestfit')

%% Compare HVA neurons vs axons
h_RFsize = zeros(1,nArea);
p_RFsize = zeros(1,nArea);
h_prefSize = zeros(1,nArea);
p_prefSize = zeros(1,nArea);
h_suppInd = zeros(1,nArea);
p_suppInd = zeros(1,nArea);
figure;
for i = 1:nArea
    [h_prefSize(i) p_prefsize(i)] = eval(strcat('ttest2(prefSize_axons_', areas(i), ',prefSize_' , areas(i), ')'));
    [h_suppInd(i) p_suppInd(i)] = eval(strcat('ttest2(suppInd_axons_', areas(i), ',suppInd_' , areas(i), ')'));
    [h_RFsize(i) p_RFsize(i)] = eval(strcat('ttest2(RFsize_axons_', areas(i), ',RFsize_' , areas(i), ')'));
    subplot(nArea,3,(i-1).*3 + 1)
    errorbar([1 2], [mean(eval(strcat('RFsize_', areas(i)))',2) mean(eval(strcat('RFsize_axons_', areas(i))),2)],[std(eval(strcat('RFsize_', areas(i)))',[],2) std(eval(strcat('RFsize_axons_', areas(i))),[],2)],'o');
    title([areas(i) 'p = ' chop(p_RFsize(i),2)])
    ylabel('RF size (deg)')
    xticks([1 2])
    xticklabels({'cells','axons'})
    xlim([0 3])
    ylim([0 30])
    subplot(nArea,3,(i-1).*3 + 2)
    errorbar([1 2], [mean(eval(strcat('prefSize_', areas(i))),2) mean(eval(strcat('prefSize_axons_', areas(i))),2)],[std(eval(strcat('prefSize_', areas(i))),[],2) std(eval(strcat('prefSize_axons_', areas(i))),[],2)],'o');
    title([areas(i) 'p = ' chop(p_prefSize(i),2)])
    ylabel('Pref size (deg)')
    xticks([1 2])
    xticklabels({'cells','axons'})
    xlim([0 3])
    ylim([0 80])
    subplot(nArea,3,(i-1).*3 + 3)
    errorbar([1 2], [mean(eval(strcat('suppInd_', areas(i))),2) mean(eval(strcat('suppInd_axons_', areas(i))),2)],[std(eval(strcat('suppInd_', areas(i))),[],2) std(eval(strcat('suppInd_axons_', areas(i))),[],2)],'o');
    title([areas(i) 'p = ' chop(p_suppInd(i),2)])
    ylabel('Supp. Index')
    xticks([1 2])
    xticklabels({'cells','axons'})
    xlim([0 3])
    ylim([-0.5 1.5])
end

print(fullfile(fout, 'Cell_Axon_HVAcomp.pdf'),'-dpdf', '-bestfit')
    

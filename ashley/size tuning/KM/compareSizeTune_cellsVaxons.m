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
load(fullfile(cell_path, 'allcells_goodfits.mat'))
load(fullfile(cell_path, 'cellBody_SI_prefSize.mat'))
load(fullfile(fout, 'axonSizeSummary.mat'))

%% Sort cell groups
RFsize_V1 = RFsize_all(intersect(goodfit_ind_all, find(areaInd==1)),:).*2;
RFsize_LM = RFsize_all(intersect(goodfit_ind_all, find(areaInd==2)),:).*2;
RFsize_AL = RFsize_all(intersect(goodfit_ind_all, find(areaInd==3)),:).*2;
RFsize_PM = RFsize_all(intersect(goodfit_ind_all, find(areaInd==4)),:).*2;

% V1_ind= find(areaInd(goodfit_ind_all)==1);
% LM_ind= find(areaInd(goodfit_ind_all)==2);
% AL_ind= find(areaInd(goodfit_ind_all)==3);
% PM_ind= find(areaInd(goodfit_ind_all)==4);
% 
% prefSize = zeros(1,length(goodfit_ind_all));
% for i = 1: length(goodfit_ind_all)
%     prefSize(1,i) = sizeFits_all(i,4).prefSize;
% end

prefSize_V1 = yPS(find(x==1))';
prefSize_LM = yPS(find(x==2))';
prefSize_AL = yPS(find(x==3))';
prefSize_PM = yPS(find(x==4))';

% suppInd = zeros(1,length(goodfit_ind_all));
% for i = 1: length(goodfit_ind_all)
%     suppInd(1,i) = sizeFits_all(i,4).suppInd;
% end

suppInd_V1 = ySI(find(x==1))';
suppInd_LM = ySI(find(x==2))';
suppInd_AL = ySI(find(x==3))';
suppInd_PM = ySI(find(x==4))';


%% Compare V1 neurons vs axons

RFsize = [RFsize_V1; RFsize_axons_LM'; RFsize_axons_AL'; RFsize_axons_PM'];
RFsize_group = [ones(size(RFsize_V1)); 2.*ones(size(RFsize_axons_LM')); 3.*ones(size(RFsize_axons_AL')); 4.*ones(size(RFsize_axons_PM'))];
[p_RFsize_anova, t_RFsize, s_RFsize] = kruskalwallis(RFsize,RFsize_group,'off');
c_RFsize = multcompare(s_RFsize,'display','off');

prefSize = [prefSize_V1'; prefSize_axons_LM'; prefSize_axons_AL'; prefSize_axons_PM'];
prefSize_group = [ones(size(prefSize_V1')); 2.*ones(size(prefSize_axons_LM')); 3.*ones(size(prefSize_axons_AL')); 4.*ones(size(prefSize_axons_PM'))];
[p_prefSize_anova, t_prefSize, s_prefSize] = kruskalwallis(prefSize,prefSize_group,'off');
c_prefSize = multcompare(s_prefSize,'display','off');

suppInd = [suppInd_V1'; suppInd_axons_LM'; suppInd_axons_AL'; suppInd_axons_PM'];
suppInd_group = [ones(size(suppInd_V1')); 2.*ones(size(suppInd_axons_LM')); 3.*ones(size(suppInd_axons_AL')); 4.*ones(size(suppInd_axons_PM'))];
[p_suppInd_anova, t_suppInd, s_suppInd] = kruskalwallis(suppInd,suppInd_group,'off');
c_suppInd = multcompare(s_suppInd,'display','off');

n1 = length(find(suppInd_V1==0)); N1 = length(suppInd_V1);
n2 = length(find(suppInd_axons_LM==0)); N2 =length(suppInd_axons_LM);
n3 = length(find(suppInd_axons_AL==0)); N3 =length(suppInd_axons_AL);
n4 = length(find(suppInd_axons_PM==0)); N4 =length(suppInd_axons_PM);
x1 = [repmat('a',N1,1); repmat('b',N2,1); repmat('c',N3,1); repmat('d',N4,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1); repmat(1,n3,1); repmat(2,N3-n3,1); repmat(1,n4,1); repmat(2,N4-n4,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2);

n1 = length(find(suppInd_V1==0)); N1 = length(suppInd_V1);
n2 = length(find(suppInd_axons_LM==0)); N2 =length(suppInd_axons_LM);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_V1vsLM] = crosstab(x1,x2);

n1 = length(find(suppInd_V1==0)); N1 = length(suppInd_V1);
n2 = length(find(suppInd_axons_AL==0)); N2 =length(suppInd_axons_AL);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_V1vsAL] = crosstab(x1,x2);

n1 = length(find(suppInd_V1==0)); N1 = length(suppInd_V1);
n2 = length(find(suppInd_axons_PM==0)); N2 =length(suppInd_axons_PM);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_V1vsPM] = crosstab(x1,x2);

n1 = length(find(suppInd_V1==0)); N1 = length(suppInd_V1);
n2 = length(find(suppInd_axons_PM==0))+length(find(suppInd_axons_AL==0))+length(find(suppInd_axons_LM==0)); N2 =length(suppInd_axons_PM)+length(suppInd_axons_AL)+length(suppInd_axons_LM);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_all] = crosstab(x1,x2);


figure;
subplot(3,1,1)
errorbar([1:4], [mean(RFsize_V1',2) mean(RFsize_axons_LM,2) mean(RFsize_axons_AL,2) mean(RFsize_axons_PM,2)],[std(RFsize_V1',[],2) std(RFsize_axons_LM,[],2) std(RFsize_axons_AL,[],2) std(RFsize_axons_PM,[],2)],'o')
ylabel('RF size (deg)')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([0 30])
title(['p = ' num2str(chop([p_RFsize_anova c_RFsize(1:3,6)'],2))])

subplot(3,1,2)
errorbar([1:4], [mean(prefSize_V1,2) mean(prefSize_axons_LM,2) mean(prefSize_axons_AL,2) mean(prefSize_axons_PM,2)],[std(prefSize_V1,[],2) std(prefSize_axons_LM,[],2) std(prefSize_axons_AL,[],2) std(prefSize_axons_PM,[],2)],'o')
ylabel('Pref size (deg)')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([0 40])
title(['p = ' num2str(chop([p_prefSize_anova c_prefSize(1:3,6)'],2))])

subplot(3,1,3)
errorbar([1:4], [mean(suppInd_V1,2) mean(suppInd_axons_LM,2) mean(suppInd_axons_AL,2) mean(suppInd_axons_PM,2)],[std(suppInd_V1,[],2) std(suppInd_axons_LM,[],2) std(suppInd_axons_AL,[],2) std(suppInd_axons_PM,[],2)],'o')
ylabel('Supp. Index')
xticks([1:4])
xticklabels({'V1 cells','V1->LM','V1->AL','V1->PM'})
xlim([0 5])
ylim([-.5 1.5])
title(['p = ' num2str(chop([p_suppInd_anova c_suppInd(1:3,6)'],2))])
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
    [p_prefSize(i), h_prefsize(i)] = eval(strcat('ranksum(prefSize_axons_', areas(i), ',prefSize_' , areas(i), ')'));
    [p_suppInd(i) h_suppInd(i)] = eval(strcat('ranksum(suppInd_axons_', areas(i), ',suppInd_' , areas(i), ')'));
    [p_RFsize(i) h_RFsize(i)] = eval(strcat('ranksum(RFsize_axons_', areas(i), ',RFsize_' , areas(i), ')'));
    subplot(nArea,3,(i-1).*3 + 1)
    errorbar([1 2], [mean(eval(strcat('RFsize_', areas(i)))',2) mean(eval(strcat('RFsize_axons_', areas(i))),2)],[std(eval(strcat('RFsize_', areas(i)))',[],2) std(eval(strcat('RFsize_axons_', areas(i))),[],2)],'o');
    hold on
    text(1.5, 28, num2str(chop(mean(eval(strcat('RFsize_', areas(i)))',2),3)));
    text(1.5, 18, num2str(chop(mean(eval(strcat('RFsize_axons_', areas(i))),2),3)));
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

[p_RFsize_LMleft, h_RFsize_LMleft] = ranksum(RFsize_axons_LM, RFsize_LM, 'tail', 'left');
[p_RFsize_ALleft, h_RFsize_ALleft] = ranksum(RFsize_axons_AL, RFsize_AL, 'tail', 'left');
[p_RFsize_PMleft, h_RFsize_PMleft] = ranksum(RFsize_axons_PM, RFsize_PM, 'tail', 'left');

[p_prefSize_LMleft, h_prefsize_LMleft] = ranksum(prefSize_axons_LM, prefSize_LM, 'tail', 'left');
[p_prefSize_ALleft, h_prefsize_ALleft] = ranksum(prefSize_axons_AL, prefSize_AL, 'tail', 'left');
[p_prefSize_PMleft, h_prefsize_PMleft] = ranksum(prefSize_axons_PM, prefSize_PM, 'tail', 'left');


n1 = length(find(suppInd_LM==0)); N1 = length(suppInd_LM);
n2 = length(find(suppInd_axons_LM==0)); N2 =length(suppInd_axons_LM);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_LM] = crosstab(x1,x2);

n1 = length(find(suppInd_AL==0)); N1 = length(suppInd_AL);
n2 = length(find(suppInd_axons_AL==0)); N2 =length(suppInd_axons_AL);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_AL] = crosstab(x1,x2);

n1 = length(find(suppInd_PM==0)); N1 = length(suppInd_PM);
n2 = length(find(suppInd_axons_PM==0)); N2 =length(suppInd_axons_PM);
x1 = [repmat('a',N1,1); repmat('b',N2,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
[tbl,chi2stat,pval_PM] = crosstab(x1,x2);

% 
% n3 = length(find(suppInd_AL==0)); N3 =length(suppInd_AL);
% n4 = length(find(suppInd_axons_AL==0)); N4 =length(suppInd_axons_AL);
% n5 = length(find(suppInd_PM==0)); N3 =length(suppInd_PM);
% n6 = length(find(suppInd_axons_PM==0)); N4 =length(suppInd_axons_PM);
% x1 = [repmat('a',N1,1); repmat('b',N4,1); repmat('c',N3,1); repmat('d',N4,1)];
% [tbl,chi2stat,pval] = crosstab(x1,x2);


print(fullfile(fout, 'Cell_Axon_HVAcomp.pdf'),'-dpdf', '-bestfit')
    

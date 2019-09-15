ds = 'FS_HVA';
rc = behavConstsAV;
eval(ds)

mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

ori_fn = fullfile(rc.caOutputDir,ds,[mouse_str '_tuningData_cells_' ds '.mat']); 
load(ori_fn)
FS_fn = fullfile(rc.caOutputDir,ds,[mouse_str '_FSData_cells_' ds '.mat']);
load(FS_fn)
decode_fn = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeStruct_cells_' ds '.mat']);
load(decode_fn)
fnout = fullfile(rc.caOutputDir,ds);

area_list = {'LM','AL','PM'};
narea = length(area_list);
nexp = size(decodeDataExpt,2);

OSI_all = cell(1,narea);
nCells_all = zeros(1,narea);
nWellfitCells_all = zeros(1,narea);
wellfitCells_all = cell(1,narea);
auROC_all = cell(1,narea);
for iexp = 1:nexp
    iarea = find(strcmp(area_list,decodeDataExpt(iexp).exptArea));
    nCells_all(:,iarea) = nCells_all(:,iarea)+oriData(iexp).nCells;
    nWellfitCells_all(:,iarea) = nWellfitCells_all(:,iarea)+oriData(iexp).nWellfitCells;
    OSI_all{iarea} = [OSI_all{iarea} oriData(iexp).OSI];
    wellfitCells_all{iarea} = [wellfitCells_all{iarea} oriData(iexp).wellfitCells];
    auROC_all{iarea} = [auROC_all{iarea}; FSData(iexp).auROC];
end

figure;
for iarea = 1:narea
    subplot(2,2,1)
    bar(iarea, nCells_all(iarea))
    hold on
    subplot(2,2,2)
    bar(iarea, nWellfitCells_all(iarea)./nCells_all(iarea))
    hold on
    [pct ci] = binofit(nWellfitCells_all(iarea), nCells_all(iarea));
    errorbar(iarea, pct, pct-ci(1), ci(2)-pct, 'k')
    subplot(2,2,3)
    cdfplot(OSI_all{iarea})
    hold on
    subplot(2,2,4)
    cdfplot(auROC_all{iarea}(:,end))
    hold on
end
subplot(2,2,1)
set(gca, 'XTick', 1:3, 'XTickLabel',area_list)
ylabel('Number of cells')
title('Total number of segmented cells')
subplot(2,2,2)
set(gca,'XTick', 1:3, 'XTickLabel',area_list)
ylim([0 1])
ylabel('Fraction of cells')
title('Fraction of well-fit cells')
subplot(2,2,3)
xlim([0 1])
xlabel('OSI')
ylabel('Cumulative fraction')
[LMvAL_OSI_H LMvAL_OSI_P] = kstest2(OSI_all{1},OSI_all{2});
[LMvPM_OSI_H LMvPM_OSI_P] = kstest2(OSI_all{1},OSI_all{3});
[PMvAL_OSI_H PMvAL_OSI_P] = kstest2(OSI_all{3},OSI_all{2});
title(num2str([LMvAL_OSI_P  LMvPM_OSI_P PMvAL_OSI_P]))
subplot(2,2,4)
vline(0.5)
xlim([0 1])
xlabel('auROC: 0 vs 90')
ylabel('Cumulative fraction')
[LMvAL_auROC_H LMvAL_auROC_P] = kstest2(auROC_all{1}(:,end),auROC_all{2}(:,end));
[LMvPM_auROC_H LMvPM_auROC_P] = kstest2(auROC_all{1}(:,end),auROC_all{3}(:,end));
[PMvAL_auROC_H PMvAL_auROC_P] = kstest2(auROC_all{3}(:,end),auROC_all{2}(:,end));
title(num2str([LMvAL_auROC_P  LMvPM_auROC_P PMvAL_auROC_P]))
suptitle(mouse_str)
print(fullfile(fnout,[mouse_str '_tuningSummary_' ds '.pdf']))

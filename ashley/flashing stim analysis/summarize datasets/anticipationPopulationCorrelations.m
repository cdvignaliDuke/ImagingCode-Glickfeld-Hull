clear all
close all
dataGroup = 'awFSAVdatasets_PM';
eval(dataGroup);
cellsetname = 'nCells';
av = behavParamsAV;
% CYC = 6;
trialLengthMs_gt = 2500;


for iexp = 1:length(expt)
    eval(dataGroup);
    imouse = find(strcmp(cellfun(@num2str,{av.mouse},'UniformOutput',0), expt(iexp).SubNum));
    col = av(imouse).col_str;
    colS = av(imouse).sec_col_str;
%%
ialign = 1;
run('divideupdatabyalignment.m')
%%
DirFolder = expt(iexp).dirtuning;
run('cellSetsLG.m')
if length(eval(cellsetname)) == 1
    cells = 1:eval(cellsetname);
else
    cells = eval(cellsetname);
end
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\summary figs' ];
exptTag = [SubNum '_' dateofexp '_'];
try 
    cd(fnout)
catch
    mkdir(['Z:\Analysis\' mouse '\two-photon imaging\'],'summary figs')
end

%%
frameratems = expt(iexp).frame_rate/1000;
cyctimems = cycTime/frameratems;
stimon_frames = input.nFramesOn;
stimoff_frames = input.nFramesOff;
stimontime = stimon_frames/frameratems;
stimofftime = stimoff_frames/frameratems;
trialresult = unique(trialOutcome);
responsedelay_frames = 3;

CYC = ceil(trialLengthMs_gt/(stimontime+stimofftime));
trialLengthFrames = cycTime*CYC;
xAxisMs = (-(prepress_frames-1):trialLengthFrames)/frameratems;

data = cycDataDFoverF_cmlvNoTarget{CYC};

V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
    
    
%% calc integral

dataInt = squeeze(trapz(data(prepress_frames:end,:,:)))';
dataInt_first = squeeze(trapz(data(prepress_frames:floor(trialLengthFrames/2)+prepress_frames,:,:)))';
dataInt_last = squeeze(trapz(data(floor(trialLengthFrames/2)+1:end,:,:)))';

dataInt_V = squeeze(mean(dataInt(V_cycInd,:),1));
dataInt_AV = squeeze(mean(dataInt(AV_cycInd,:),1));

% sort cells by magnitude of response to each trial type
[meanIntV meanIntSortInd_V] = sort(dataInt_V);
[meanIntAV meanIntSortInd_AV] = sort(dataInt_AV);

%% sort cells by orientation preference
[oriPrefSort oriPrefCellsSortInd] = sort(oriPref_ind);
oriGroupBorders = cumsum(cell2mat(cellfun(@length,cellsPref,'Unif',false)));
oriGroupLabels = cellfun(@num2str,num2cell(Oris),'Unif',false);

%% noise correlation of integral of response of all neurons during anticipation period
cellSortOrder = eval('oriPrefCellsSortInd');
cellGroupTicks = eval('oriGroupBorders');
cellGroupLabels = eval('oriGroupLabels');

% whole trial
nsCorrIntHeatMaps_whole = figure;
colormap(brewermap([],'*RdBu'))
suptitle([mouse '-' date '; noise correlation, whole trial'])

nsCorrIntMat_V = tril(corrcoef(dataInt(V_cycInd,cellSortOrder)),-1);
nsCorrIntMean_V{iexp} = mean(nsCorrIntMat_V(nsCorrIntMat_V ~= 0));

figure(nsCorrIntHeatMaps_whole);
subplot(2,2,1)
imagesc(nsCorrIntMat_V)
title(['vis tr - ' num2str(nsCorrIntMean_V{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)


nsCorrIntMat_AV = tril(corrcoef(dataInt(AV_cycInd,cellSortOrder)),-1);
nsCorrIntMean_AV{iexp} = mean(nsCorrIntMat_AV(nsCorrIntMat_V ~= 0));

figure(nsCorrIntHeatMaps_whole);
subplot(2,2,2)
imagesc(nsCorrIntMat_AV);
title(['aud tr - ' num2str(nsCorrIntMean_AV{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)

% first half
nsCorrIntHeatMaps_1half = figure;
colormap(brewermap([],'*RdBu'))
suptitle([mouse '-' date '; noise correlation, 1st half'])

nsCorrIntMat1_V = tril(corrcoef(dataInt_first(V_cycInd,cellSortOrder)),-1);
nsCorrIntMean1_V{iexp} = mean(nsCorrIntMat1_V(nsCorrIntMat1_V ~= 0));

figure(nsCorrIntHeatMaps_1half);
subplot(2,2,1)
imagesc(nsCorrIntMat1_V)
title(['vis tr - ' num2str(nsCorrIntMean1_V{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)


nsCorrIntMat1_AV = tril(corrcoef(dataInt_first(AV_cycInd,cellSortOrder)),-1);
nsCorrIntMean1_AV{iexp} = mean(nsCorrIntMat1_AV(nsCorrIntMat1_V ~= 0));

figure(nsCorrIntHeatMaps_1half);
subplot(2,2,2)
imagesc(nsCorrIntMat1_AV);
title(['aud tr - ' num2str(nsCorrIntMean1_AV{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)

% second half
nsCorrIntHeatMaps_2half = figure;
colormap(brewermap([],'*RdBu'))
suptitle([mouse '-' date '; noise correlation, 2nd half'])

nsCorrIntMat2_V = tril(corrcoef(dataInt_last(V_cycInd,cellSortOrder)),-1);
nsCorrIntMean2_V{iexp} = mean(nsCorrIntMat2_V(nsCorrIntMat2_V ~= 0));

figure(nsCorrIntHeatMaps_2half);
subplot(2,2,1)
imagesc(nsCorrIntMat2_V)
title(['vis tr - ' num2str(nsCorrIntMean2_V{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)


nsCorrIntMat2_AV = tril(corrcoef(dataInt_last(AV_cycInd,cellSortOrder)),-1);
nsCorrIntMean2_AV{iexp} = mean(nsCorrIntMat2_AV(nsCorrIntMat2_V ~= 0));

figure(nsCorrIntHeatMaps_2half);
subplot(2,2,2)
imagesc(nsCorrIntMat2_AV);
title(['aud tr - ' num2str(nsCorrIntMean2_AV{iexp})])
colorbar
caxis([-1 1])
axis square
for i = 1:length(cellGroupTicks)
    vline(cellGroupTicks(i),'k')
end
for i = 1:length(cellGroupTicks)
    hline(cellGroupTicks(i),'k')
end
set(gca,'XTick',cellGroupTicks)
set(gca,'XTickLabel',cellGroupLabels)

%% plot noise correlation as a function of mean stimulus modulation across a pair of neurons
dataIntPairs_V = tril(bsxfun(@plus, repmat(dataInt_V,[length(dataInt_V), 1]),dataInt_V')/2,-1);
dataIntPairs_AV = tril(bsxfun(@plus, repmat(dataInt_AV,[length(dataInt_AV), 1]),dataInt_AV')/2,-1);

%************* whole trial
figure(nsCorrIntHeatMaps_whole);
subplot(2,2,3)
scatter(dataIntPairs_V(dataIntPairs_V ~= 0),nsCorrIntMat_V(nsCorrIntMat_V ~= 0),'g.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('vis tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_V(dataIntPairs_V ~= 0), nsCorrIntMat_V(nsCorrIntMat_V ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_V(dataIntPairs_V ~= 0)), max(dataIntPairs_V(dataIntPairs_V ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

figure(nsCorrIntHeatMaps_whole);
subplot(2,2,4)
scatter(dataIntPairs_AV(dataIntPairs_AV ~= 0),nsCorrIntMat_AV(nsCorrIntMat_AV ~= 0),'k.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('aud tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_AV(dataIntPairs_AV ~= 0), nsCorrIntMat_AV(nsCorrIntMat_AV ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_AV(dataIntPairs_AV ~= 0)), max(dataIntPairs_AV(dataIntPairs_AV ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%********** first half
figure(nsCorrIntHeatMaps_1half);
subplot(2,2,3)
scatter(dataIntPairs_V(dataIntPairs_V ~= 0),nsCorrIntMat1_V(nsCorrIntMat1_V ~= 0),'g.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('vis tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_V(dataIntPairs_V ~= 0), nsCorrIntMat1_V(nsCorrIntMat1_V ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_V(dataIntPairs_V ~= 0)), max(dataIntPairs_V(dataIntPairs_V ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

figure(nsCorrIntHeatMaps_1half);
subplot(2,2,4)
scatter(dataIntPairs_AV(dataIntPairs_AV ~= 0),nsCorrIntMat1_AV(nsCorrIntMat1_AV ~= 0),'k.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('aud tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_AV(dataIntPairs_AV ~= 0), nsCorrIntMat1_AV(nsCorrIntMat1_AV ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_AV(dataIntPairs_AV ~= 0)), max(dataIntPairs_AV(dataIntPairs_AV ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%*************8 second half
figure(nsCorrIntHeatMaps_2half);
subplot(2,2,3)
scatter(dataIntPairs_V(dataIntPairs_V ~= 0),nsCorrIntMat2_V(nsCorrIntMat2_V ~= 0),'g.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('vis tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_V(dataIntPairs_V ~= 0), nsCorrIntMat2_V(nsCorrIntMat2_V ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_V(dataIntPairs_V ~= 0)), max(dataIntPairs_V(dataIntPairs_V ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

figure(nsCorrIntHeatMaps_2half);
subplot(2,2,4)
scatter(dataIntPairs_AV(dataIntPairs_AV ~= 0),nsCorrIntMat2_AV(nsCorrIntMat2_AV ~= 0),'k.')
hold on
vline(0,'k:')
hold on
hline(0,'k:')
ylim([-1 1])
title('aud tr - across pairs of neurons')
ylabel('ns corr')
xlabel('mean resp')

coeffs = polyfit(dataIntPairs_AV(dataIntPairs_AV ~= 0), nsCorrIntMat2_AV(nsCorrIntMat2_AV ~= 0), 1);
% Get fitted values
fittedX = linspace(min(dataIntPairs_AV(dataIntPairs_AV ~= 0)), max(dataIntPairs_AV(dataIntPairs_AV ~= 0)), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);

%% plot avg ns corr by ori preference


%% save figs

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

figure(nsCorrIntHeatMaps_whole);
print(fullfile(fnout,[exptTag 'nsCorrHeatMapsFRCorr_' cellsetname]), '-dpdf');
figure(nsCorrIntHeatMaps_1half);
print(fullfile(fnout,[exptTag 'nsCorrHeatMapsFRCorr_1half_' cellsetname]), '-dpdf');
figure(nsCorrIntHeatMaps_2half);
print(fullfile(fnout,[exptTag 'nsCorrHeatMapsFRCorr_2half_' cellsetname]), '-dpdf');

end

%% plot average ns corr for each experiment
for iexp = 1:size(expt,2)
   imouse = find(strcmp(cellfun(@num2str,{av.mouse},'UniformOutput',0), expt(iexp).SubNum));
   cStr{iexp} = av(imouse).col_str;
   mouseMat(iexp) = av(imouse).mouse;
end
nsCorrMeanSummary = figure;
suptitle(dataGroup)
mice = unique(mouseMat);

for iMs = 1:length(mice)
    imouse = find(cell2mat({av.mouse}) == mice(iMs));
    msInd = find(mouseMat == mice(iMs));
    msMean_V = mean(cell2mat(nsCorrIntMean2_V(msInd)));
    msMean_AV = mean(cell2mat(nsCorrIntMean2_AV(msInd)));
    msSte_V = std(cell2mat(nsCorrIntMean2_V(msInd)))/(length(msInd));
    msSte_AV = std(cell2mat(nsCorrIntMean2_AV(msInd)))/(length(msInd));
    
    colorCell = {[av(imouse).col_str 'o'],av(imouse).col_str,av(imouse).col_str};
%     errorbarxy(cell2mat(nsCorrIntMean_V(msInd)),cell2mat(nsCorrIntMean_AV(msInd)),[],[],colorCell)
    plot(cell2mat(nsCorrIntMean2_V(msInd)),cell2mat(nsCorrIntMean2_AV(msInd)),[av(imouse).col_str 'o'])
    hold on
    errorbarxy(msMean_V,msMean_AV,msSte_V,msSte_AV,colorCell)
    hold on
    msColMean(iMs) = scatter(msMean_V,msMean_AV,[av(imouse).col_str 'o'],'filled');
    summaryLegend{iMs} = num2str(mice(iMs));
    
end
plot([-1:0.1:1],[-1:0.1:1],'k--')
hold on
xlim([-1 1]);
ylim([-1 1]);
axis square
xlabel('Vis Tr Corr')
ylabel('Aud Tr Corr')
legend(msColMean,summaryLegend,'Location','southeast')
title('mean ns corr (2nd half)')

%% save summary fig
figure(nsCorrMeanSummary)
print(fullfile('Z:\Analysis\FSAV Summaries',dataGroup,'nsCorrMeanSummary'),'-dpdf')

ialign = 1;
awFSAVdatasets;
for iexp = 1:size(expt,2)
% iexp = 5;
divideupdatabyalignment
%%

fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetNoiseCorr'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetNoiseCorr')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetNoiseCorr')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% find sets of cells
DirFolder = expt(iexp).dirtuning;
run('cellSets.m')

%%

%%
cells = 1:nCells;
cellgroupname = 'all cells';
figBaseName = 'press2pre-target_integral_noiseCorr_allcells_success';
%% set data
CYC = 6;

cycInd = find(tCyclesOn == cycles(CYC));
successIx = find(ismember(cycInd,find(strcmp(trialOutcome,'success'))));
cycBlock2 = block2(cycInd);
trsV = intersect(find(cycBlock2 == 0),successIx);
trsAV = intersect(find(cycBlock2 == 1),successIx);

tempdfoverf = cycDataDFoverF{CYC};
tempDataIntegral = squeeze(trapz(tempdfoverf(30:end,:,:)));
% tempDataIntegral(tempDataIntegral < 0) = 0;


%sort data by ori pref
[sortedOriPref, sortedOriPref_ind] = sort(oriPref_ind(cells,:));
dataSelectCells = tempDataIntegral(cells,:);
dataSortedByOriPref = dataSelectCells(sortedOriPref_ind,:);

dataSortedMean_V = mean(dataSortedByOriPref(:,trsV),2);
dataSortedMean_AV = mean(dataSortedByOriPref(:,trsAV),2);

%% find mean FR between pairs of neurons
allPairsMeanFR_V = bsxfun(@plus,(ones(length(dataSortedMean_V),1)*dataSortedMean_V'),dataSortedMean_V)/2;
allPairsMeanFR_AV = bsxfun(@plus,(ones(length(dataSortedMean_AV),1)*dataSortedMean_AV'),dataSortedMean_AV)/2;;

allPairsMeanFR_V = tril(allPairsMeanFR_V,-1);
allPairsMeanFR_AV = tril(allPairsMeanFR_AV,-1);

allPairsMeanFR_Vv = allPairsMeanFR_V(allPairsMeanFR_V  ~= 0);
allPairsMeanFR_AVv = allPairsMeanFR_AV(allPairsMeanFR_AV  ~= 0);
%% noise correlations last 3cyles of long trials

%find correlation coefficient (r) for cells
rSC_V= corrcoef(dataSortedByOriPref(:,trsV)');
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(dataSortedByOriPref(:,trsAV)');
rSC_sub = rSC_AV-rSC_V;

rSC_V = tril(rSC_V,-1);
rSC_AV = tril(rSC_AV,-1);
rSC_sub = tril(rSC_sub,-1);

avgrSC_V = mean(mean(rSC_V(rSC_V  ~= 0),2),1);
avgrSC_AV = mean(mean(rSC_AV(rSC_AV  ~= 0),2),1);
avgrSC_sub = mean(mean(rSC_sub(rSC_sub  ~= 0),2),1);


%% bin rSC by firing rate
nbins = 5;
rSC_Vv = rSC_V(rSC_V  ~= 0);
rSC_AVv = rSC_AV(rSC_AV  ~= 0);

FR_linbins = linspace(min(cat(1,allPairsMeanFR_Vv,allPairsMeanFR_AVv)),max(cat(1,allPairsMeanFR_Vv,allPairsMeanFR_AVv)),nbins);
FR_logbins = logspace(min(cat(1,allPairsMeanFR_Vv,allPairsMeanFR_AVv)),max(cat(1,allPairsMeanFR_Vv,allPairsMeanFR_AVv)),nbins);

[nVperbin Vbin_ind] = histc(allPairsMeanFR_Vv,FR_linbins);
[nAVperbin AVbin_ind] = histc(allPairsMeanFR_Vv,FR_linbins);

rSCbinMean_V = zeros(1,nbins);
rSCbinMean_AV = zeros(1,nbins);
rSCbinErr_V = zeros(1,nbins);
rSCbinErr_AV = zeros(1,nbins);
for ibin = 1:nbins
    ind = find(Vbin_ind == ibin);
    rSCbinMean_V(1,ibin) = mean(rSC_Vv(ind));
    rSCbinMean_AV(1,ibin) = mean(rSC_AVv(ind));
    rSCbinErr_V(1,ibin) = std(rSC_Vv(ind))/sqrt(length(ind));
    rSCbinErr_AV(1,ibin) = std(rSC_AVv(ind))/sqrt(length(ind));
    clear ind
end

yminmax = [(min(cat(2,rSCbinMean_V,rSCbinMean_AV))-.2) (max(cat(2,rSCbinMean_V,rSCbinMean_AV))+.2)]; 

figure;
errorbar(FR_linbins,rSCbinMean_V,rSCbinErr_V,'g')
hold on
errorbar(FR_linbins,rSCbinMean_AV,rSCbinErr_AV,'k')
hold on
xlabel('FR bin edge')
ylabel('nCorr')
ylim(yminmax)
title({[mouse '; ' date '; ' cellgroupname];['ns corr as fxn of FR (dF/F integral)'];[num2str(length(trsV)) ' vis, ' num2str(length(trsAV)) 'aud; ' num2str(length(cells)) ' cells']})
set(gca,'XTick',FR_linbins)
set(gca,'XTickLabel',FR_linbins)
hold on
strV = strcat(repmat({'\leftarrow n='},nbins,1),arrayfun(@num2str,nVperbin,'UniformOutput',false));
text(FR_linbins,rSCbinMean_V,strV)
hold on
strAV = strcat(repmat({'n='},nbins,1),arrayfun(@num2str,nAVperbin,'UniformOutput',false));
text(FR_linbins,rSCbinMean_AV,strAV)


%% plot noise correlations first 3cycles vs last 3cyles of long trials
oriGroupBorders = zeros(size(unique(sortedOriPref)))';
runsum = 0;
for i = 1:length(oriGroupBorders)
    x = sum(sortedOriPref == i);
    oriGroupBorders(1,i) = x+runsum;
    runsum = runsum+x;
end

%plot noise correlation
figure;
suptitle({[mouse '; ' date '; ' cellgroupname 'noise corr of press to pre-target integral'];[num2str(length(trsV)) ' vis, ' num2str(length(trsAV)) 'aud; ' num2str(length(cells)) ' cells']})
colormap(brewermap([],'*RdBu'))
subplot(1,3,1);
imagesc(rSC_V)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V)])
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)
subplot(1,3,2);
imagesc(rSC_AV)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV)])
hold on
colorbar
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)
subplot(1,3,3);
imagesc(rSC_sub)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub)])
hold on
colorbar
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)

% print([fnout ['\' figBaseName '_heatmap' '.pdf']], '-dpdf')
end

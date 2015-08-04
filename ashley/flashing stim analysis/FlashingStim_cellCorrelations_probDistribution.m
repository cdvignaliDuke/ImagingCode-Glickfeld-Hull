edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

dirFolder = '006';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')

%%
figure;
for i = 1:length(cycles(4:end))
    icyc = cycles(i+3);
dataTrialStart = cycDataDFoverF_cmlvNoTarget{icyc};
V_cycInd = cycV_ind{icyc};
% A_cycInd = cycA_ind{icyc};
AV_cycInd = cycAV_ind{icyc};


V_cellAvgResp = squeeze(mean(dataTrialStart(:,baselineStimRespIndex_V,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataTrialStart(:,cellsPrefZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataTrialStart(:,baselineStimRespIndex_V,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);
rSC_sub = rSC_AV-rSC_V;

avgrSC_V = mean(mean(rSC_V,2),1);
avgrSC_AV = mean(mean(rSC_AV,2),1);
avgrSC_sub = mean(mean(rSC_sub,2),1);

%vector for correlations witout 1s (e.g. cell#1xcell#1) or duplicates
rSC_Vsingle = tril(rSC_V,-1);
rSC_AVsingle = tril(rSC_AV,-1);
rSC_subSingle = tril(rSC_sub,-1);

rSC_Vvector = rSC_Vsingle(:);
rSC_Vvector = rSC_Vvector(rSC_Vvector ~= 0);
rSC_AVvector = rSC_AVsingle(:);
rSC_AVvector = rSC_AVvector(rSC_AVvector ~= 0);
rSC_subVector = rSC_subSingle(:);
rSC_subVector = rSC_subVector(rSC_subVector ~= 0);

avgrSC_V = mean(rSC_Vvector);
avgrSC_AV = mean(rSC_AVvector);
avgrSC_sub = mean(rSC_subVector);

%histogram for each correlation matrix
subplot(2,4,i)
xRange = -1:0.1:1;
rSC_Vbins = hist(rSC_Vvector,xRange);
rSC_AVbins = hist(rSC_AVvector,xRange);
plot(xRange,rSC_Vbins./numel(rSC_Vbins),'g');
hold on
plot(xRange,rSC_AVbins./numel(rSC_AVbins),'k');
hold on
vline(avgrSC_V,'g');
hold on
vline(avgrSC_AV,'k');
hold on
xlabel('corrcoef')
ylabel('prob')

% %cumulative distribution plots
% subplot(2,4,i)
% p1 = cdfplot(rSC_Vvector);
% set(p1,'color','g')
% hold on
% hline(avgrSC_V,'g');
% hold on
% p2 = cdfplot(rSC_AVvector);
% set(p2,'color','k')
% hold on
% hline(avgrSC_AV,'k');

% hist(rSC_Vvector,100);
% alpha(0.25)
% hold on
% hist(rSC_AVvector,100);
% h = findobj(gca,'Type','patch');
% set(h(2),'FaceColor','g','EdgeColor','none')
% set(h(1),'FaceColor','k','EdgeColor','none')
% alpha(0.25)
% xlim([-1 1])
% vline(0,'k:')
end


%%

figure;
colormap(brewermap([],'*RdBu'))
subplot(1,3,1);
imagesc(rSC_Vsingle)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V)])
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
% subplot(1,3,2);
% imagesc(rSC_A)
% title('rSC_A')
% hold on
% colorbar
% caxis([-1 1]);
% axis('square')
subplot(1,3,2);
imagesc(rSC_AVsingle)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV)])
hold on
colorbar
caxis([-1 1]);
axis('square')
subplot(1,3,3);
imagesc(rSC_subSingle)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub)])
hold on
colorbar
caxis([-1 1]);
axis('square')
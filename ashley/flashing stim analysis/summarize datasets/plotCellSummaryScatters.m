close all
clear all
% iexp = 1;
% fnouttemp = ['Z:\Analysis\temp figs\150924retreatposter']; %
awFSAVdatasets
av = behavParamsAV;
analysisName = {'PreTarget_respIntegral_Cells'; 'PreTarget_respAmp2basestim1_Cells'; 'PreTarget_respAmp2basestimlast_Cells';'PreTarget_std_Cells'};
figName = {'integralScatterFig';'respfirststimScatterFig';'resplaststimScatterFig';'stdScatterFig'};

CYC = 6;

integralScatterFig = figure;
respfirststimScatterFig = figure;
resplaststimScatterFig = figure;
stdScatterFig = figure;

for iexp = 1:size(expt,2)
    imouse = find(strcmp(cellfun(@num2str,{av.mouse},'UniformOutput',0), expt(iexp).SubNum));
    col = av(imouse).col_str;
    colS = av(imouse).sec_col_str;
%%
ialign = 1;
run('divideupdatabyalignment.m')
%%
DirFolder = expt(iexp).dirtuning;
run('cellSetsLG.m')
cells = 1:nCells;
%%

fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' dateofexp '\' runstr];

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

%%
frameratems = expt(iexp).frame_rate/1000;
cyctimems = cycTime/frameratems;
stimon_frames = input.nFramesOn;
stimoff_frames = input.nFramesOff;
stimontime = stimon_frames/frameratems;
stimofftime = stimoff_frames/frameratems;
trialresult = unique(trialOutcome);
responsedelay_frames = 3;

vs1R = [];
vsendR = [];
for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    vs1R{icyc} = bsxfun(@minus, squeeze(mean(tempdata(prepress_frames+responsedelay_frames+1:prepress_frames+responsedelay_frames+stimon_frames,:,:),1)), squeeze(mean(tempdata(prepress_frames:prepress_frames+responsedelay_frames,:,:),1)));
    vsendR{icyc} = bsxfun(@minus, squeeze(mean(tempdata(prepress_frames+(cycTime*(icyc-1))+responsedelay_frames+1:prepress_frames+(cycTime*(icyc-1))+responsedelay_frames+stimon_frames,:,:),1)), squeeze(mean(tempdata(prepress_frames+(cycTime*(icyc-1)):prepress_frames+(cycTime*(icyc-1))+responsedelay_frames,:,:),1)));
    
end
dataIntegral = cellfun(@squeeze, cellfun(@trapz,cycDataDFoverF_cmlvNoTarget,'UniformOutput', 0),'UniformOutput', 0);
%% resp integral scatter

figure(integralScatterFig);
hold on
tempdata = dataIntegral{CYC};
V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
% V_cycInd = cycV_ind{CYC};
% AV_cycInd = cycAV_ind{CYC};
scatter(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2),50,[col '.']);
% errorbarxy(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {['.' col], col, col});
hold on

meanIntegral_V{iexp} = mean(mean(tempdata(cells,V_cycInd),2),1);
meanIntegral_AV{iexp} = mean(mean(tempdata(cells,AV_cycInd),2),1);

errorbarxy(meanIntegral_V{iexp},meanIntegral_AV{iexp}, std(unique(reshape(tempdata(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,V_cycInd),1,[])))), std(unique(reshape(tempdata(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on


str1 = ['\leftarrow ' num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];
text(meanIntegral_V{iexp},meanIntegral_AV{iexp},str1)

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on

xlim([-3 7])
ylim([-3 7])
% xlim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
% ylim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

mouse_str = cellfun(@num2str,{av.mouse},'UniformOutput',0);
title([{['summary of mice - ' strjoin(mouse_str)] ; 'response integral, aud vs vis'; ['success trials' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'}])


%% resp 2 1st cyc scatter

figure(respfirststimScatterFig);
hold on
tempdata = vs1R{1};
V_cycInd = intersect(cycV_ind{1},find(strcmp(cycTrialOutcome{1},'success')));
AV_cycInd = intersect(cycAV_ind{1},find(strcmp(cycTrialOutcome{1},'success')));

scatter(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2),50,[col '.']);
% errorbarxy(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {['.' col], col, col});
hold on

meanVS1R_V{iexp} = mean(mean(tempdata(cells,V_cycInd),2),1);
meanVS1R_AV{iexp} = mean(mean(tempdata(cells,AV_cycInd),2),1);

errorbarxy(meanVS1R_V{iexp},meanVS1R_AV{iexp}, std(unique(reshape(tempdata(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,V_cycInd),1,[])))), std(unique(reshape(tempdata(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on

xlim([-0.02 0.04])
ylim([-0.02 0.04])
% xlim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
% ylim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

mouse_str = cellfun(@num2str,{av.mouse},'UniformOutput',0);
title([{['summary of mice - ' strjoin(mouse_str)] ; 'respo to first stim, aud vs vis'; ['success trials']; 'gr = vis, blk = aud'}])

%% resp 2 last cyc scatter

figure(resplaststimScatterFig);
hold on
tempdata = vsendR{CYC};
V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));

scatter(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2),50,[col '.']);
% errorbarxy(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {['.' col], col, col});
hold on

meanVSendR_V{iexp} = mean(mean(tempdata(cells,V_cycInd),2),1);
meanVSendR_AV{iexp} = mean(mean(tempdata(cells,AV_cycInd),2),1);

errorbarxy(meanVSendR_V{iexp},meanVSendR_AV{iexp}, std(unique(reshape(tempdata(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,V_cycInd),1,[])))), std(unique(reshape(tempdata(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on

xlim([-0.02 0.04])
ylim([-0.02 0.04])
% xlim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
% ylim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')

mouse_str = cellfun(@num2str,{av.mouse},'UniformOutput',0);
title([{['summary of mice - ' strjoin(mouse_str)] ; 'respo to last stim, aud vs vis'; ['success trials' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'}])

%% resp to last stim std
figure(stdScatterFig);

tempdata = vsendR{CYC};
V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
% V_cycInd = cycV_ind{CYC};
% AV_cycInd = cycAV_ind{CYC};
stdV = std(tempdata(cells,V_cycInd),[],2);
stdAV = std(tempdata(cells,AV_cycInd),[],2);
scatter(stdV,stdAV,50,[col '.'])
hold on
errorbarxy(mean(stdV), mean(stdAV),std(stdV)/sqrt(length(stdV)),std(stdAV)/sqrt(length(stdAV)),{['o' colS], colS, colS})
% errorbarxy(std(tempdata(cells,V_cycInd),[],2),std(tempdata(cells,AV_cycInd),[],2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {'ob', 'b', 'b'});
hold on
plot([-1:0.01:1],[-1:0.01:1],'k-')
hold on
xlim([0 0.1]);
ylim([0 0.1]);
axis square
xlabel('Vis Tr Resp (dF/F)')
ylabel('Aud Tr Resp (dF/F)')
mouse_str = cellfun(@num2str,{av.mouse},'UniformOutput',0);
title([{['summary of mice - ' strjoin(mouse_str)] ; 'std, aud vs vis'; ['success trials' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'}])

end


for ianalysis = 1:4
    figNameBase = analysisName{ianalysis};
    figure(eval(genvarname(figName{ianalysis})))
if ~isempty(fnouttemp)
    print([fnouttemp '\'  date '-' mouse dateofexp '_' figNameBase '_audVSvis_cellssummary'], '-dpdf');
else
    print([fnout '\'  date '_' figNameBase '_audVSvis'], '-dpdf');
end
    
end



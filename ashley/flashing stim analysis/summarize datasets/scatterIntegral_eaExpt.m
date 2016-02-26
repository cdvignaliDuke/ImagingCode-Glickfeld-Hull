close all
clear all
% iexp = 1;
fnouttemp = []; %
awFSAVdatasets
cellsetname = 'nCells';

av = behavParamsAV;
analysisName = {'PreTarget_respIntegral_Cells'; 'PreTarget_respAmp2basestim1_Cells'; 'PreTarget_respAmp2basestimlast_Cells';'PreTarget_std_Cells'};
figName = {'integralScatterFig';'respfirststimScatterFig';'resplaststimScatterFig';'stdScatterFig'};

% CYC = 6;
trialLengthMs_gt = 2500;

integralScatterFig_Vstd = figure;
integralScatterFig_AVstd = figure;
% respfirststimScatterFig = figure;
% resplaststimScatterFig = figure;
% stdScatterFig = figure;

subpC = 3;
if size(expt,2) >= subpC;
    subpR = ceil(size(expt,2)/subpC);
else
    subpR = 1;
end
mouse_str = unique(cellfun(@num2str,{expt.SubNum},'UniformOutput',0));


%%
iexp = 4;
% for iexp = 1:size(expt,2)
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

fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' dateofexp '\' runstr];

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

%% resp integral 

CYC = ceil(trialLengthMs_gt/(stimontime+stimofftime));
trialLengthFrames = cycTime*CYC;


data = cycDataDFoverF_cmlvNoTarget{CYC};

tempdata = dataIntegral{CYC};

% find integral of first and last half of response
dataInt = squeeze(trapz(data(prepress_frames:end,:,:)));
dataInt_first = squeeze(trapz(data(prepress_frames:floor(trialLengthFrames/2)+prepress_frames,:,:)));
dataInt_last = squeeze(trapz(data(floor(trialLengthFrames/2)+1:end,:,:)));

V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));

cellsInt_V = mean(dataInt(cells,V_cycInd),2);
cellsInt_AV = mean(dataInt(cells,AV_cycInd),2);

cellsIntFirst_V = mean(dataInt_first(cells,V_cycInd),2);
cellsIntFirst_AV = mean(dataInt_first(cells,AV_cycInd),2);

cellsIntLast_V = mean(dataInt_last(cells,V_cycInd),2);
cellsIntLast_AV = mean(dataInt_last(cells,AV_cycInd),2);

meanInt_V = mean(mean(dataInt(cells,V_cycInd),2),1);
meanInt_AV = mean(mean(dataInt(cells,AV_cycInd),2),1);

stdInt_V = std(dataInt(cells,V_cycInd),[],2);
stdInt_AV = std(dataInt(cells,AV_cycInd),[],2);

CoEV_V = stdInt_V/meanInt_V;

%% scatter
figure(integralScatterFig_Vstd);
suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'response integral, aud vs vis (Vste)'; ['success trials-' num2str(CYC) ' cycles']; [';n= ' num2str(length(cells)) '; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});

%plot whole integral response aud vs vis trials
subplot(1,3,1)
% scatter(cellsInt_V,cellsInt_AV,50,stdInt_V/(sqrt(length(V_cycInd))),'.');
scatter(cellsInt_V,cellsInt_AV,50,stdInt_V,'.');
hold on
errorbarxy(meanInt_V,meanInt_AV, std(unique(reshape(dataInt(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(dataInt(cells,V_cycInd),1,[])))), std(unique(reshape(dataInt(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(dataInt(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on
scatter(meanInt_V,meanInt_AV,'ko','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
ylim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V)]);
title('whole trial')

%first half of trial
subplot(1,3,2)
scatter(cellsIntFirst_V,cellsIntFirst_AV,50,std(dataInt_first(cells,V_cycInd),[],2),'.');
hold on
errorbarxy(mean(cellsIntFirst_V),mean(cellsIntFirst_AV), std(cellsIntFirst_V)/(sqrt(length(cellsIntFirst_V))),std(cellsIntFirst_AV)/(sqrt(length(cellsIntFirst_AV))), {['o' colS], colS, colS});
hold on
scatter(mean(cellsIntFirst_V),mean(cellsIntFirst_AV),'ro','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
ylim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V)]);
title('first half')

%second half of trial
subplot(1,3,3)
scatter(cellsIntLast_V,cellsIntLast_AV,50,std(dataInt_last(cells,V_cycInd),[],2),'.');
hold on
errorbarxy(mean(cellsIntLast_V),mean(cellsIntLast_AV), std(cellsIntLast_V)/(sqrt(length(cellsIntLast_V))),std(cellsIntLast_AV)/(sqrt(length(cellsIntLast_AV))), {['o' colS], colS, colS});
hold on
scatter(mean(cellsIntLast_V),mean(cellsIntLast_AV),'ro','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
ylim([floor(min(mean(dataInt_last,2)))-(floor(min(mean(dataInt_last,2)))/2) ceil(max(mean(dataInt_last,2)))+(ceil(max(mean(dataInt_last,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V)]);
title('second half')


%%

figure(integralScatterFig_Vstd);
if iexp == 1
    suptitle({[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; ['summary of mice - ' strjoin(mouse_str)] ; 'response integral, aud vs vis (Vstd)'; ['success trials-' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'})
end
hold on
subplot(subpR,subpC,iexp)

% scatter(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2),50,[col '.']);
scatter(cellsInt_V{iexp},cellsInt_AV{iexp},50,stdInt_V{iexp},'.');
hold on
% colorbar;
% caxis([0 10])
% % errorbarxy(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {['.' col], col, col});
% hold on


errorbarxy(meanInt_V{iexp},meanInt_AV{iexp}, std(unique(reshape(tempdata(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,V_cycInd),1,[])))), std(unique(reshape(tempdata(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on


% str1 = ['\leftarrow ' num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];
% text(meanIntegral_V{iexp},meanIntegral_AV{iexp},str1)

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on

xlim([-10 10])
ylim([-10 10])
% xlim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
% ylim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 10]);
title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];[expt(iexp).SubNum ';n= ' num2str(length(cells)) ';' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']})

figure(integralScatterFig_AVstd);
hold on
if iexp == 1
    suptitle({[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; ['summary of mice - ' strjoin(mouse_str)] ; 'response integral, aud vs vis (AVstd)'; ['success trials-' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'})
end

tempdata = dataIntegral{CYC};
V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));

cellsInt_V{iexp} = mean(tempdata(cells,V_cycInd),2);
cellsInt_AV{iexp} = mean(tempdata(cells,AV_cycInd),2);

meanInt_V{iexp} = mean(mean(tempdata(cells,V_cycInd),2),1);
meanInt_AV{iexp} = mean(mean(tempdata(cells,AV_cycInd),2),1);

stdInt_V{iexp} = std(tempdata(cells,V_cycInd),[],2);
stdInt_AV{iexp} = std(tempdata(cells,AV_cycInd),[],2);

subplot(subpR,subpC,iexp)

% scatter(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2),50,[col '.']);
scatter(cellsInt_V{iexp},cellsInt_AV{iexp},50,stdInt_AV{iexp},'.');
hold on
% colorbar;
% caxis([0 10]);
% % errorbarxy(mean(tempdata(cells,V_cycInd),2),mean(tempdata(cells,AV_cycInd),2), std(tempdata(cells,V_cycInd),[],2)/sqrt(size(tempdata(cells,V_cycInd),2)), std(tempdata(cells,AV_cycInd),[],2)/sqrt(size(tempdata(cells,AV_cycInd),2)), {['.' col], col, col});
% hold on


errorbarxy(meanInt_V{iexp},meanInt_AV{iexp}, std(unique(reshape(tempdata(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,V_cycInd),1,[])))), std(unique(reshape(tempdata(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(tempdata(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on


% str1 = ['\leftarrow ' num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];
% text(meanIntegral_V{iexp},meanIntegral_AV{iexp},str1)

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on

xlim([-10 10])
ylim([-10 10])
% xlim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
% ylim([floor(min(mean(tempdata,2)))-(floor(min(mean(tempdata,2)))/2) ceil(max(mean(tempdata,2)))+(ceil(max(mean(tempdata,2)))/2)]);
axis square

xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 10]);
title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];[expt(iexp).SubNum ';n= ' num2str(length(cells)) ';' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']})

% end

figure(integralScatterFig_Vstd);
% hold on
% suptitle({[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; ['summary of mice - ' strjoin(mouse_str)] ; 'response integral, aud vs vis (Vstd)'; ['success trials-' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'})
% colorbar;
% caxis([0 10]);

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

% print([fnout '\'  date '_PreTarget_respIntegral_audVSvis_Vstd_' cellsetname], '-dpdf');

figure(integralScatterFig_AVstd);
% hold on
% suptitle({[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; ['summary of mice - ' strjoin(mouse_str)] ; 'response integral, aud vs vis (AVstd)'; ['success trials-' num2str(CYC) ' cycles']; 'gr = vis, blk = aud'})
% colorbar;
% caxis([0 10]);
% print([fnout '\'  date '_PreTarget_respIntegral_audVSvis_AVstd_' cellsetname], '-dpdf');

%% modulation index V - AV/max(all)
maxInt = max([cell2mat(cellfun(@max,cellsInt_V,'UniformOutput',false)) cell2mat(cellfun(@max,cellsInt_AV,'UniformOutput',false))]);
mi = cellfun(@(x) x/maxInt,cellfun(@minus,cellsInt_V,cellsInt_AV,'UniformOutput',false),'UniformOutput',false);
figure;
for iexp = 1:size(expt,2)
    subplot(subpR,subpC,iexp)
    hist(mi{iexp},100)
    title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];[expt(iexp).SubNum ';n= ' num2str(length(cells)) ';' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']})
end

suptitle({[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; ['summary of mice - ' strjoin(mouse_str)] ; 'response integral, modulation index'})

% print([fnout '\'  date '_PreTarget_respIntegral_MIhist_' cellsetname], '-dpdf');
%% save fig

save(fullfile(fnout,'respIntegral.mat'),'cellsIntegral_V','cellsIntegral_AV','stdIntegral_V','stdIntegral_AV')
% 
% for ianalysis = 1:4
%     figNameBase = analysisName{ianalysis};
%     figure(eval(genvarname(figName{ianalysis})))
% if ~isempty(fnouttemp)
%     print([fnouttemp '\'  date '-' mouse dateofexp '_' figNameBase '_' num2str(CYC) 'cyc_audVSvis_cellssummary'], '-dpdf');
% else
%     print([fnout '\'  date '_' figNameBase '_audVSvis'], '-dpdf');
% end
%     
% end



close all
clear all
% iexp = 1;
fnouttemp = []; %
dataGroup = 'awFSAVdatasets';
eval(dataGroup);
cellsetname = 'nCells';

av = behavParamsAV;

% CYC = 6;
trialLengthMs_gt = 2500;

% integralScatterFig_Vstd = figure;
% integralScatterFig_AVstd = figure;


% subpC = 3;
% if size(expt,2) >= subpC;
%     subpR = ceil(size(expt,2)/subpC);
% else
%     subpR = 1;
% end
% mouse_str = unique(cellfun(@num2str,{expt.SubNum},'UniformOutput',0));


%%
% iexp = 4;
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

%% plot response trace for random cells, avg all cells
avgAllCellsTrace = figure;
randCellsAvgTrace = figure;
randCellsAllTrials = figure;
dataResp_V = mean(data(:,:,V_cycInd),3);
dataResp_AV = mean(data(:,:,AV_cycInd),3);
dataSte_V = std(data(:,:,V_cycInd),[],3);
dataSte_AV = std(data(:,:,AV_cycInd),[],3);

randCellSet = randperm(length(cells),16);

figure(randCellsAvgTrace)
suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'average response of rand cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
for i = 1:16
    subplot(4,4,i)
    plotCell = cells(randCellSet(i));
    shadedErrorBar(xAxisMs,dataResp_V(:,plotCell),dataSte_V(:,plotCell),'g',1);
    hold on
    shadedErrorBar(xAxisMs,dataResp_AV(:,plotCell),dataSte_AV(:,plotCell),'k',1);
    hold on
    vline(0,'k')
    for icyc = 1:CYC
        vline(cyctimems*icyc,'k:')
        hold on
    end
    title(['cell# ' num2str(plotCell)]);
end

figure(randCellsAllTrials)
suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'all trials of rand cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
for i = 1:16
    subplot(4,4,i)
    plotCell = cells(randCellSet(i));
    plot(xAxisMs,squeeze(data(:,plotCell,V_cycInd)),'g');
    hold on
    plot(xAxisMs,squeeze(data(:,plotCell,AV_cycInd)),'k');
    hold on
    vline(0,'k')
    hold on
    for icyc = 1:CYC
        vline(cyctimems*icyc,'k:')
        hold on
    title(['cell# ' num2str(plotCell)]);
    end

end

%% responsive cells
dataPercentile = 0.7;
dataVals = sort(data(:));
taskRespCutoff = dataVals(ceil(length(dataVals)*dataPercentile));


baselineRespCells = find(any(dataResp_V(1:prepress_frames,:) > taskRespCutoff,1));
stimRespCells_V = find(any(dataResp_V(prepress_frames+1:size(data,1),:) > taskRespCutoff,1));
stimRespCells_AV = find(any(dataResp_AV(prepress_frames+1:size(data,1),:) > taskRespCutoff,1));
drivenCellsInd_V = setdiff(stimRespCells_V,baselineRespCells);
drivenCellsInd_AV = setdiff(stimRespCells_AV,baselineRespCells);
drivenCellsInd = unique([drivenCellsInd_V drivenCellsInd_AV]);

%% plot traces of responsive cells
respCellsTrace = figure;
respCellsAvgTrace = figure;
respCellsAllTrials = figure;

if length(drivenCellsInd) > 16
    randCellSet = drivenCellsInd(sort(randperm(length(drivenCellsInd),16)));
    
    figure(respCellsAvgTrace)
    suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'average response of driven cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
    for i = 1:16
        subplot(4,4,i)
        plotCell = cells(randCellSet(i));
        shadedErrorBar(xAxisMs,dataResp_V(:,plotCell),dataSte_V(:,plotCell),'g',1);
        hold on
        shadedErrorBar(xAxisMs,dataResp_AV(:,plotCell),dataSte_AV(:,plotCell),'k',1);
        hold on
        vline(0,'k')
        for icyc = 1:CYC
            vline(cyctimems*icyc,'k:')
            hold on
        end
        title(['cell# ' num2str(plotCell)]);
    end
    
    figure(respCellsAllTrials)
    suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'all trials of driven cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
    for i = 1:16
        subplot(4,4,i)
        plotCell = cells(randCellSet(i));
        plot(xAxisMs,squeeze(data(:,plotCell,V_cycInd)),'g');
        hold on
        plot(xAxisMs,squeeze(data(:,plotCell,AV_cycInd)),'k');
        hold on
        vline(0,'k')
        hold on
        for icyc = 1:CYC
            vline(cyctimems*icyc,'k:')
            hold on
        end
        title(['cell# ' num2str(plotCell)]);
    end
    
else

    figure(respCellsAvgTrace)
    suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'average response of driven cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
    for i = 1:length(drivenCellsInd)
        subplot(4,4,i)
        plotCell = cells(drivenCellsInd(i));
        shadedErrorBar(xAxisMs,dataResp_V(:,plotCell),dataSte_V(:,plotCell),'g',1);
        hold on
        shadedErrorBar(xAxisMs,dataResp_AV(:,plotCell),dataSte_AV(:,plotCell),'k',1);
        hold on
        vline(0,'k')
        for icyc = 1:CYC
            vline(cyctimems*icyc,'k:')
            hold on
        end
        title(['cell# ' num2str(plotCell)]);
        
    end
    
    figure(respCellsAllTrials)
    suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'all trials of driven cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});
    for i = 1:length(drivenCellsInd)
        subplot(4,4,i)
        plotCell = cells(drivenCellsInd(i));
        plot(xAxisMs,squeeze(data(:,plotCell,V_cycInd)),'g');
        hold on
        plot(xAxisMs,squeeze(data(:,plotCell,AV_cycInd)),'k');
        hold on
        vline(0,'k')
        hold on
        for icyc = 1:CYC
            vline(cyctimems*icyc,'k:')
            hold on
        end
        title(['cell# ' num2str(plotCell)]);
    end
end

%% plot average trace
avgAllCellsTrace = figure;
avgRespCellsTrace = figure;

meanAll_V = mean(mean(data(:,:,V_cycInd),3),2);
steAll_V = std(mean(data(:,:,V_cycInd),3),[],2)/sqrt(length(V_cycInd));
meanAll_AV = mean(mean(data(:,:,AV_cycInd),3),2);
steAll_AV = std(mean(data(:,:,AV_cycInd),3),[],2)/sqrt(length(AV_cycInd));

meanResp_V = mean(mean(data(:,drivenCellsInd,V_cycInd),3),2);
steResp_V = std(mean(data(:,drivenCellsInd,V_cycInd),3),[],2)/sqrt(length(V_cycInd));
meanResp_AV = mean(mean(data(:,drivenCellsInd,AV_cycInd),3),2);
steResp_AV = std(mean(data(:,drivenCellsInd,AV_cycInd),3),[],2)/sqrt(length(AV_cycInd));


figure(avgAllCellsTrace)
shadedErrorBar(xAxisMs,meanAll_V,steAll_V,'g',1);
hold on
shadedErrorBar(xAxisMs,meanAll_AV,steAll_AV,'k',1);
hold on
vline(0,'k')
hold on
for icyc = 1:CYC
    vline(cyctimems*icyc,'k:')
    hold on
end
title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'all trials'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});

figure(avgRespCellsTrace)
shadedErrorBar(xAxisMs,meanResp_V,steResp_V,'g',1);
hold on
shadedErrorBar(xAxisMs,meanResp_AV,steResp_AV,'k',1);
hold on
vline(0,'k')
hold on
for icyc = 1:CYC
    vline(cyctimems*icyc,'k:')
    hold on
end
title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'all trials - responsive cells'; ['success trials-' num2str(CYC) ' cycles; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});

%% resp integral 

% find integral of first and last half of response
dataInt = squeeze(trapz(data(prepress_frames:end,:,:)));
dataInt_first = squeeze(trapz(data(prepress_frames:floor(trialLengthFrames/2)+prepress_frames,:,:)));
dataInt_last = squeeze(trapz(data(floor(trialLengthFrames/2)+1:end,:,:)));

cellsInt_V{iexp} = mean(dataInt(cells,V_cycInd),2);
cellsInt_AV{iexp} = mean(dataInt(cells,AV_cycInd),2);

cellsIntFirst_V{iexp} = mean(dataInt_first(cells,V_cycInd),2);
cellsIntFirst_AV{iexp} = mean(dataInt_first(cells,AV_cycInd),2);

cellsIntLast_V{iexp} = mean(dataInt_last(cells,V_cycInd),2);
cellsIntLast_AV{iexp} = mean(dataInt_last(cells,AV_cycInd),2);

meanInt_V{iexp} = mean(mean(dataInt(cells,V_cycInd),2),1);
meanInt_AV{iexp} = mean(mean(dataInt(cells,AV_cycInd),2),1);

stdInt_V{iexp} = std(dataInt(cells,V_cycInd),[],2);
stdInt_AV{iexp} = std(dataInt(cells,AV_cycInd),[],2);

CoEV_V{iexp} = stdInt_V{iexp}/meanInt_V{iexp};

%% plot scatters of aud vis vis trials
integralScatterFig_Vstd = figure;
% integralScatterFig_AVstd = figure;
figure(integralScatterFig_Vstd);
suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'response integral, aud vs vis (Vstdev)'; ['success trials-' num2str(CYC) ' cycles']; [';n= ' num2str(length(cells)) '; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});

%plot whole integral response aud vs vis trials
subplot(1,3,1)

scatter(cellsInt_V{iexp},cellsInt_AV{iexp},50,stdInt_V{iexp},'.');
hold on
errorbarxy(meanInt_V{iexp},meanInt_AV{iexp}, std(unique(reshape(dataInt(cells,V_cycInd),1,[])))/sqrt(length(unique(reshape(dataInt(cells,V_cycInd),1,[])))), std(unique(reshape(dataInt(cells,AV_cycInd),1,[])))/sqrt(length(unique(reshape(dataInt(cells,AV_cycInd),1,[])))), {['o' colS], colS, colS});
hold on
scatter(meanInt_V{iexp},meanInt_AV{iexp},'ro','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V{iexp})]);
title('whole trial')

%first half of trial
subplot(1,3,2)
scatter(cellsIntFirst_V{iexp},cellsIntFirst_AV{iexp},50,std(dataInt_first(cells,V_cycInd),[],2),'.');
hold on
errorbarxy(mean(cellsIntFirst_V{iexp}),mean(cellsIntFirst_AV{iexp}), std(cellsIntFirst_V{iexp})/(sqrt(length(cellsIntFirst_V{iexp}))),std(cellsIntFirst_AV{iexp})/(sqrt(length(cellsIntFirst_AV{iexp}))), {['o' colS], colS, colS});
hold on
scatter(mean(cellsIntFirst_V{iexp}),mean(cellsIntFirst_AV{iexp}),'ro','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V{iexp})]);
title('first half')

%second half of trial
subplot(1,3,3)
scatter(cellsIntLast_V{iexp},cellsIntLast_AV{iexp},50,std(dataInt_last(cells,V_cycInd),[],2),'.');
hold on
errorbarxy(mean(cellsIntLast_V{iexp}),mean(cellsIntLast_AV{iexp}), std(cellsIntLast_V{iexp})/(sqrt(length(cellsIntLast_V{iexp}))),std(cellsIntLast_AV{iexp})/(sqrt(length(cellsIntLast_AV{iexp}))), {['o' colS], colS, colS});
hold on
scatter(mean(cellsIntLast_V{iexp}),mean(cellsIntLast_AV{iexp}),'ro','filled')
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
% xlim([-10 10])
% ylim([-10 10])
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
colorbar;
caxis([0 max(stdInt_V{iexp})]);
title('second half')

%% plot integral scatter, labelling each orientation pref a different color
integralScatterFig_oriPref = figure;

oriColor = {'g';'b';'r';'k'};
figure(integralScatterFig_oriPref);

suptitle({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date) '; ' expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'response integral, aud vs vis (Vstdev)'; ['success trials-' num2str(CYC) ' cycles']; ['n= ' num2str(length(cells)) '; ' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']});

% whole trial
subplot(1,3,1)
for i = 1:length(Oris)
    oriCells = intersect(cells,cellsPref{i});
    oriLegend{i} = num2str(Oris(i));
    cellsScatter(i) = scatter(cellsInt_V{iexp}(oriCells,:),cellsInt_AV{iexp}(oriCells,:),50,'.',oriColor{i});
    hold on
    scatter(mean(cellsInt_V{iexp}(oriCells,:),1),mean(cellsInt_AV{iexp}(oriCells,:),1),[oriColor{i} 'o'],'filled')
    hold on
    steX = std(cellsInt_V{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    steY = std(cellsInt_AV{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    errorbarxy(mean(cellsInt_V{iexp}(oriCells,:),1),mean(cellsInt_AV{iexp}(oriCells,:),1),steX, steY,{oriColor{i} oriColor{i} oriColor{i}})
    hold on
end
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(cellsScatter,oriLegend,'Location','southeast')
title('whole trial')
% first half
subplot(1,3,2)
for i = 1:length(Oris)
    oriCells = intersect(cells,cellsPref{i});
    oriLegend{i} = num2str(Oris(i));
    cellsScatter(i) = scatter(cellsIntFirst_V{iexp}(oriCells,:),cellsIntFirst_AV{iexp}(oriCells,:),50,'.',oriColor{i});
    hold on
    scatter(mean(cellsIntFirst_V{iexp}(oriCells,:),1),mean(cellsIntFirst_AV{iexp}(oriCells,:),1),[oriColor{i} 'o'],'filled')
    hold on
    steX = std(cellsIntFirst_V{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    steY = std(cellsIntFirst_AV{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    errorbarxy(mean(cellsIntFirst_V{iexp}(oriCells,:),1),mean(cellsIntFirst_AV{iexp}(oriCells,:),1),steX, steY,{oriColor{i} oriColor{i} oriColor{i}})
    hold on
end
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(cellsScatter,oriLegend,'Location','southeast')
title('first half')
% second half
subplot(1,3,3)
for i = 1:length(Oris)
    oriCells = intersect(cells,cellsPref{i});
    oriLegend{i} = num2str(Oris(i));
    cellsScatter(i) = scatter(cellsIntLast_V{iexp}(oriCells,:),cellsIntLast_AV{iexp}(oriCells,:),50,'.',oriColor{i});
    hold on
    scatter(mean(cellsIntLast_V{iexp}(oriCells,:),1),mean(cellsIntLast_AV{iexp}(oriCells,:),1),[oriColor{i} 'o'],'filled')
    hold on
    steX = std(cellsIntLast_V{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    steY = std(cellsIntLast_AV{iexp}(oriCells,:),[],1)/sqrt(length(oriCells));
    errorbarxy(mean(cellsIntLast_V{iexp}(oriCells,:),1),mean(cellsIntLast_AV{iexp}(oriCells,:),1),steX, steY,{oriColor{i} oriColor{i} oriColor{i}})
    hold on
end

plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
ylim([floor(min(mean(dataInt,2)))-(floor(min(mean(dataInt,2)))/2) ceil(max(mean(dataInt,2)))+(ceil(max(mean(dataInt,2)))/2)]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(cellsScatter,oriLegend,'Location','southeast')
title('second half')

%% calc modulation index for each cell, plot histogram
maxInt = max([cell2mat(cellfun(@max,cellsInt_V,'UniformOutput',false)) cell2mat(cellfun(@max,cellsInt_AV,'UniformOutput',false))]);
mi = cellfun(@(x) x/maxInt,cellfun(@minus,cellsInt_V,cellsInt_AV,'UniformOutput',false),'UniformOutput',false);
MIhist = figure;
hist(mi{iexp},100)
title({[num2str(expt(iexp).SubNum) '-' num2str(expt(iexp).date)];[expt(1).img_loc{1} '-' expt(1).img_strct{1}] ; 'response integral, modulation index';['n= ' num2str(length(cells)) ';' num2str(length(V_cycInd)) 'vis tr, ' num2str(length(AV_cycInd)) 'aud tr']})


%% save figures and data
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

figure(integralScatterFig_Vstd);
print([fnout '\'  exptTag 'PreTarget_respIntegral_audVSvis_Vstd_' cellsetname], '-dpdf');
figure(integralScatterFig_oriPref);
print([fnout '\'  exptTag 'PreTarget_respIntegral_audVSvis_oriPref_' cellsetname], '-dpdf');
figure(MIhist);
print([fnout '\'  exptTag 'PreTarget_respIntegral_MI_hist_' cellsetname], '-dpdf');
figure(respCellsAvgTrace);
print([fnout '\'  exptTag 'PreTarget_drivenCellsAvg_' cellsetname], '-dpdf');
figure(respCellsAllTrials);
print([fnout '\'  exptTag 'PreTarget_drivenCellsTrials_' cellsetname], '-dpdf');
figure(avgAllCellsTrace);
print([fnout '\'  exptTag 'PreTarget_meanAllTrace_' cellsetname], '-dpdf');
figure(avgRespCellsTrace);
print([fnout '\'  exptTag 'PreTarget_meanResponsiveCellsTrace_' cellsetname], '-dpdf');
end
%% summarize integral scatter across mice
% clear expt
eval(dataGroup);
meanInt_V = cell2mat(cellfun(@mean,cellsInt_V,'UniformOutput',false));
steInt_V = cell2mat(cellfun(@std,cellsInt_V,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsInt_V,'UniformOutput',false)));
meanInt_AV = cell2mat(cellfun(@mean,cellsInt_AV,'UniformOutput',false));
steInt_AV = cell2mat(cellfun(@std,cellsInt_AV,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsInt_V,'UniformOutput',false)));

meanIntFirst_V = cell2mat(cellfun(@mean,cellsIntFirst_V,'UniformOutput',false));
steIntFirst_V = cell2mat(cellfun(@std,cellsIntFirst_V,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsIntFirst_V,'UniformOutput',false)));
meanIntFirst_AV = cell2mat(cellfun(@mean,cellsIntFirst_AV,'UniformOutput',false));
steIntFirst_AV = cell2mat(cellfun(@std,cellsIntFirst_AV,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsIntFirst_AV,'UniformOutput',false)));

meanIntLast_V = cell2mat(cellfun(@mean,cellsIntLast_V,'UniformOutput',false));
steIntLast_V = cell2mat(cellfun(@std,cellsIntLast_V,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsIntLast_V,'UniformOutput',false)));
meanIntLast_AV = cell2mat(cellfun(@mean,cellsIntLast_AV,'UniformOutput',false));
steIntLast_AV = cell2mat(cellfun(@std,cellsIntLast_AV,'UniformOutput',false))./sqrt(cell2mat(cellfun(@length,cellsIntLast_AV,'UniformOutput',false)));


for iexp = 1:size(expt,2)
   imouse = find(strcmp(cellfun(@num2str,{av.mouse},'UniformOutput',0), expt(iexp).SubNum));
   cStr{iexp} = av(imouse).col_str;
   mouseMat(iexp) = av(imouse).mouse;
end

integralSummaryFig = figure;
suptitle(dataGroup)
mice = unique(mouseMat);

subplot(1,3,1)
for iMs = 1:length(mice)
    imouse = find(cell2mat({av.mouse}) == mice(iMs));
    msInd = find(mouseMat == mice(iMs));
    msMean_V = mean(meanInt_V(msInd));
    msMean_AV = mean(meanInt_AV(msInd));
    msSte_V = std(meanInt_V(msInd))/sqrt(length(meanInt_V(msInd)));
    msSte_AV = std(meanInt_AV(msInd))/sqrt(length(meanInt_AV(msInd)));
    
    colorCell = {[av(imouse).col_str 'o'],av(imouse).col_str,av(imouse).col_str};
    errorbarxy(meanInt_V(msInd),meanInt_AV(msInd),steInt_V(msInd),steInt_AV(msInd),colorCell)
    hold on
    errorbarxy(msMean_V,msMean_AV,msSte_V,msSte_AV,colorCell)
    hold on
    msColMean(iMs) = scatter(msMean_V,msMean_AV,[av(imouse).col_str 'o'],'filled');
    summaryLegend{iMs} = num2str(mice(iMs));
    
end
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
ylim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(msColMean,summaryLegend,'Location','southeast')
title('whole trial')

subplot(1,3,2)
for iMs = 1:length(mice)
    imouse = find(cell2mat({av.mouse}) == mice(iMs));
    msInd = find(mouseMat == mice(iMs));
    msMean_V = mean(meanIntFirst_V(msInd));
    msMean_AV = mean(meanIntFirst_AV(msInd));
    msSte_V = std(meanIntFirst_V(msInd))/sqrt(length(meanIntFirst_V(msInd)));
    msSte_AV = std(meanIntFirst_AV(msInd))/sqrt(length(meanIntFirst_AV(msInd)));
    
    colorCell = {[av(imouse).col_str 'o'],av(imouse).col_str,av(imouse).col_str};
    errorbarxy(meanIntFirst_V(msInd),meanIntFirst_AV(msInd),steIntFirst_V(msInd),steIntFirst_AV(msInd),colorCell)
    hold on
    errorbarxy(msMean_V,msMean_AV,msSte_V,msSte_AV,colorCell)
    hold on
    msColMean(iMs) = scatter(msMean_V,msMean_AV,[av(imouse).col_str 'o'],'filled');
    summaryLegend{iMs} = num2str(mice(iMs));
    
end
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
ylim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(msColMean,summaryLegend,'Location','southeast')
title('first half')

subplot(1,3,3)
for iMs = 1:length(mice)
    imouse = find(cell2mat({av.mouse}) == mice(iMs));
    msInd = find(mouseMat == mice(iMs));
    msMean_V = mean(meanIntLast_V(msInd));
    msMean_AV = mean(meanIntLast_AV(msInd));
    msSte_V = std(meanIntLast_V(msInd))/sqrt(length(meanIntLast_V(msInd)));
    msSte_AV = std(meanIntLast_AV(msInd))/sqrt(length(meanIntLast_AV(msInd)));
    
    colorCell = {[av(imouse).col_str 'o'],av(imouse).col_str,av(imouse).col_str};
    errorbarxy(meanIntLast_V(msInd),meanIntLast_AV(msInd),steIntLast_V(msInd),steIntLast_AV(msInd),colorCell)
    hold on
    errorbarxy(msMean_V,msMean_AV,msSte_V,msSte_AV,colorCell)
    hold on
    msColMean(iMs) = scatter(msMean_V,msMean_AV,[av(imouse).col_str 'o'],'filled');
    summaryLegend{iMs} = num2str(mice(iMs));
    
end
plot([-10:0.1:20],[-10:0.1:20],'k--')
hold on
xlim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
ylim([floor(min([meanInt_V meanInt_AV])) ceil(max([meanInt_V meanInt_AV]))]);
axis square
xlabel('Vis Tr Resp')
ylabel('Aud Tr Resp')
legend(msColMean,summaryLegend,'Location','southeast')
title('second half')

%save fig
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

figure(integralSummaryFig);
print(fullfile('Z:\Analysis\FSAV Summaries\',dataGroup,'integralSummaryAcrossMice'), '-dpdf');

%% save integral variables only if all cells were analyzed
if length(eval(cellsetname)) == 1
    save(fullfile('Z:\Analysis\FSAV Summaries\',dataGroup,'cellsIntegral.mat'),'cellsInt_V','cellsInt_AV','cellsIntFirst_V','cellsIntFirst_AV','cellsIntLast_V','cellsIntLast_AV');
end

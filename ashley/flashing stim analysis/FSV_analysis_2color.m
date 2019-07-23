clear all
close all
ds = 'FSAV_V1_SOM';
cellsOrDendrites = 1;
doLoadPreviousAnalysis = false;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_V1_SOM')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','visual only','behavior',...
        [titleStr '_']); 
    bxExpt = true;
elseif strcmp(ds,'FSAV_V1_SOM_naive')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','visual only','naive',...
        [titleStr '_']);     
    bxExpt = false;
end


%%
rewwin = 32:39;
nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 7;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
firstWinFr = (3:44)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

% minCellN_SIFRmatch = 36;

tar2analyze = [22.5,90];

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data

    antiDataExpt = struct;
    antiDataExpt_passive = struct;
    tarDataExpt = struct;
    tarDataExpt_passive = struct;
    decodeDataExpt = struct;
%     oriTuningExpt = struct;
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end

            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            tarDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            tarDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            
            for itag = 1:2
                dd = d.tag(itag).av(visualTrials).align(alignStart);
                if ~isfield(dd,'respTC')
                    continue
                end
                maxCycles = max(dd.nCycles);
                cycTC = cell(1,maxCycles);
    %             longTC = [];
                if bxExpt
                    hits = strcmp(dd.outcome,'success');
                    misses = strcmp(dd.outcome,'ignore');
                else
                    hits = true(1,length(dd.outcome));
                    misses = false(1,length(dd.outcome));
                end
                for icyc = 1:maxCycles
                    tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
                    cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                    cycTC{icyc} = tc(...
                        (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                end
                longTC = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                    dd.nCycles >= nCycLongTC & (hits | misses));

                antiDataExpt(exptN).tag(itag).name = d.tag(itag).name;
                antiDataExpt(exptN).tag(itag).longTC = longTC;
                antiDataExpt(exptN).tag(itag).cycTC = cycTC;
                
                dd = d.tag(itag).av(visualTrials).align(alignTarget);
                easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(2));
                hardTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(1));
                tarDataExpt(exptN).tag(itag).targetTC = {hardTarTC,easyTarTC};
                
                nc = size(dd.respTC,2);
                rt = dd.reactTime;
                hits = strcmp(dd.outcome,'success');
                misses = strcmp(dd.outcome,'ignore');
                for itar = 1:2
                    rt(misses & dd.ori == tar2analyze(itar)) = ...
                        mean(rt(hits & dd.ori == tar2analyze(itar)));
                end
                rt_fr = round(rt/1000*frameRateHz);
                tarDataExpt(exptN).tag(itag).rewTC = cell(1,2);
                rewTC = nan(nBaselineFr+nFrames1s,nc,length(rt));
                for itr = 1:length(rt)
                    rewTC(:,:,itr) = dd.respTC(...
                        (1:(nBaselineFr+nFrames1s))+rt_fr(itr),:,itr);
                end
                rewTC_bl = rewTC - mean(rewTC(basewin_0,:,:),1);
                tarDataExpt(exptN).tag(itag).rewTC = {...
                    rewTC_bl(:,:,hits & dd.ori == tar2analyze(1)),...
                    rewTC_bl(:,:,hits & dd.ori == tar2analyze(2));...
                    rewTC_bl(:,:,misses & dd.ori == tar2analyze(1)),...
                    rewTC_bl(:,:,misses & dd.ori == tar2analyze(2))};
            end
%             if bxExpt
%                 d = mousePass(imouse).expt(iexp);
%                 antiDataExpt_passive(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
%                 antiDataExpt_passive(exptN).exptCycLengthFr = cycLengthFr;
%                 tarDataExpt_passive(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
%                 tarDataExpt_passive(exptN).exptCycLengthFr = cycLengthFr;
%                 for itag = 1:2
%                     cycTC = cell(1,maxCycles);
%         %             longTC = [];
%                     dd = d.tag(itag).av(visualTrials).align(alignStart);
%                     hits = true(1,length(dd.outcome));
%                     misses = false(1,length(dd.outcome));
%                     for icyc = 1:maxCycles
%                         tc = dd.respTC(:,:,dd.nCycles >= icyc & (hits | misses));
%                         cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
%                         cycTC{icyc} = tc(...
%                             (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
%                     end
%                     longTC = dd.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
%                         dd.nCycles >= nCycLongTC & (hits | misses));
% 
%                     antiDataExpt_passive(exptN).tag(itag).name = d.tag(itag).name;
%                     antiDataExpt_passive(exptN).tag(itag).longTC = longTC;
%                     antiDataExpt_passive(exptN).tag(itag).cycTC = cycTC;
% 
%                     dd = d.tag(itag).av(visualTrials).align(alignTarget);
%                     easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori > 60);
%                     tarDataExpt_passive(exptN).tag(itag).easyTarTC = easyTarTC;
%                 end
%             end
            
            
            decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            if isfield(d.tag(1).av(visualTrials).align(alignStart),'respTC') && ...
                    isfield(d.tag(2).av(visualTrials).align(alignStart),'respTC')
                dd = d.tag(1).av(visualTrials).align(alignFA);
                tc = dd.respTC;
                [~,nc,nt] = size(tc);
                trOut = repmat({'fa'},[1,nt]);
                trStim = repmat(0,[1,nt]);
                trResp = squeeze(mean(tc(respwin,:,:),1) - mean(tc(basewin_0,:,:),1));
                dd = d.tag(1).av(visualTrials).align(alignCR);
                tc = dd.respTC;
                nt = size(tc,3);
                trOut = cat(2,trOut,repmat({'cr'},[1,nt]));
                trStim = cat(2,trStim,repmat(0,[1,nt]));
                trResp = cat(2,trResp,...
                    squeeze(mean(tc(respwin,:,:),1) - mean(tc(basewin_0,:,:),1)));
                dd = d.tag(1).av(visualTrials).align(alignTarget);
                tc = dd.respTC;
                n = size(tc,3);
                trOut = cat(2,trOut,dd.outcome);
                trOut(strcmp(trOut,'success')) = {'h'};
                trOut(strcmp(trOut,'ignore')) = {'m'};
                trStim = cat(2,trStim,dd.ori);
                trResp = cat(2,trResp,...
                    squeeze(mean(tc(respwin_target,:,:),1) - mean(tc(basewin_0_target,:,:),1)));
                cellTagID = false(1,nc);
%                 cellInd = respCellsExpt(exptN).tag(1).lateCycRespCells | ...
%                     respCellsExpt(exptN).tag(1).lateCycRespCells;
                tc = cat(3,d.tag(2).av(visualTrials).align(alignFA).respTC,...
                    d.tag(2).av(visualTrials).align(alignCR).respTC,...
                    d.tag(2).av(visualTrials).align(alignTarget).respTC);
                [~,nc,nt] = size(tc);
                trResp = cat(1,trResp,...
                    squeeze(mean(tc(respwin,:,:),1) - mean(tc(basewin_0,:,:),1)));
                cellTagID = cat(2,cellTagID,true(1,nc));
%                 cellInd = cat(1,cellInd,respCellsExpt(exptN).tag(2).lateCycRespCells | ...
%                     respCellsExpt(exptN).tag(2).lateCycRespCells);                
                
                decodeDataExpt(exptN).resp = trResp;
                decodeDataExpt(exptN).trOut = trOut;
                decodeDataExpt(exptN).trStim = trStim;
                decodeDataExpt(exptN).cellTag = cellTagID;
%                 decodeDataExpt(exptN).cellInd = cellInd';
            end
            
        end
    end

    shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
%     isShortCycExpt = [];
    respCellsExpt = struct;
    antiAnalysis = struct;
    antiAnalysis_passive = struct;
    tarAnalysis = struct;
    tarAnalysis_passive = struct;
    cellInfo = struct;
    for itag = 1:2
        tagNames = cellfun(@(y) y.name,...
            cellfun(@(x) x(itag),{antiDataExpt.tag},'unif',0),'unif',0);
        ind = cellfun(@(x) ~isempty(x),tagNames);
        antiAnalysis.tag(itag).name = unique(tagNames(ind));
        for iexp = 1:nexp
            if isempty(antiDataExpt(iexp).tag(itag).name)
                continue
            end
            antiAnalysis.tag(itag).name = antiDataExpt(iexp).tag(itag).name;
            firstTC = antiDataExpt(iexp).tag(itag).cycTC{1};
            longTC = antiDataExpt(iexp).tag(itag).longTC;
            tarTC = tarDataExpt(iexp).tag(itag).targetTC;

            firstRespCells = ttest(...
                squeeze(mean(firstTC(respwin,:,:),1)),...
                squeeze(mean(firstTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateWinRespCells = ttest(...
                squeeze(mean(longTC(lateWinFr,:,:),1)),...
                squeeze(mean(longTC(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha);
            lateWinSuppCells = ttest(...
                squeeze(mean(longTC(lateWinFr,:,:),1)),...
                squeeze(mean(longTC(basewin_0,:,:),1)),...
                'dim',2,'tail','left','alpha',cellGroupsAlpha);
            tarRespCells = sum(cell2mat(cellfun(@(x) ttest(...
                squeeze(mean(x(respwin,:,:),1)),...
                squeeze(mean(x(basewin_0,:,:),1)),...
                'dim',2,'tail','right','alpha',cellGroupsAlpha),tarTC,'unif',0)),2)>0;

            respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
            respCellsExpt(iexp).tag(itag).firstRespCells = firstRespCells;
            respCellsExpt(iexp).tag(itag).lateWinRespCells = lateWinRespCells;
            respCellsExpt(iexp).tag(itag).lateWinSuppCells = lateWinSuppCells;
            respCellsExpt(iexp).tag(itag).targetRespCells = tarRespCells;
        end
        
            antiAnalysis.tag(itag).longTC = [];
            antiAnalysis.tag(itag).longTCErr = [];
            antiAnalysis.tag(itag).cycTC = cell(1,nCycles);
            antiAnalysis.tag(itag).cycTCErr = cell(1,nCycles);
            antiAnalysis.tag(itag).lateCycTC = [];
            antiAnalysis.tag(itag).lateCycTCErr = [];
            for iexp = 1:nexp
                if isempty(antiDataExpt(iexp).tag(itag).name)
                    continue
                end
                longTC = antiDataExpt(iexp).tag(itag).longTC;
                cycTC = antiDataExpt(iexp).tag(itag).cycTC(1:nCycles);
                lateCycTC = [];
                for icyc = 5:length(cycTC)
                    lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
                end

                antiAnalysis.tag(itag).longTC = cat(2,antiAnalysis.tag(itag).longTC,...
                    mean(longTC,3));
                antiAnalysis.tag(itag).longTCErr = cat(2,antiAnalysis.tag(itag).longTCErr,...
                    ste(longTC,3));
                antiAnalysis.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.tag(itag).cycTC,...
                    cycTC,'unif',0);
                antiAnalysis.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.tag(itag).cycTCErr,...
                    cycTC,'unif',0);
                antiAnalysis.tag(itag).lateCycTC = cat(2,antiAnalysis.tag(itag).lateCycTC,...
                    mean(lateCycTC,3));
                antiAnalysis.tag(itag).lateCycTCErr = cat(2,antiAnalysis.tag(itag).lateCycTCErr,...
                    ste(lateCycTC,3));

                lateCycRespCells = ttest(...
                    squeeze(mean(lateCycTC(respwin,:,:),1)),...
                    squeeze(mean(lateCycTC(basewin_0,:,:),1)),...
                    'dim',2,'tail','right','alpha',cellGroupsAlpha);

                respCellsExpt(iexp).tag(itag).lateCycRespCells = lateCycRespCells;
            end
        
        if bxExpt
            antiAnalysis_passive.tag(itag).longTC = [];
            antiAnalysis_passive.tag(itag).longTCErr = [];
            antiAnalysis_passive.tag(itag).cycTC = cell(1,nCycles);
            antiAnalysis_passive.tag(itag).cycTCErr = cell(1,nCycles);
            antiAnalysis_passive.tag(itag).lateCycTC = [];
            antiAnalysis_passive.tag(itag).lateCycTCErr = [];
            for iexp = 1:nexp
            longTC = antiDataExpt_passive(iexp).tag(itag).longTC;
            cycTC = antiDataExpt_passive(iexp).tag(itag).cycTC(1:nCycles);
            lateCycTC = [];
            for icyc = 5:length(cycTC)
                lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
            end

            antiAnalysis_passive.tag(itag).longTC = cat(2,antiAnalysis_passive.tag(itag).longTC,...
                mean(longTC,3));
            antiAnalysis_passive.tag(itag).longTCErr = cat(2,antiAnalysis_passive.tag(itag).longTCErr,...
                ste(longTC,3));
            antiAnalysis_passive.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis_passive.tag(itag).cycTC,...
                cycTC,'unif',0);
            antiAnalysis_passive.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis_passive.tag(itag).cycTCErr,...
                cycTC,'unif',0);
            antiAnalysis_passive.tag(itag).lateCycTC = cat(2,antiAnalysis_passive.tag(itag).lateCycTC,...
                mean(lateCycTC,3));
            antiAnalysis_passive.tag(itag).lateCycTCErr = cat(2,antiAnalysis_passive.tag(itag).lateCycTCErr,...
                ste(lateCycTC,3));
            end
        end
        
        tarAnalysis.tag(itag).tc = cell(1,2);
        tarAnalysis.tag(itag).rewTC = cell(2,2);
        tarAnalysis_passive.tag(itag).easyTC = [];
        for iexp = 1:nexp
            if isempty(antiDataExpt(iexp).tag(itag).name)
                continue
            end
            htc = tarDataExpt(iexp).tag(itag).targetTC{1};
            etc = tarDataExpt(iexp).tag(itag).targetTC{2};
            tarAnalysis.tag(itag).tc{1} = cat(2,tarAnalysis.tag(itag).tc{1},...
                mean(htc,3));
            tarAnalysis.tag(itag).tc{2} = cat(2,tarAnalysis.tag(itag).tc{2},...
                mean(etc,3));
            htc = tarDataExpt(iexp).tag(itag).rewTC(:,1);
            etc = tarDataExpt(iexp).tag(itag).rewTC(:,2);
            tarAnalysis.tag(itag).rewTC{1,1} = cat(2,tarAnalysis.tag(itag).rewTC{1,1},...
                mean(htc{1},3));
            tarAnalysis.tag(itag).rewTC{1,2} = cat(2,tarAnalysis.tag(itag).rewTC{1,2},...
                mean(htc{2},3));
            tarAnalysis.tag(itag).rewTC{2,1} = cat(2,tarAnalysis.tag(itag).rewTC{2,1},...
                mean(etc{1},3));
            tarAnalysis.tag(itag).rewTC{2,2} = cat(2,tarAnalysis.tag(itag).rewTC{2,2},...
                mean(etc{2},3));
            if bxExpt
                tc = tarDataExpt_passive(iexp).tag(itag).easyTarTC;
                tarAnalysis_passive.tag(itag).easyTC = cat(2,tarAnalysis_passive.tag(itag).easyTC,...
                    mean(tc,3));
            end
        end
        
        cellInfo.tag(itag).name = antiDataExpt(1).tag(itag).name;
        cellInfo.tag(itag).isFlex = [];
        cellInfo.tag(itag).firstRespCells = [];
        cellInfo.tag(itag).lateWinRespCells = [];
        cellInfo.tag(itag).lateWinSuppCells = [];
        cellInfo.tag(itag).targetRespCells = [];
        cellInfo.tag(itag).lateCycRespCells = [];
        cellInfo.tag(itag).aurocReward = [];
        for iexp = 1:nexp
            if isempty(antiDataExpt(iexp).tag(itag).name)
                continue
            end
            if itag == 2 && isempty(antiDataExpt(iexp).tag(1).name)
                isFlex = true(length(respCellsExpt(iexp).tag(itag).firstRespCells),1);
            elseif itag == 2
                isFlex = false(length(respCellsExpt(iexp).tag(itag).firstRespCells),1);                
            else
                isFlex = [];
            end
            
            htc = tarDataExpt(iexp).tag(itag).rewTC(:,1);    
            nc = size(htc{1},2);
            aurocrew = nan(nc,1);
            for ic = 1:nc
                aurocrew(ic) = roc_gh(squeeze(mean(htc{2}(rewwin,ic,:))),...
                    squeeze(mean(htc{1}(rewwin,ic,:))));
            end
            
            cellInfo.tag(itag).isFlex = cat(1,cellInfo.tag(itag).isFlex,isFlex);
            cellInfo.tag(itag).firstRespCells = cat(1,cellInfo.tag(itag).firstRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).firstRespCells}')));
            cellInfo.tag(itag).lateWinRespCells = cat(1,cellInfo.tag(itag).lateWinRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinRespCells}')));
            cellInfo.tag(itag).lateWinSuppCells = cat(1,cellInfo.tag(itag).lateWinSuppCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinSuppCells}')));
            cellInfo.tag(itag).targetRespCells = cat(1,cellInfo.tag(itag).targetRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).targetRespCells}')));
            cellInfo.tag(itag).lateCycRespCells = cat(1,cellInfo.tag(itag).lateCycRespCells,...
                logical(cell2mat({respCellsExpt(iexp).tag(itag).lateCycRespCells}')));
            cellInfo.tag(itag).aurocReward = cat(1,cellInfo.tag(itag).aurocReward,aurocrew);
        end
        cellInfo.tag(itag).minRespCells = (mean(antiAnalysis.tag(itag).cycTC{1}(respwin,:),1) > ...
            minRespThreshold)';
    end
%     cellInfo.tag(itag).isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
%     cellInfo.tag(itag).oriResp = cell2mat({oriTuningExpt.oriResp}');
%     cellInfo.tag(itag).oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
%     cellInfo.tag(itag).oriFit = cell2mat({oriTuningExpt.fit})';
%     cellInfo.tag(itag).oriPref = cell2mat({oriTuningExpt.oriPref})';
%     cellInfo.tag(itag).hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';

for iexp = 1:nexp
        cellInd = respCellsExpt(iexp).tag(1).lateCycRespCells | ...
            respCellsExpt(iexp).tag(1).lateCycRespCells;
        cellInd = cat(1,cellInd,respCellsExpt(iexp).tag(2).lateCycRespCells | ...
            respCellsExpt(iexp).tag(2).lateCycRespCells); 
        decodeDataExpt(iexp).cellInd = cellInd';
end

save([fnout 'adaptAnalysis'],'antiAnalysis','tarAnalysis','antiAnalysis_passive','cellInfo','respCellsExpt')
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.01 0.12];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];

tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.tag(1).longTC,1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_rew = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = size(antiAnalysis.tag(itag).cycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 4.2];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.8 0.8];
siLimSum = [-0.5 2.5];

% lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);

%%
tarAdaptFig = figure;
for itag = 1:2
    ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells;
    cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),antiAnalysis.tag(itag).cycTC,'unif',0);

    adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);

    setFigParams4Print('landscape')
    figure
    if bxExpt
        suptitle(sprintf('behav, %s, Min. & First Stim Resp. Neurons',...
            cellInfo.tag(itag).name))
    else
        suptitle(sprintf('naive, %s, Min. & First Stim Resp. Neurons',...
            cellInfo.tag(itag).name))
    end
    for icyc = 1:nCycles
        subplot(2,nCycles,icyc)
        bl = mean(antiAnalysis.tag(itag).cycTC{icyc}(basewin_0,ind),1);
        y = antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind) - bl;
        yerr = ste(antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind),2);
        shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
        figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
        figYAxis([],'dF/F',cycTCLim)
        figAxForm
        hline(0,'k:')
        hold on
        vline(respWinTT,'k--');
        if icyc == 1
            title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
                num2str(length(ind))))
        else
            title(sprintf('Stim #%s',num2str(icyc)))
        end
    end

    subplot 223
    y = cellfun(@(x) mean(x(ind)),cycResp);
    yerr = cellfun(@(x) ste(x(ind),2),cycResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    title(sprintf('First Stim Resp. Cells (%s/%s)',num2str(sum(ind)),num2str(length(ind))))
    subplot 224
    y = cellfun(@(x) mean(x(ind)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'Norm dF/F',[0 1.5])
    figAxForm
    title(sprintf('All Resp. Cells, n=%s', num2str(sum(ind))))

    print([fnout 'adaptation_anti_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
    
    ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells & ...
        cellInfo.tag(itag).targetRespCells;
    tarResp = cellfun(@(x) ...
        mean(x(respwin_target,:),1) - mean(x(basewin_0_target,:),1),...
        cat(2, antiAnalysis.tag(itag).cycTC{1}, ...
        tarAnalysis.tag(itag).tc),'unif',0);

    adaptResp = cellfun(@(x) x./tarResp{1},tarResp,'unif',0);
    
    figure(tarAdaptFig)
    subplot(1,2,itag)
    y = cellfun(@(x) mean(x(ind)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
    errorbar(1:3,y,yerr,'.')
    figXAxis([],'Stim (deg)',[0 4],1:3,[0,tar2analyze])
    figYAxis([],'Norm dF/F',[0 4])
    figAxForm
    title(sprintf('Target & First Stim. Resp. %s Cells, n=%s',...
        antiAnalysis.tag(itag).name,num2str(sum(ind))))
    print([fnout 'adaptation_target_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
end

%%

hmLim = [-0.3 0.3];
figure
suptitle('Anticipation: Any Task-Modulated Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
%     ind = logical(cellInfo.tag(itag).lateWinRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
    ind = cellInfo.tag(itag).lateWinRespCells | ...
        cellInfo.tag(itag).lateWinSuppCells | ...
        cellInfo.tag(itag).firstRespCells | ...
        cellInfo.tag(itag).lateCycRespCells | ...
        cellInfo.tag(itag).targetRespCells;
    lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
    lateWinResp = mean(lateWinTC(lateWinFr,:),1);
    [~,lateWinSortInd] = sort(lateWinResp);
    hm = flipud(lateWinTC(:,lateWinSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
    if itag == 2
        subplot(2,2,3)
        ind = cellInfo.tag(itag).isFlex & (cellInfo.tag(itag).lateWinRespCells | ...
            cellInfo.tag(itag).lateWinSuppCells | ...
            cellInfo.tag(itag).firstRespCells | ...
            cellInfo.tag(itag).lateCycRespCells);
        lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
        lateWinResp = mean(lateWinTC(lateWinFr,:),1);
        [~,lateWinSortInd] = sort(lateWinResp);
        hm = flipud(lateWinTC(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('%s Cells,flex-GCaMP6s',antiAnalysis.tag(itag).name))
        subplot(2,2,4)
        ind = ~cellInfo.tag(itag).isFlex & (cellInfo.tag(itag).lateWinRespCells | ...
            cellInfo.tag(itag).lateWinSuppCells | ...
            cellInfo.tag(itag).firstRespCells | ...
            cellInfo.tag(itag).lateCycRespCells);
        lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
        lateWinResp = mean(lateWinTC(lateWinFr,:),1);
        [~,lateWinSortInd] = sort(lateWinResp);
        hm = flipud(lateWinTC(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('%s Cells,flex-tdTom',antiAnalysis.tag(itag).name))
    end
end

print([fnout '_anti_heatmapsAllCells'],'-dpdf','-fillpage')


figure
suptitle('Target Resp Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    for it = 1:2
        if it == 1
            offset = 0;
        else
            offset = 2;
        end
        subplot(2,2,itag+offset)
        ind = logical(cellInfo.tag(itag).targetRespCells);
    %     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
    %     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
        tc = tarAnalysis.tag(itag).tc{it}(:,ind);
        if it == 1
            resp = mean(tc(respwin,:),1);
            [~,lateWinSortInd] = sort(resp);
        else
            tc1 = tarAnalysis.tag(itag).tc{it}(:,ind);
            resp = mean(tc1(respwin,:),1);
            [~,lateWinSortInd] = sort(resp);
        end
            
        hm = flipud(tc(:,lateWinSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        imagesc(hm(:,:))
        hold on
        figXAxis([],'Time from Target (ms)',[],ttLabelFr_target,ttLabel_target)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('%s Cells, Target=%s',antiAnalysis.tag(itag).name,...
            num2str(tar2analyze(it))))
    end
end

print([fnout '_tar_heatmapsTarCells'],'-dpdf','-fillpage')

%% decoding analysis
decodeAnalysis = struct;
for iexp = 1:nexp
    rng(0)
    if sum(decodeDataExpt(iexp).cellTag) == 0
        continue
    end
    cellInd_notag = decodeDataExpt(iexp).cellInd & ~decodeDataExpt(iexp).cellTag;
    cellInd_tag = decodeDataExpt(iexp).cellInd & decodeDataExpt(iexp).cellTag;
    cellInd = false(length(cellInd_notag),1);
    if sum(cellInd_notag) > maxCellN
        ind = find(cellInd_notag);
        cellSampleID = randsample(ind,maxCellN);
        cellInd(cellSampleID) = true;
    else
        cellInd = cellInd_notag;
    end
    maxTagCellN = round(0.2*sum(cellInd));
    if sum(cellInd_tag) > maxTagCellN
        ind = find(cellInd_tag);
        cellSampleID = randsample(ind,maxTagCellN);
        cellInd(cellSampleID) = true;
    else
        cellInd(cellInd_tag) = true;
    end
    fprintf('Expt %s: %s\n',num2str(iexp),decodeDataExpt(iexp).exptName)
    decodeAnalysis(iexp).exptName = decodeDataExpt(iexp).exptName;
    decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
    decodeAnalysis(iexp).cellInd = cellInd;
    decodeAnalysis(iexp).cellTag = decodeDataExpt(iexp).cellTag;
    
    respAllCells = zscore(decodeDataExpt(iexp).resp)';
    trOut = decodeDataExpt(iexp).trOut;
    
    fprintf('Select Trials...\n')
    trStimID = discretize(decodeDataExpt(iexp).trStim,oriBins);
    nStimPerBin = histcounts(trStimID);
    minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
    respStimSort = cell(1,nStimBins);
    trOutStimSort = cell(1,nStimBins);
    for istim = 1:nStimBins
        ind = find(trStimID == istim);
        if length(ind) >= minTrN_mdl
            if istim == 1
                matchTrialsInd = [];
                if sum(nStimPerBin >= minTrN_mdl) == 2
                    n = minBinN;
                elseif minBinN == nStimPerBin(istim)
                    error('not enough FA/CR trials')
                else
                    n = (nStimBins-1).*minBinN;
                    if n > length(ind)
                        error('not enough FA/CR trials')
                    end
                end
                indSample = randsample(ind,n);
                matchTrialsInd = cat(2,matchTrialsInd,indSample);
            else
                indSample = randsample(ind,minBinN);
                matchTrialsInd = cat(2,matchTrialsInd,indSample);
            end
            respStimSort{istim} = respAllCells(indSample,cellInd);
            trOutStimSort{istim} = trOut(indSample);
        end
    end
    nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
    for istim = 1:nStimBins
        if istim == 1
            stimSortInd = cell(1,nStimBins);
            stimSortInd{istim} = 1:nMatchedTrials;
        else
            stimSortInd{istim} = ...
                (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
        end
    end
    
    resp = respAllCells(matchTrialsInd,cellInd);
    [detectTrInd, targetTrInd] = getStimAndBehaviorYs(...
        decodeDataExpt(iexp).trOut(matchTrialsInd));
    
    fprintf('Run Model...\n')
    C = eye(size(resp,2));
    p=1;
    [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
    [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');

    dv_detect = mean(detectTrInd);
    dv_target = mean(targetTrInd);

    pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
    pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);

    pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
    pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
    
    pctCorrDetect_xStim_train = nan(1,nStimBins);
    pctCorrDetect_xStim_ho = nan(1,nStimBins);
    pctCorrTarget_xStim_train = nan(1,nStimBins);
    pctCorrTarget_xStim_ho = nan(1,nStimBins);
    for istim = 1:nStimBins
        if nStimPerBin(istim) >= minTrN_mdl
            [detectStimInd, targetStimInd] = getStimAndBehaviorYs(...
                trOutStimSort{istim});
            pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
                detectGLM,respStimSort{istim},detectStimInd,dv_detect);
            pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                resp,detectTrInd,stimSortInd{istim},dv_detect);
            pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
                targetGLM,respStimSort{istim},targetStimInd,dv_target);
            pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
                resp,targetTrInd,stimSortInd{istim},dv_target);
        end
    end
    
    
    fprintf('Add Results to Struct...\n')
    decodeAnalysis(iexp).dvDetect = dv_detect;
    decodeAnalysis(iexp).dvTarget = dv_target;
    decodeAnalysis(iexp).detectWeight = detectGLM.beta(2:end);
    decodeAnalysis(iexp).targetWeight = targetGLM.beta(2:end);
    decodeAnalysis(iexp).pctCorrectAllDetect_train = pctCorrectDetect_train;
    decodeAnalysis(iexp).pctCorrectAllDetect_holdout = pctCorrectDetect_ho;
    decodeAnalysis(iexp).pctCorrectAllTarget_train = pctCorrectTarget_train;
    decodeAnalysis(iexp).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
    decodeAnalysis(iexp).pctCorrectXStimDetect_train = pctCorrDetect_xStim_train;
    decodeAnalysis(iexp).pctCorrectXStimDetect_holdout = pctCorrDetect_xStim_ho;
    decodeAnalysis(iexp).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
    decodeAnalysis(iexp).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;

end
ind = cellfun(@(x) ~isempty(x),{decodeAnalysis.exptName});
decodeAnalysis = decodeAnalysis(ind);

%%
nexp_dc = size(decodeAnalysis,2);
targetWeight = [];
detectWeight = [];
cellTag = [];
pctCorrect_target = nan(nexp_dc,4);
pctCorrect_detect = nan(nexp_dc,4);
for iexp = 1:nexp_dc
    targetWeight = cat(1,targetWeight,decodeAnalysis(iexp).targetWeight);
    detectWeight = cat(1,detectWeight,decodeAnalysis(iexp).detectWeight);
    ind = false(1,sum(decodeAnalysis(iexp).cellInd));
    cellTag = cat(2,cellTag,...
        decodeAnalysis(iexp).cellTag(decodeAnalysis(iexp).cellInd));
    
    pctCorrect_target(iexp,1:3) = decodeAnalysis(iexp).pctCorrectXStimTarget_holdout;
    pctCorrect_target(iexp,end) = decodeAnalysis(iexp).pctCorrectAllTarget_holdout;
    pctCorrect_detect(iexp,1:3) = decodeAnalysis(iexp).pctCorrectXStimDetect_holdout;
    pctCorrect_detect(iexp,end) = decodeAnalysis(iexp).pctCorrectAllDetect_holdout;
end

binEdges = -2:0.2:2;
figure
subplot 221
histogram(targetWeight(~cellTag),binEdges,'Normalization','probability')
hold on
histogram(targetWeight(cellTag==1),binEdges,'Normalization','probability')
title('Stimulus Model Weights')
figXAxis([],'Weight',[-2 2])
figYAxis([],'Fraction of Cells',[0 0.5])
h(1) = vline(mean(targetWeight(~cellTag)),'b');
h(2) = vline(mean(targetWeight(cellTag==1)),'r');
legend(h,{antiAnalysis.tag(1).name,antiAnalysis.tag(2).name},'location','northeast')
figAxForm
subplot 222
histogram(detectWeight(~cellTag),binEdges,'Normalization','probability')
hold on
histogram(detectWeight(cellTag==1),binEdges,'Normalization','probability')
title('Choice Model Weights')
figXAxis([],'Weight',[-2 2])
figYAxis([],'Fraction of Cells',[0 0.5])
h(1) = vline(mean(detectWeight(~cellTag)),'b');
h(2) = vline(mean(detectWeight(cellTag==1)),'r');
legend(h,{antiAnalysis.tag(1).name,antiAnalysis.tag(2).name},'location','northeast')
figAxForm
stimLabel = {'0';'HT';'ET';'All'};
subplot 223
y = pctCorrect_target;
yerr = ste(pctCorrect_target,1);
plot(1:4,y,'k-');
hold on
errorbar(1:4,mean(y),yerr,'.')
figXAxis([],'Stim',[0 5], 1:4, stimLabel)
figYAxis([],'Fraction Correct',[0 1])
hline(0.5,'k--')
figAxForm
title('Stimulus Model Performance')
subplot 224
y = pctCorrect_detect;
yerr = ste(pctCorrect_detect,1);
plot(1:4,y,'k-');
hold on
errorbar(1:4,mean(y),yerr,'.')
figXAxis([],'Stim',[0 5], 1:4, stimLabel)
figYAxis([],'Fraction Correct',[0 1])
hline(0.5,'k--')
figAxForm
title('Choice Model Performance')

print([fnout '_decodeModelPerf&Weights'],'-dpdf','-fillpage')

%% reward analysis
binEdges = 0:0.1:1;
hmLim = [-0.2 0.2];
tcLim = [-0.08 0.08];
rewHM = figure;
suptitle('Hard Target, Reward: Target Resp. Cells')
colormap(brewermap([],'*RdBu'));
rewAnal = figure;
suptitle('Hard Target, Reward: Target Resp. Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
%     ind = cellInfo.tag(itag).lateWinRespCells | ...
%         cellInfo.tag(itag).lateWinSuppCells | ...
%         cellInfo.tag(itag).firstRespCells | ...
%         cellInfo.tag(itag).lateCycRespCells | ...
%         cellInfo.tag(itag).targetRespCells;
    ind = cellInfo.tag(itag).targetRespCells==1;
    figure(rewHM)
    subplot(2,2,itag)
    rewTC = tarAnalysis.tag(itag).rewTC{1,1}(:,ind);
    rewResp = mean(rewTC(rewwin,:),1);
    [~,rewWinSortInd] = sort(rewResp);
    hm = flipud(rewTC(:,rewWinSortInd)');
    imagesc(hm)
    hold on
    figXAxis([],'Time from Reward (ms)',[],ttLabelFr_rew,ttLabel_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Hits, %s Cells',antiAnalysis.tag(itag).name))
    subplot(2,2,itag+2)
    rewTC = tarAnalysis.tag(itag).rewTC{1,2}(:,ind);
    hm = flipud(rewTC(:,rewWinSortInd)');
    imagesc(hm)
    hold on
    figXAxis([],'Time from Reward (ms)',[],ttLabelFr_rew,ttLabel_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Misses, %s Cells',antiAnalysis.tag(itag).name))
    
    figure(rewAnal)
    subplot(2,2,itag)
    y = tarAnalysis.tag(itag).rewTC{1,1}(:,ind);
    yerr = ste(y,2);
    h=shadedErrorBar_chooseColor(tt_targetTC,mean(y,2),yerr,[0 0 0]);
    hold on
    y = tarAnalysis.tag(itag).rewTC{1,2}(:,ind);
    yerr = ste(y,2);
    if ~all(isnan(yerr(:)))
        shadedErrorBar_chooseColor(tt_targetTC,mean(y,2),yerr,[.75 0 0]);
    end
    figXAxis([],'Time from Reward',[tt_targetTC(1) tt_targetTC(end)],ttLabel_target,ttLabel_target)
    figYAxis([],'dF/F',tcLim)
    hline(0,'k--')
    figAxForm
    title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
    
    auroc = cellInfo.tag(itag).aurocReward;
    subplot(2,2,3)
    hold on
%     if ~all(isnan(auroc(ind)))
%         histogram(auroc(ind),binEdges,'Normalization','probability');
%     end
%     figXAxis([],'Reward auROC (Miss vs. Hit)',[0 1])
%     figYAxis([],'Fraction of Cells,',[0 0.5])
%     figAxForm
    if ~all(isnan(auroc(ind)))
        cdfplot(auroc(ind));
    end
    figXAxis([],'Reward auROC (Miss vs. Hit)',[0 1])
    figYAxis([],'Fraction of Cells,',[0 1])
    figAxForm
    legend({antiAnalysis.tag(itag).name,antiAnalysis.tag(itag).name},...
        'location','northwest')    
end
figure(rewHM)
print([fnout '_reward_heatmaps'],'-dpdf','-fillpage')
figure(rewAnal)
print([fnout '_reward_TC'],'-dpdf','-fillpage')
%% cluster analysis
allCells = [];
allCells_pass = [];
cellTag = [];
for itag = 1:2
    ind = cellInfo.tag(itag).lateWinRespCells | ...
        cellInfo.tag(itag).lateWinSuppCells | ...
        cellInfo.tag(itag).firstRespCells | ...
        cellInfo.tag(itag).lateCycRespCells;
    lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
    lateWinResp = mean(lateWinTC(lateWinFr,:),1);
    [~,lateWinSortInd] = sort(lateWinResp);
    lateWinTC_pass = antiAnalysis_passive.tag(itag).longTC(:,ind);
    allCells = cat(2,allCells,lateWinTC(:,lateWinSortInd));
    allCells_pass = cat(2,allCells_pass,lateWinTC_pass(:,lateWinSortInd));
    cellTag = cat(2,cellTag,repmat({cellInfo.tag(itag).name},1,sum(ind)));
end

nc = size(allCells,2);
allCells_norm = allCells./max(allCells,[],1);
allCells_pass_norm = allCells_pass./max(allCells,[],1);

% allCells_pca = pca(allCells_norm);

nCluster = 3;
[clustID, clustCentroid] = kmeans(allCells_norm',nCluster,'MaxIter',1000000,'Replicates',1000);
[~,clustID_pass] = pdist2(clustCentroid,allCells_pass_norm','euclidean','smallest',1);

figure
subplot 221
histogram(clustID,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 222
histogram(clustID_pass,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 223
for ic = 1:nCluster
    cInd = clustID == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
tagInd = strcmp(cellTag,'SOM+');
plot(ones(1,sum(tagInd)).*(nCluster+1),find(tagInd),'k.','MarkerSize',10)
figXAxis([],'Cluster #',[0 (nCluster+2)],1:(nCluster+1),...
    cat(2,cellfun(@num2str,num2cell(1:nCluster),'unif',0),{'SOM+'}))
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
subplot 224
for ic = 1:nCluster
    cInd = clustID_pass == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
figXAxis([],'Cluster #',[0 (nCluster+1)])
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
print([fnout 'clusterAnalysisBxPass_clusterID'],'-dpdf','-fillpage')

clusterTC = cell(1,nCluster);
clusterTC_pass = cell(1,nCluster);
for ic = 1:nCluster
    cInd = clustID == ic;
    clusterTC{ic} = allCells(:,cInd);
    cInd = clustID_pass == ic;
    clusterTC_pass{ic} = allCells_pass(:,cInd);    
end

figure
for ic = 1:nCluster
    subplot(2,nCluster,ic)
    y = mean(clusterTC{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC{ic},2))))
    
    subplot(2,nCluster,ic+nCluster)
    y = mean(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC_pass{ic},2))))
end
print([fnout 'clusterAnalysisBxPass_clusterTC'],'-dpdf','-fillpage')

switchClustID = clustID'-clustID_pass ~= 0;

figure
subplot 121
y = sum(switchClustID & ~tagInd)./sum(~tagInd);
plot(1,y,'.','MarkerSize',20)
hold on
y = sum(switchClustID & tagInd)./sum(tagInd);
plot(2,y,'.','MarkerSize',20);
figXAxis([],'',[0 3],1:2,{'SOM-','SOM+'})
figYAxis([],'Fraction Switch Cluster',[0 0.2])
figAxForm
subplot 122
hold on
for ic = 1:nCluster
    cInd = clustID' == ic & switchClustID;
    n = histcounts(clustID_pass(cInd),1:nCluster+1);
    scatter(ic.*(ones(1,nCluster)),1:nCluster,(n*5)+1)
end
figXAxis([],'Behavior Cluster #',[0 nCluster+1])
figYAxis([],'Passive Cluster #',[0 nCluster+1])
figAxForm
title(sprintf('Cluster Switch Cells, n=%s',num2str(sum(switchClustID))))
print([fnout 'clusterAnalysisBxPass_nClusterSwitch'],'-dpdf','-fillpage')

%%
    figure
    suptitle('Target Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).targetRespCells);
    tarResp_bx = mean(tarAnalysis.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis.tag(itag).easyTC(basewin_0,:),1);
    tarResp_pass = mean(tarAnalysis_passive.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis_passive.tag(itag).easyTC(basewin_0,:),1);

    tarRespLim = [-0.01 0.12];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = tarResp_bx(ind);
        y = tarResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(tarRespLim,tarRespLim,'k--')
        figXAxis([],'Behavior',tarRespLim)
        figYAxis([],'Passive',tarRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end

figure
    suptitle('Late Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).lateCycRespCells);
    lateResp_bx = mean(antiAnalysis.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis.tag(itag).lateCycTC(basewin_0,:),1);
    lateResp_pass = mean(antiAnalysis_passive.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis_passive.tag(itag).lateCycTC(basewin_0,:),1);

    lateRespLim = [-0.01 0.1];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = lateResp_bx(ind);
        y = lateResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(lateRespLim,lateRespLim,'k--')
        figXAxis([],'Behavior',lateRespLim)
        figYAxis([],'Passive',lateRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end
clear all
close all
ds = 'FSV_V1';
cellsOrDendrites = 1;
doLoadPreviousAnalysis = false;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
% titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));

fnout = fullfile(rc.caOutputDir,ds,[ds '_']);

%%
% rewwin = 32:39;
% nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 8;
lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
respwin_opt = respwin;
respwin_target_opt = respwin_target;
% firstWinFr = (3:44)+nBaselineFr;
% minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;
% 
% oriBinSize = 45;
% orientations = 0:oriBinSize:(180-oriBinSize);
% oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
% nOri = length(orientations);
% 
% nMovWin = 15;
% movWinLabelFr = 30:(30+nMovWin-1);
% movWinLabelFr_target = 30:(30+nMovWin-1);
% % movWinLabelMs = 
% 
% % minCellN_SIFRmatch = 36;
% 
% tar2analyze = [22.5,90];
% 
% trOutType = {'h';'m';'fa';'cr'};
% trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
% tagName = {'SOM-','SOM+'};

avName = {'Visual';'Auditory'};
%% pool experiment data
    mice = struct;
    antiDataExpt = struct;
    tarDataExpt = struct;
    decodeDataExpt = struct;
%     antiDataExpt_passive = struct;
%     tarDataExpt_passive = struct;
    decodeDataExpt_passive = struct;
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
            indicator = expt(exptID).indicator{2};
            if contains(indicator,'f')
                respwin_opt = respwin;
                respwin_target_opt = respwin_target;
            else
                respwin_opt = respwin+1;
                respwin_target_opt = respwin_target+1;
            end
            stimSize = mouse(imouse).expt(iexp).info.fsavSize;
            bxExpt = logical(d.isBehav);
            if bxExpt
                if logical(mousePass(imouse).expt(iexp).passExptDone)
                    passExpt = true;
                else
                    passExpt = false;
                end
            else
                passExpt = false;
            end
            attnExpt = d.hasAttn;
            if attnExpt == 1
                attnExpt = true;
            else
                attnExpt = false;
            end
            if iexp == 1
                mice(imouse).name = mouse(imouse).mouse_name;
                mice(imouse).behav = bxExpt;
                mice(imouse).attention = attnExpt;
                mice(imouse).passive = passExpt;
                tagtypes = {d.tag.name};cellfun(@(x) ~isempty(x),{d.tag.name});
                mice(imouse).tagTypes = tagtypes(cellfun(@(x) ~isempty(x),{d.tag.name}));
            end
            
            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            antiDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            antiDataExpt(exptN).isBehav = bxExpt;
            antiDataExpt(exptN).hasAttn = attnExpt;
            antiDataExpt(exptN).passExpt = passExpt;
            antiDataExpt(exptN).indicator = indicator;
            antiDataExpt(exptN).stimSize = stimSize;
            antiDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            tarDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            tarDataExpt(exptN).exptCycLengthFr = cycLengthFr;
            
            if contains(antiDataExpt(exptN).indicator,'f')
                respwin_opt = respwin;
                respwin_target_opt = respwin_target;
            else
                respwin_opt = respwin+1;
                respwin_target_opt = respwin_target+1;
            end
            
            decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            for iav = 1:2
                decodeDataExpt(exptN).av(iav).name = avName{iav};
                for itag = 1:2
                    decodeDataExpt(exptN).av(iav).tag(itag).name = d.tag(itag).name;
                    dd = d.tag(itag).av(iav);
                    trOut = [];
                    trStim = [];
                    trResp = [];
                    if isempty(dd.align(1).respTC)
                        continue
                    end
                    for ialign = 2:4
                        if ialign==4
                            rw = respwin_target_opt;
                            bw = basewin_0_target;
                        else
                            rw = respwin_opt;
                            bw = basewin_0;
                        end
                            
                        dc = dd.align(ialign);
                        r = squeeze(mean(dc.respTC(rw,:,:),1)) - squeeze(mean(dc.respTC(bw,:,:),1));
                        trResp = cat(2,trResp,r);
                        if ialign == 2
                            trOut = cat(2,trOut,repmat({'fa'},[1,size(r,2)]));
                            trStim = cat(2,trStim,repmat(0,[1,size(r,2)]));
                        elseif ialign == 3
                            trOut = cat(2,trOut,repmat({'cr'},[1,size(r,2)]));
                            trStim = cat(2,trStim,repmat(0,[1,size(r,2)]));
                        else
                            outs = dc.outcome;
                            outs(strcmp(outs,'success')) = {'h'};
                            outs(strcmp(outs,'ignore')) = {'m'};
                            trOut = cat(2,trOut,outs);
                            if iav == 1
                                trStim = cat(2,trStim,dc.ori);
                            else
                                trStim = cat(2,trStim,dc.amp);
                            end
                        end
                    end
                    decodeDataExpt(exptN).av(iav).tag(itag).trResp = trResp;
                    decodeDataExpt(exptN).av(iav).tag(itag).trOut = trOut;
                    decodeDataExpt(exptN).av(iav).tag(itag).trStim = trStim;
                end
            end
            
            for iav = 1:2
                for itag = 1:2
                    
                    dd = d.tag(itag).av(iav).align(alignStart);
                    if isempty(dd.respTC)
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

                    antiDataExpt(exptN).av(iav).name = avName{iav};
                    antiDataExpt(exptN).av(iav).tag(itag).name = d.tag(itag).name;
                    antiDataExpt(exptN).av(iav).tag(itag).longTC = longTC;
                    antiDataExpt(exptN).av(iav).tag(itag).cycTC = cycTC;
                    
                    if passExpt
                        dp = mousePass(imouse).expt(iexp).tag(itag).av(iav).align(alignStart);
                        if isempty(dp.respTC)
                            continue
                        end
                        maxCycles = max(dp.nCycles);
                        cycTC = cell(1,maxCycles);
            %             longTC = [];
                        if bxExpt
                            hits = strcmp(dp.outcome,'success');
                            misses = strcmp(dp.outcome,'ignore');
                        else
                            hits = true(1,length(dp.outcome));
                            misses = false(1,length(dp.outcome));
                        end
                        for icyc = 1:maxCycles
                            tc = dp.respTC(:,:,dp.nCycles >= icyc & (hits | misses));
                            cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                            cycTC{icyc} = tc(...
                                (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                        end
                        longTC = dp.respTC(1:(longTrialLengthFr+nBaselineFr),:,...
                            dp.nCycles >= nCycLongTC & (hits | misses));
                        
                        antiDataExpt(exptN).av(iav).tag(itag).longTC_pass = longTC;
                        antiDataExpt(exptN).av(iav).tag(itag).cycTC_pass = cycTC;
                    end

                    dd = d.tag(itag).av(iav).align(alignTarget);
                    if iav == 1
                        targets = dd.ori;
                        binnedTarget = discretize(targets,oriBins);
                    elseif iav == 2
                        targets = dd.amp;
                        binnedTarget = discretize(targets,ampBins);
                    end
    %                 easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(2));
    %                 hardTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(1));
                    easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,binnedTarget == 3);
                    hardTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,binnedTarget == 2);
                    
                    tarDataExpt(exptN).av(iav).name = avName{iav};
                    tarDataExpt(exptN).av(iav).meanTarget = ...
                        [mean(targets(binnedTarget == 2)), mean(targets(binnedTarget == 3))];
                    tarDataExpt(exptN).av(iav).tag(itag).targetTC = {hardTarTC,easyTarTC};
                    
                    if passExpt
                        dd = mousePass(imouse).expt(iexp).tag(itag).av(iav).align(alignTarget);
                        if iav == 1
                            targets = dd.ori;
                            binnedTarget = discretize(targets,oriBins);
                        elseif iav == 2
                            targets = dd.amp;
                            binnedTarget = discretize(targets,ampBins);
                        end
        %                 easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(2));
        %                 hardTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,dd.ori == tar2analyze(1));
                        easyTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,binnedTarget == 3);
                        hardTarTC = dd.respTC(1:(nBaselineFr+nFrames1s),:,binnedTarget == 2);
                        tarDataExpt(exptN).av(iav).tag(itag).targetTC_pass = {hardTarTC,easyTarTC};
                    
                    end
                end
            end

        end
    end

%     shortCycExptInd = cell2mat({antiDataExpt.exptCycLengthFr}) == 11;
%     isShortCycExpt = [];

tagName = [];
indicatorName = [];
stimSize = [];
for iexp = 1:nexp
    tagName = cat(2,tagName,{antiDataExpt(iexp).av(1).tag(1).name},...
        {antiDataExpt(iexp).av(1).tag(2).name});
    indicatorName  = cat(2,indicatorName,{antiDataExpt(iexp).indicator});
    stimSize = double(unique(cat(2,stimSize,antiDataExpt(iexp).stimSize)));
end
tagName = unique(tagName(cellfun(@(x) ~isempty(x),tagName)));
indicatorName = unique(indicatorName);
ngcamp = length(indicatorName);
nsize = length(stimSize);
tagName_unique = tagName;
tagName_unique(strcmp(tagName,'ANY')|strcmp(tagName,'EMX+')) = {'pPYR'};
tagName_unique = unique(tagName_unique);
ntag = length(tagName_unique);
isBehav = cell2mat({antiDataExpt.isBehav});

cellInfo = struct;
antiAnalysis = struct;
decodeCellStruct = struct;
for iav = 1:2
    antiAnalysis.av(iav).name = avName{iav};
    for itag = 1:ntag
        antiAnalysis.av(iav).tag(itag).name = tagName_unique{itag};
        antiAnalysis.av(iav).tag(itag).longTC = [];
        antiAnalysis.av(iav).tag(itag).longTCErr = [];
        antiAnalysis.av(iav).tag(itag).cycTC = cell(1,nCycles);
        antiAnalysis.av(iav).tag(itag).cycTCErr = cell(1,nCycles);
        antiAnalysis.av(iav).tag(itag).lateCycTC = [];
        antiAnalysis.av(iav).tag(itag).lateCycTCErr = [];
        antiAnalysis.av(iav).tag(itag).longTC_pass = [];
        antiAnalysis.av(iav).tag(itag).cycTC_pass = cell(1,nCycles);
        antiAnalysis.av(iav).tag(itag).lateCycTC_pass = [];
        if iav == 1
            cellInfo.tag(itag).name = tagName_unique{itag};
            cellInfo.tag(itag).behav = [];
            cellInfo.tag(itag).attention = [];
            cellInfo.tag(itag).passive = [];
            cellInfo.tag(itag).indicatorID = [];
            cellInfo.tag(itag).stimSize = [];
            cellInfo.tag(itag).firstResp = [];
            cellInfo.tag(itag).lateCycResp = [];
            cellInfo.tag(itag).lateWinResp = [];
            cellInfo.tag(itag).lateWinSupp = [];
            cellInfo.tag(itag).minRespCells = [];
            cellInfo.tag(itag).targetResp = [];
            cellInfo.tag(itag).lateCycSI = [];
            cellInfo.tag(itag).lateCycSI_pass = [];
            cellInfo.tag(itag).aurocTD = [];
            cellInfo.tag(itag).aurocTD_test = [];
            cellInfo.tag(itag).aurocTD_pass = [];
            cellInfo.tag(itag).aurocTD_test_pass = [];
            decodeCellStruct(iexp).tag(itag).name = tagName_unique(itag);
        end
        for iexp = 1:nexp
            if contains(antiDataExpt(iexp).indicator,'f')
                respwin_opt = respwin;
                respwin_target_opt = respwin_target;
            else
                respwin_opt = respwin+1;
                respwin_target_opt = respwin_target+1;
            end
            for iexptag = 1:2
                if strcmp(tagName_unique(itag),'pPYR')
                    tagInd = false(1,ntag);
                    tagInd(itag) = strcmp(...
                        antiDataExpt(iexp).av(iav).tag(iexptag).name,'ANY')|...
                        strcmp(antiDataExpt(iexp).av(iav).tag(iexptag).name,'EMX+');
                else
                    tagInd = strcmp(tagName_unique,antiDataExpt(iexp).av(iav).tag(iexptag).name);
                end
                if tagInd(itag)
                    fprintf('%s, expt %s\n', tagName_unique{tagInd},num2str(iexp))
                    longTC = antiDataExpt(iexp).av(iav).tag(iexptag).longTC;
                    longTC_AV = cat(3,...
                        antiDataExpt(iexp).av(visualTrials).tag(iexptag).longTC,...
                        antiDataExpt(iexp).av(auditoryTrials).tag(iexptag).longTC);
                    cycTC = antiDataExpt(iexp).av(iav).tag(iexptag).cycTC(1:nCycles);
                    cycTC_AV = cellfun(@(x,y) cat(3,x,y),...
                        antiDataExpt(iexp).av(visualTrials).tag(iexptag).cycTC(1:nCycles),...
                        antiDataExpt(iexp).av(auditoryTrials).tag(iexptag).cycTC(1:nCycles),'unif',0);
                    firstTC_AV = cycTC_AV{1};
                    lateCycTC = [];
                    lateCycTC_AV = [];
                    for icyc = 5:length(cycTC)
                        lateCycTC = cat(3,lateCycTC,...
                            cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
                        lateCycTC_AV = cat(3,lateCycTC_AV,...
                            cycTC_AV{icyc} - mean(cycTC_AV{icyc}(basewin_0,:,:),1));
                    end
                    
                    tarTC_mat = tarDataExpt(iexp).av(iav).tag(iexptag).targetTC;
                    allTargetTC = [];
                    for itar = 1:length(tarTC_mat)
                        allTargetTC = cat(3,allTargetTC,tarTC_mat{itar});
                    end
                        
                    if iav == 1
                        firstRespCells = ttest(...
                            squeeze(mean(firstTC_AV(respwin_opt,:,:),1)),...
                            squeeze(mean(firstTC_AV(basewin_0,:,:),1)),...
                            'dim',2,'tail','right','alpha',cellGroupsAlpha);
                        lateCycRespCells = ttest(...
                            squeeze(mean(lateCycTC_AV(respwin_opt,:,:),1)),...
                            squeeze(mean(lateCycTC_AV(basewin_0,:,:),1)),...
                            'dim',2,'tail','right','alpha',cellGroupsAlpha);
                        lateWinRespCells = ttest(...
                            squeeze(mean(longTC_AV(lateWinFr,:,:),1)),...
                            squeeze(mean(longTC_AV(basewin_0,:,:),1)),...
                            'dim',2,'tail','right','alpha',cellGroupsAlpha);
                        lateWinSuppCells = ttest(...
                            squeeze(mean(longTC_AV(lateWinFr,:,:),1)),...
                            squeeze(mean(longTC_AV(basewin_0,:,:),1)),...
                            'dim',2,'tail','left','alpha',cellGroupsAlpha);
                        allTarRespCells = ttest(squeeze(mean(...
                            allTargetTC(respwin_target,:,:),1)),...
                            squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
                            'dim',2,'tail','right','alpha',cellGroupsAlpha);
                        eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
                            ttest(squeeze(mean(x(respwin_target,:,:),1)),...
                            squeeze(mean(x(basewin_0_target,:,:),1)),...
                            'dim',2,'tail','right','alpha',cellGroupsAlpha),...
                            tarTC_mat,'unif',0)),2) > 0;
                        decodeCellStruct(iexp).tag(itag).cellInd = allTarRespCells | eaTarRespCells | lateCycRespCells;

                        if antiDataExpt(iexp).isBehav
                            behavCells = true(size(firstRespCells));
                            if antiDataExpt(iexp).hasAttn
                                attnCells = true(size(firstRespCells));
                            else
                                attnCells = false(size(firstRespCells));
                            end
                            if antiDataExpt(iexp).passExpt
                                passCells = true(size(firstRespCells));
                            else
                                passCells = false(size(firstRespCells));
                            end
                        else
                            behavCells = false(size(firstRespCells));
                            attnCells = false(size(firstRespCells));
                            passCells = false(size(firstRespCells));
                        end
                        minRespCells = mean(mean(...
                            lateCycTC(respwin_opt,:,:),1),3)' > minRespThreshold;
                        indicatorID = ones(size(firstRespCells)).*find(strcmp(indicatorName,antiDataExpt(iexp).indicator));
                        stimSizeExpt = ones(size(firstRespCells)).*double(antiDataExpt(iexp).stimSize);
                                                
                        aurocTD = nan(length(allTarRespCells),1);
                        aurocTD_test = nan(length(allTarRespCells),1);
                        tarRespAll = cellfun(@(x) ...
                            squeeze(mean(x(respwin_target,:,:),1)) - ...
                            squeeze(mean(x(basewin_0_target,:,:),1)),...
                            cat(2,tarTC_mat,{allTargetTC}),'unif',0);
                        lateCycResp = squeeze(mean(lateCycTC(respwin_opt,:,:),1));
                        for icell = 1:length(allTarRespCells)
                            aurocTD_temp = cellfun(@(x) ....
                                roc_gh(lateCycResp(icell,:),x(icell,:)),...
                                tarRespAll);
                            aurocTD_test_temp = cellfun(@(x) ....
                                ranksum(lateCycResp(icell,:),x(icell,:)),...
                                tarRespAll) < 0.05;
                            aurocTT = roc_gh(tarRespAll{1}(icell,:),...
                                tarRespAll{2}(icell,:));
                            aurocTT_test = ranksum(tarRespAll{1}(icell,:),...
                                tarRespAll{2}(icell,:)) < 0.05;
                            if ~(eaTarRespCells(icell) | allTarRespCells(icell))
                                aurocTD(icell) = aurocTD_temp(3);
                                aurocTD_test(icell) = aurocTD_test_temp(3);
                            elseif ~aurocTT_test
                                aurocTD(icell) = aurocTD_temp(3);
                                aurocTD_test(icell) = aurocTD_test_temp(3);
                            elseif aurocTT_test && aurocTT < 0.5
                                aurocTD(icell) = aurocTD_temp(1);
                                aurocTD_test(icell) = aurocTD_test_temp(1);
                            elseif aurocTT_test && aurocTT >= 0.5
                                aurocTD(icell) = aurocTD_temp(2);
                                aurocTD_test(icell) = aurocTD_test_temp(2);
                            end                            
                        end
                        
                        cellInfo.tag(itag).behav = cat(1,...
                            cellInfo.tag(itag).behav,behavCells);
                        cellInfo.tag(itag).passive = cat(1,...
                            cellInfo.tag(itag).passive,passCells);
                        cellInfo.tag(itag).attention = cat(1,...
                            cellInfo.tag(itag).attention,attnCells); 
                        cellInfo.tag(itag).indicatorID = cat(1,...
                            cellInfo.tag(itag).indicatorID,indicatorID);  
                        cellInfo.tag(itag).stimSize = cat(1,...
                            cellInfo.tag(itag).stimSize,stimSizeExpt);                            
                        cellInfo.tag(itag).firstResp = cat(1,...
                            cellInfo.tag(itag).firstResp,...
                            firstRespCells);
                        cellInfo.tag(itag).lateCycResp = cat(1,...
                            cellInfo.tag(itag).lateCycResp,...
                            lateCycRespCells);
                        cellInfo.tag(itag).lateWinResp = cat(1,...
                            cellInfo.tag(itag).lateWinResp,...
                            lateWinRespCells);
                        cellInfo.tag(itag).lateWinSupp = cat(1,...
                            cellInfo.tag(itag).lateWinSupp,...
                            lateWinSuppCells);
                        cellInfo.tag(itag).targetResp = cat(1,....
                            cellInfo.tag(itag).targetResp,...
                            eaTarRespCells | allTarRespCells);
                        cellInfo.tag(itag).minRespCells = cat(1,...
                            cellInfo.tag(itag).minRespCells,...
                            minRespCells);
                        cellInfo.tag(itag).aurocTD = cat(1,...
                            cellInfo.tag(itag).aurocTD, aurocTD);
                        cellInfo.tag(itag).aurocTD_test = cat(1,...
                            cellInfo.tag(itag).aurocTD_test, aurocTD);
                    else
                        cycTC_Vis = antiDataExpt(iexp).av(visualTrials).tag(iexptag).cycTC(1:nCycles);
                        lateCycTC_Vis = [];
                        for icyc = 5:length(cycTC)
                            lateCycTC_Vis = cat(3,lateCycTC_Vis,...
                                cycTC_Vis{icyc} - mean(cycTC_Vis{icyc}(basewin_0,:,:),1));
                        end
                        lateCycSI = getSelectivityIndex(...
                            squeeze(mean(lateCycTC_Vis(respwin_opt,:,:),1))',...
                            squeeze(mean(lateCycTC(respwin_opt,:,:),1))');
                        cellInfo.tag(itag).lateCycSI = cat(1,...
                            cellInfo.tag(itag).lateCycSI,...
                            lateCycSI');
                    end
                    
                    antiAnalysis.av(iav).tag(itag).longTC = cat(2,antiAnalysis.av(iav).tag(itag).longTC,...
                        mean(longTC,3));
                    antiAnalysis.av(iav).tag(itag).longTCErr = cat(2,antiAnalysis.av(iav).tag(itag).longTCErr,...
                        ste(longTC,3));
                    antiAnalysis.av(iav).tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.av(iav).tag(itag).cycTC,...
                        cycTC,'unif',0);
                    antiAnalysis.av(iav).tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.av(iav).tag(itag).cycTCErr,...
                        cycTC,'unif',0);
                    antiAnalysis.av(iav).tag(itag).lateCycTC = cat(2,antiAnalysis.av(iav).tag(itag).lateCycTC,...
                        mean(lateCycTC,3));
                    antiAnalysis.av(iav).tag(itag).lateCycTCErr = cat(2,antiAnalysis.av(iav).tag(itag).lateCycTCErr,...
                        ste(lateCycTC,3));
                    if antiDataExpt(iexp).passExpt
                        longTC = antiDataExpt(iexp).av(iav).tag(iexptag).longTC_pass;
                        cycTC = antiDataExpt(iexp).av(iav).tag(iexptag).cycTC_pass(1:nCycles);
                        lateCycTC = [];
                        for icyc = 5:length(cycTC)
                            lateCycTC = cat(3,lateCycTC,...
                                cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
                        end                        
                        
                        fprintf('n pass = %s\n', num2str(size(longTC,2)))
                        antiAnalysis.av(iav).tag(itag).longTC_pass = cat(2,antiAnalysis.av(iav).tag(itag).longTC_pass,...
                            mean(longTC,3));
                        antiAnalysis.av(iav).tag(itag).cycTC_pass = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.av(iav).tag(itag).cycTC_pass,...
                            cycTC,'unif',0);
                        antiAnalysis.av(iav).tag(itag).lateCycTC_pass = cat(2,antiAnalysis.av(iav).tag(itag).lateCycTC_pass,...
                            mean(lateCycTC,3));
                        
                        if iav == 1
                            tarTC_mat = tarDataExpt(iexp).av(iav).tag(iexptag).targetTC_pass;
                            eaTarRespCells = sum(cell2mat(cellfun(@(x) ...
                                ttest(squeeze(mean(x(respwin_target,:,:),1)),...
                                squeeze(mean(x(basewin_0_target,:,:),1)),...
                                'dim',2,'tail','right','alpha',cellGroupsAlpha),...
                                tarTC_mat,'unif',0)),2) > 0;
                            allTargetTC = [];
                            for itar = 1:length(tarTC_mat)
                                allTargetTC = cat(3,allTargetTC,tarTC_mat{itar});
                            end
                            allTarRespCells = ttest(squeeze(mean(...
                                allTargetTC(respwin_target,:,:),1)),...
                                squeeze(mean(allTargetTC(basewin_0_target,:,:),1)),...
                                'dim',2,'tail','right','alpha',cellGroupsAlpha);
                            
                            aurocTD = nan(length(allTarRespCells),1);
                            aurocTD_test = nan(length(allTarRespCells),1);
                            tarRespAll = cellfun(@(x) ...
                                squeeze(mean(x(respwin_target,:,:),1)) - ...
                                squeeze(mean(x(basewin_0_target,:,:),1)),...
                                cat(2,tarTC_mat,{allTargetTC}),'unif',0);
                            lateCycResp = squeeze(mean(lateCycTC(respwin_opt,:,:),1));
                            for icell = 1:length(allTarRespCells)
                                aurocTD_temp = cellfun(@(x) ....
                                    roc_gh(lateCycResp(icell,:),x(icell,:)),...
                                    tarRespAll);
                                aurocTD_test_temp = cellfun(@(x) ....
                                    ranksum(lateCycResp(icell,:),x(icell,:)),...
                                    tarRespAll) < 0.05;
                                aurocTT = roc_gh(tarRespAll{1}(icell,:),...
                                    tarRespAll{2}(icell,:));
                                aurocTT_test = ranksum(tarRespAll{1}(icell,:),...
                                    tarRespAll{2}(icell,:)) < 0.05;
                                if ~(eaTarRespCells(icell) | allTarRespCells(icell))
                                    aurocTD(icell) = aurocTD_temp(3);
                                    aurocTD_test(icell) = aurocTD_test_temp(3);
                                elseif ~aurocTT_test
                                    aurocTD(icell) = aurocTD_temp(3);
                                    aurocTD_test(icell) = aurocTD_test_temp(3);
                                elseif aurocTT_test && aurocTT < 0.5
                                    aurocTD(icell) = aurocTD_temp(1);
                                    aurocTD_test(icell) = aurocTD_test_temp(1);
                                elseif aurocTT_test && aurocTT >= 0.5
                                    aurocTD(icell) = aurocTD_temp(2);
                                    aurocTD_test(icell) = aurocTD_test_temp(2);
                                end                            
                            end
                            cellInfo.tag(itag).aurocTD_pass = cat(1,...
                                cellInfo.tag(itag).aurocTD_pass, aurocTD);
                            cellInfo.tag(itag).aurocTD_test_pass = cat(1,...
                                cellInfo.tag(itag).aurocTD_pass, aurocTD);
                        else
                            cycTC_Vis = antiDataExpt(iexp).av(visualTrials).tag(iexptag).cycTC_pass(1:nCycles);
                            lateCycTC_Vis = [];
                            for icyc = 5:length(cycTC)
                                lateCycTC_Vis = cat(3,lateCycTC_Vis,...
                                    cycTC_Vis{icyc} - mean(cycTC_Vis{icyc}(basewin_0,:,:),1));
                            end
                            lateCycSI = getSelectivityIndex(...
                                squeeze(mean(lateCycTC_Vis(respwin_opt,:,:),1))',...
                                squeeze(mean(lateCycTC(respwin_opt,:,:),1))');
                            cellInfo.tag(itag).lateCycSI_pass = cat(1,...
                                cellInfo.tag(itag).lateCycSI_pass,...
                                lateCycSI');                            
                        end
                    else
                        if iav == 1
                            cellInfo.tag(itag).aurocTD_pass = cat(1,...
                                cellInfo.tag(itag).aurocTD_pass, ...
                                nan(size(antiDataExpt(iexp).av(visualTrials).tag(iexptag).longTC,2),1));
                            cellInfo.tag(itag).aurocTD_test_pass = cat(1,...
                                cellInfo.tag(itag).aurocTD_test_pass,...
                                nan(size(antiDataExpt(iexp).av(visualTrials).tag(iexptag).longTC,2),1));
                        else
                            cellInfo.tag(itag).lateCycSI_pass = cat(1,...
                                cellInfo.tag(itag).lateCycSI_pass,...
                                nan(size(antiDataExpt(iexp).av(visualTrials).tag(iexptag).longTC,2),1));
                        end
                    end
                end
            end
        end
    end
end

decodeAnalysis = struct;
iav = 1
for iexp = 1:nexp
    fprintf('Expt %s\n', num2str(iexp))
    if isBehav(iexp)
        decodeAnalysis(iexp).behavior = true;
    else
        decodeAnalysis(iexp).behavior = false;
    end
    dc = decodeDataExpt(iexp).av(iav);
    decodeAnalysis(iexp).av(iav).name = avName{iav};
    for itag = 1:ntag
        decodeAnalysis(iexp).av(iav).tag(itag).name = tagName_unique{itag};
    end
    for iexptag = 1:2
        if strcmp(dc.tag(iexptag).name,'ANY') | strcmp(dc.tag(iexptag).name,'EMX')
            tagInd = strcmp(tagName_unique,'pPYR');
        else
            tagInd = strcmp(tagName_unique,dc.tag(iexptag).name);
        end
        if ~any(tagInd)
            continue
        end
        rng(0) % cells randomized and trials randomized
        cellInd = decodeCellStruct(iexp).tag(tagInd).cellInd;
        if sum(cellInd) > maxCellN
            ind = find(cellInd);
            cellSampleID = randsample(ind,maxCellN);
            cellInd = false(length(cellInd),1);
            cellInd(cellSampleID) = true;
        end
        respAllCells = zscore(dc.tag(iexptag).trResp');
        trOut = dc.tag(iexptag).trOut;
        if iav == 1
            trStimID = discretize(dc.tag(iexptag).trStim,oriBins);
        else
            trStimID = [];
        end
        nStimPerBin = histcounts(trStimID);
        minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
        for istim = 1:nStimBins
            ind = find(trStimID == istim);
            if length(ind) >= minTrN_mdl
                if istim == 1
                    matchTrialsInd = [];
                    if sum(nStimPerBin >= minTrN_mdl) == 2
                        n = minBinN;
                    elseif minBinN == nStimPerBin(istim)
                        error('Expt %s, not enough FA/CR trials',num2str(iexp))
                    else
                        n = (nStimBins-1).*minBinN;
                        if n > length(ind)
                            error('Expt %s, not enough FA/CR trials',num2str(iexp))
                        end
                    end
                    indSample = randsample(ind,n);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                else
                    indSample = randsample(ind,minBinN);
                    matchTrialsInd = cat(2,matchTrialsInd,indSample);
                end
            end
        end
        
        resp = respAllCells(matchTrialsInd,cellInd);
        [detectTrInd, targetTrInd] = getStimAndBehaviorYs(trOut(matchTrialsInd));
        C = eye(size(resp,2));
        p=1;
        [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
        [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
        dv_detect = mean(detectTrInd);
        dv_target = mean(targetTrInd);
        pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);
        pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
        
        decodeAnalysis(iexp).av(iav).tag(tagInd).pctCorrectDetect = pctCorrectDetect_ho;
        decodeAnalysis(iexp).av(iav).tag(tagInd).pctCorrectTarget = pctCorrectTarget_ho;
    end
end
%     tarAnalysis = struct;
%     antiAnalysis_passive = struct;
%     tarAnalysis_passive = struct;
%     cellInfo = struct;
%     for itag = 1:2
%         for iexp = 1:nexp
%             if isempty(antiDataExpt(iexp).tag(itag).name)
%                 continue
%             end
%             antiAnalysis.tag(itag).name = antiDataExpt(iexp).tag(itag).name;
%             firstTC = antiDataExpt(iexp).tag(itag).cycTC{1};
%             longTC = antiDataExpt(iexp).tag(itag).longTC;
%             tarTC = tarDataExpt(iexp).tag(itag).targetTC;
% 
%             firstRespCells = ttest(...
%                 squeeze(mean(firstTC(respwin,:,:),1)),...
%                 squeeze(mean(firstTC(basewin_0,:,:),1)),...
%                 'dim',2,'tail','right','alpha',cellGroupsAlpha);
%             lateWinRespCells = ttest(...
%                 squeeze(mean(longTC(lateWinFr,:,:),1)),...
%                 squeeze(mean(longTC(basewin_0,:,:),1)),...
%                 'dim',2,'tail','right','alpha',cellGroupsAlpha);
%             lateWinSuppCells = ttest(...
%                 squeeze(mean(longTC(lateWinFr,:,:),1)),...
%                 squeeze(mean(longTC(basewin_0,:,:),1)),...
%                 'dim',2,'tail','left','alpha',cellGroupsAlpha);
%             tarRespCells = sum(cell2mat(cellfun(@(x) ttest(...
%                 squeeze(mean(x(respwin,:,:),1)),...
%                 squeeze(mean(x(basewin_0,:,:),1)),...
%                 'dim',2,'tail','right','alpha',cellGroupsAlpha),tarTC,'unif',0)),2)>0;
% 
%             respCellsExpt(iexp).exptName = antiDataExpt(iexp).exptName;
%             respCellsExpt(iexp).tag(itag).firstRespCells = firstRespCells;
%             respCellsExpt(iexp).tag(itag).lateWinRespCells = lateWinRespCells;
%             respCellsExpt(iexp).tag(itag).lateWinSuppCells = lateWinSuppCells;
%             respCellsExpt(iexp).tag(itag).targetRespCells = tarRespCells;
%         end
%         
%             antiAnalysis.tag(itag).longTC = [];
%             antiAnalysis.tag(itag).longTCErr = [];
%             antiAnalysis.tag(itag).cycTC = cell(1,nCycles);
%             antiAnalysis.tag(itag).cycTCErr = cell(1,nCycles);
%             antiAnalysis.tag(itag).lateCycTC = [];
%             antiAnalysis.tag(itag).lateCycTCErr = [];
%             for iexp = 1:nexp
%                 if isempty(antiDataExpt(iexp).tag(itag).name)
%                     continue
%                 end
%                 longTC = antiDataExpt(iexp).tag(itag).longTC;
%                 cycTC = antiDataExpt(iexp).tag(itag).cycTC(1:nCycles);
%                 lateCycTC = [];
%                 for icyc = 5:length(cycTC)
%                     lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
%                 end
% 
%                 antiAnalysis.tag(itag).longTC = cat(2,antiAnalysis.tag(itag).longTC,...
%                     mean(longTC,3));
%                 antiAnalysis.tag(itag).longTCErr = cat(2,antiAnalysis.tag(itag).longTCErr,...
%                     ste(longTC,3));
%                 antiAnalysis.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis.tag(itag).cycTC,...
%                     cycTC,'unif',0);
%                 antiAnalysis.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis.tag(itag).cycTCErr,...
%                     cycTC,'unif',0);
%                 antiAnalysis.tag(itag).lateCycTC = cat(2,antiAnalysis.tag(itag).lateCycTC,...
%                     mean(lateCycTC,3));
%                 antiAnalysis.tag(itag).lateCycTCErr = cat(2,antiAnalysis.tag(itag).lateCycTCErr,...
%                     ste(lateCycTC,3));
% 
%                 lateCycRespCells = ttest(...
%                     squeeze(mean(lateCycTC(respwin,:,:),1)),...
%                     squeeze(mean(lateCycTC(basewin_0,:,:),1)),...
%                     'dim',2,'tail','right','alpha',cellGroupsAlpha);
% 
%                 respCellsExpt(iexp).tag(itag).lateCycRespCells = lateCycRespCells;
%             end
%         
% %         if bxExpt
% %              antiAnalysis_passive.tag(itag).name = {unique(tagNames(ind))};
% %             antiAnalysis_passive.tag(itag).longTC = [];
% %             antiAnalysis_passive.tag(itag).longTCErr = [];
% %             antiAnalysis_passive.tag(itag).cycTC = cell(1,nCycles);
% %             antiAnalysis_passive.tag(itag).cycTCErr = cell(1,nCycles);
% %             antiAnalysis_passive.tag(itag).lateCycTC = [];
% %             antiAnalysis_passive.tag(itag).lateCycTCErr = [];
% %             for iexp = 1:nexp
% %                 if isempty(antiDataExpt_passive(iexp).tag)
% %                     longTC = 
% %             longTC = antiDataExpt_passive(iexp).tag(itag).longTC;
% %             cycTC = antiDataExpt_passive(iexp).tag(itag).cycTC(1:nCycles);
% %             lateCycTC = [];
% %             for icyc = 5:length(cycTC)
% %                 lateCycTC = cat(3,lateCycTC,cycTC{icyc} - mean(cycTC{icyc}(basewin_0,:,:),1));
% %             end
% % 
% %             antiAnalysis_passive.tag(itag).longTC = cat(2,antiAnalysis_passive.tag(itag).longTC,...
% %                 mean(longTC,3));
% %             antiAnalysis_passive.tag(itag).longTCErr = cat(2,antiAnalysis_passive.tag(itag).longTCErr,...
% %                 ste(longTC,3));
% %             antiAnalysis_passive.tag(itag).cycTC = cellfun(@(x,y) cat(2,x,mean(y,3)),antiAnalysis_passive.tag(itag).cycTC,...
% %                 cycTC,'unif',0);
% %             antiAnalysis_passive.tag(itag).cycTCErr = cellfun(@(x,y) cat(2,x,ste(y,3)),antiAnalysis_passive.tag(itag).cycTCErr,...
% %                 cycTC,'unif',0);
% %             antiAnalysis_passive.tag(itag).lateCycTC = cat(2,antiAnalysis_passive.tag(itag).lateCycTC,...
% %                 mean(lateCycTC,3));
% %             antiAnalysis_passive.tag(itag).lateCycTCErr = cat(2,antiAnalysis_passive.tag(itag).lateCycTCErr,...
% %                 ste(lateCycTC,3));
% %             end
% %         end
%         
%         tarAnalysis.tag(itag).tc = cell(1,2);
%         tarAnalysis.tag(itag).rewTC = cell(2,2);
%         tarAnalysis_passive.tag(itag).easyTC = [];
%         for iexp = 1:nexp
%             if isempty(antiDataExpt(iexp).tag(itag).name)
%                 continue
%             end
%             htc = tarDataExpt(iexp).tag(itag).targetTC{1};
%             etc = tarDataExpt(iexp).tag(itag).targetTC{2};
%             tarAnalysis.tag(itag).tc{1} = cat(2,tarAnalysis.tag(itag).tc{1},...
%                 mean(htc,3));
%             tarAnalysis.tag(itag).tc{2} = cat(2,tarAnalysis.tag(itag).tc{2},...
%                 mean(etc,3));
%             htc = tarDataExpt(iexp).tag(itag).rewTC(:,1);
%             etc = tarDataExpt(iexp).tag(itag).rewTC(:,2);
%             tarAnalysis.tag(itag).rewTC{1,1} = cat(2,tarAnalysis.tag(itag).rewTC{1,1},...
%                 mean(htc{1},3));
%             tarAnalysis.tag(itag).rewTC{1,2} = cat(2,tarAnalysis.tag(itag).rewTC{1,2},...
%                 mean(htc{2},3));
%             tarAnalysis.tag(itag).rewTC{2,1} = cat(2,tarAnalysis.tag(itag).rewTC{2,1},...
%                 mean(etc{1},3));
%             tarAnalysis.tag(itag).rewTC{2,2} = cat(2,tarAnalysis.tag(itag).rewTC{2,2},...
%                 mean(etc{2},3));
%             if bxExpt
% %                 tarAnalysis_passive.tag(itag).name = unique(tagNames(ind));
% %                 tc = tarDataExpt_passive(iexp).tag(itag).easyTarTC;
% %                 tarAnalysis_passive.tag(itag).easyTC = cat(2,tarAnalysis_passive.tag(itag).easyTC,...
% %                     mean(tc,3));
%             end
%         end
%         
%         cellInfo.tag(itag).name = antiDataExpt(1).tag(itag).name;
%         cellInfo.tag(itag).isFlex = [];
%         cellInfo.tag(itag).firstRespCells = [];
%         cellInfo.tag(itag).lateWinRespCells = [];
%         cellInfo.tag(itag).lateWinSuppCells = [];
%         cellInfo.tag(itag).targetRespCells = [];
%         cellInfo.tag(itag).lateCycRespCells = [];
%         cellInfo.tag(itag).aurocReward = [];
%         for iexp = 1:nexp
%             if isempty(antiDataExpt(iexp).tag(itag).name)
%                 continue
%             end
%             if itag == 2 && isempty(antiDataExpt(iexp).tag(1).name)
%                 isFlex = true(length(respCellsExpt(iexp).tag(itag).firstRespCells),1);
%             elseif itag == 2
%                 isFlex = false(length(respCellsExpt(iexp).tag(itag).firstRespCells),1);                
%             else
%                 isFlex = [];
%             end
%             
%             htc = tarDataExpt(iexp).tag(itag).rewTC(:,1);    
%             nc = size(htc{1},2);
%             aurocrew = nan(nc,1);
%             for ic = 1:nc
%                 aurocrew(ic) = roc_gh(squeeze(mean(htc{2}(rewwin,ic,:))),...
%                     squeeze(mean(htc{1}(rewwin,ic,:))));
%             end
%             
%             cellInfo.tag(itag).isFlex = cat(1,cellInfo.tag(itag).isFlex,isFlex);
%             cellInfo.tag(itag).firstRespCells = cat(1,cellInfo.tag(itag).firstRespCells,...
%                 logical(cell2mat({respCellsExpt(iexp).tag(itag).firstRespCells}')));
%             cellInfo.tag(itag).lateWinRespCells = cat(1,cellInfo.tag(itag).lateWinRespCells,...
%                 logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinRespCells}')));
%             cellInfo.tag(itag).lateWinSuppCells = cat(1,cellInfo.tag(itag).lateWinSuppCells,...
%                 logical(cell2mat({respCellsExpt(iexp).tag(itag).lateWinSuppCells}')));
%             cellInfo.tag(itag).targetRespCells = cat(1,cellInfo.tag(itag).targetRespCells,...
%                 logical(cell2mat({respCellsExpt(iexp).tag(itag).targetRespCells}')));
%             cellInfo.tag(itag).lateCycRespCells = cat(1,cellInfo.tag(itag).lateCycRespCells,...
%                 logical(cell2mat({respCellsExpt(iexp).tag(itag).lateCycRespCells}')));
%             cellInfo.tag(itag).aurocReward = cat(1,cellInfo.tag(itag).aurocReward,aurocrew);
%         end
%         cellInfo.tag(itag).minRespCells = (mean(antiAnalysis.tag(itag).cycTC{1}(respwin,:),1) > ...
%             minRespThreshold)';
%     end
% %     cellInfo.tag(itag).isTuned = logical(cell2mat({oriTuningExpt.isTuned}))';
% %     cellInfo.tag(itag).oriResp = cell2mat({oriTuningExpt.oriResp}');
% %     cellInfo.tag(itag).oriRespErr = cell2mat({oriTuningExpt.oriRespErr}');
% %     cellInfo.tag(itag).oriFit = cell2mat({oriTuningExpt.fit})';
% %     cellInfo.tag(itag).oriPref = cell2mat({oriTuningExpt.oriPref})';
% %     cellInfo.tag(itag).hwhm = hwhmFromOriFit(cellInfo.oriFit(:,1:180)',1:180)';
% 
% for iexp = 1:nexp
%         cellInd = respCellsExpt(iexp).tag(1).lateCycRespCells | ...
%             respCellsExpt(iexp).tag(1).lateCycRespCells;
%         cellInd = cat(1,cellInd,respCellsExpt(iexp).tag(2).lateCycRespCells | ...
%             respCellsExpt(iexp).tag(2).lateCycRespCells); 
%         decodeDataExpt(iexp).cellInd = cellInd';
% end

% save([fnout 'adaptAnalysis'],'antiAnalysis','tarAnalysis','antiAnalysis_passive','cellInfo','respCellsExpt')

%% plotting params
cycTCLim = [-0.005 0.05];
cycRespLim = [0 0.05];
adaptRespLim = [0 1.4];
hmLim = [-0.3 0.3];
longTCLim = [-0.005, 0.2];
winScatLim = [-0.005 1.3];
respScatLim = [-0.005 0.08];
siLim = [-10 10];

% respTCLim = [-0.005 0.05];
% cycTCLim = [-0.01 0.12];
% cycTCLim_minRespCells = [-0.005 0.025];
% scatLim_win = [-0.2 0.6];
% scatLim_cyc = [-0.035 0.085];
% hmLim = [-0.1 0.1];
% exCellTCLim = [-0.02 0.15];
% oriRespLim = [-0.05 0.15];
% siLim = [-10 10];
% siOriLim = [-3 3];
% oriBarLim_win = [0 0.08];
% oriBarLim_resp = [0 0.04];
% oriLim_taskResp = [-0.005 0.035];
% oriNLim = [0 120];
% oriTCLim = [-0.005 0.08];
% targetTCLim = [-0.015 0.08];
% outTCLim = [-0.005 0.04];
% firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];
% suppTCLim = [-0.05 0.005];
% suppScatLim_win = [-0.2 0.1];
% suppScatLim_cyc = [-0.015 0.015];
% 
tcStartFrame = 26;
cycTCEndTimeMs = 350;
% cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
% ttLabel_target = -1000:250:900;
% preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.av(1).tag(1).longTC,1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_rew = ((ttLabel_target./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)+1);
% ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr_target)+1);
% 
nFr_cyc = size(antiAnalysis.av(1).tag(itag).cycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
% respWinTT_target = (...
%     [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
%     .*(1000/frameRateHz);
% baseWinTT = (...
%     [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
%     .*(1000/frameRateHz);
% 
% movWinLabelFr = 30:(30+nMovWin-1);
% movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% 
% weightLim = [-3 4.2];
% binnedWeightLim = [-0.4 0.4];
% weightLimSum = [-0.8 0.8];
% siLimSum = [-0.5 2.5];
% 
% % lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);
indicatorColors = brewermap(ngcamp+2,'YlGn');
indicatorColors = indicatorColors(3:end,:);
%% heatmaps of responses of each cell type
%attn mice hm
itag = 2
    hmfig = figure;
    colormap(brewermap([],'*RdBu'));
    suptitle({antiAnalysis.av(visualTrials).tag(itag).name;...
        'Anticipation: Any Anticipation-Modulated Cells';...
        'Visual Trials'})
    subplot 221
%     ind = logical(cellInfo.tag(itag).lateWinRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
    ind = true(length(cellInfo.tag(itag).firstResp),1);
    bxInd = cellInfo.tag(itag).behav & cellInfo.tag(itag).attention;
    passInd = cellInfo.tag(itag).passive == 1;
    lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(:,ind & bxInd);
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
    title('Behaving')
    print([fnout '_anti_heatmaps_attn' antiAnalysis.av(visualTrials).tag(itag).name],'-dpdf','-fillpage')
    
    figure
    colormap gray
    subplot 221
    ind2 = ind & bxInd;
    ind_label = ind2(ind2) & (cellInfo.tag(itag).lateWinResp(ind2)|...
        cellInfo.tag(itag).firstResp(ind2)|cellInfo.tag(itag).lateCycResp(ind2));
    labelSort = imcomplement(flipud(ind_label(lateWinSortInd)));
    imagesc(labelSort)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    title('Behaving')
    print([fnout '_anti_heatmaps_attn_label_' antiAnalysis.av(visualTrials).tag(itag).name],'-dpdf','-fillpage')
    
setFigParams4Print('landscape')
for itag = 1:ntag
    hmfig = figure;
    colormap(brewermap([],'*RdBu'));
    celllabelfig = figure;
    colormap gray;
    suptitle({antiAnalysis.av(visualTrials).tag(itag).name;...
        'Anticipation: Any Anticipation-Modulated Cells';...
        'Visual Trials'})
    figure(hmfig)
    subplot 221
%     ind = logical(cellInfo.tag(itag).lateWinRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
    ind = cellInfo.tag(itag).lateWinResp | ...
        cellInfo.tag(itag).lateWinSupp | ...
        cellInfo.tag(itag).firstResp | ...
        cellInfo.tag(itag).lateCycResp;
    bxInd = cellInfo.tag(itag).behav;
    passInd = cellInfo.tag(itag).passive == 1;
    lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(:,ind & bxInd);
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
    title('Behaving')
    figure(celllabelfig)
    subplot 221
    ind2 = ind & bxInd;
    ind_label = ind2(ind2) & (cellInfo.tag(itag).lateWinResp(ind2)|...
        cellInfo.tag(itag).firstResp(ind2)|cellInfo.tag(itag).lateCycResp(ind2))
    labelSort = imcomplement(flipud(ind_label(lateWinSortInd)));
    imagesc(labelSort)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    title('Behaving')
    
    figure(hmfig)
    subplot 223
    lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(:,ind & passInd);
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
    title('Behaving (with Passive Condition)')
    
    subplot 224
    lateWinTC_pass = antiAnalysis.av(visualTrials).tag(itag).longTC_pass(:,ind(passInd));
    hm = flipud(lateWinTC_pass(:,lateWinSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title('Passive')
    
    subplot 222
    lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(:,ind & ~bxInd);
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
    title('Naive')
    print([fnout '_anti_heatmaps_' antiAnalysis.av(visualTrials).tag(itag).name],'-dpdf','-fillpage')
    
    figure(celllabelfig)
    print([fnout '_anti_heatmaps_label_' antiAnalysis.av(visualTrials).tag(itag).name],'-dpdf','-fillpage')
end

%% indicator type analysis
cycTC1_indicator_fig = figure;
suptitle({['Naive Stim 1' ];...
    'Anticipation: Any Anticipation-Modulated Cells';...
    'Visual Trials'})
cycTC2_indicator_fig = figure;
suptitle({['Naive Stim 2' ];...
    'Anticipation: Any Anticipation-Modulated Cells';...
    'Visual Trials'})
for itag = 1:ntag
    ind = cellInfo.tag(itag).lateWinResp | ...
        cellInfo.tag(itag).lateWinSupp | ...
        cellInfo.tag(itag).firstResp | ...
        cellInfo.tag(itag).lateCycResp;
    bxInd = cellInfo.tag(itag).behav;
    figure
    suptitle({['Naive ' antiAnalysis.av(visualTrials).tag(itag).name];...
        'Anticipation: Any Anticipation-Modulated Cells';...
        'Visual Trials'})
    colormap(brewermap([],'*RdBu'));
    [nrows,ncols] = optimizeSubplotDim(ngcamp);
    for igcamp = 1:ngcamp        
        gcampInd = ind & ~bxInd & cellInfo.tag(itag).indicatorID == igcamp;
        subplot(nrows,ncols,igcamp)
        lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(:,gcampInd);
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
        title(indicatorName(igcamp))
    end
    print([fnout '_anti_heatmaps_' antiAnalysis.av(visualTrials).tag(itag).name '_indicator'],'-dpdf','-fillpage')

    figure(cycTC1_indicator_fig)
    [nrows,ncols] = optimizeSubplotDim(ntag);
    subplot(nrows,ncols,itag)
    nc = zeros(1,ngcamp);
    clear L
    for igcamp = 1:ngcamp
        hold on
        gcampInd = ind & ~bxInd & cellInfo.tag(itag).indicatorID == igcamp;
        if sum(gcampInd) > 0
            bl = mean(antiAnalysis.av(visualTrials).tag(itag).cycTC{1}...
                (basewin_0,gcampInd),1);
            y = antiAnalysis.av(visualTrials).tag(itag).cycTC{1}(tcStartFrame:end,gcampInd) - bl;
            yerr = ste(antiAnalysis.av(visualTrials).tag(itag).cycTC{1}...
                (tcStartFrame:end,gcampInd),2);
            h = shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,indicatorColors(igcamp,:));
            L(igcamp) = h.mainLine;
        end
        nc(igcamp) = sum(gcampInd);
    end
    figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    hline(0,'k:')
    vline(respWinTT,'k--');
    title(antiAnalysis.av(visualTrials).tag(itag).name)
    legend(L(nc>0), strcat(indicatorName(nc > 0)',...
        cellfun(@(x) [': n=' num2str(x)],num2cell(nc(nc > 0))','unif',0)),...
        'location','northeast')
    
    figure(cycTC2_indicator_fig)
    [nrows,ncols] = optimizeSubplotDim(ntag);
    subplot(nrows,ncols,itag)
    nc = zeros(1,ngcamp);
    clear L
    for igcamp = 1:ngcamp
        hold on
        gcampInd = ind & ~bxInd & cellInfo.tag(itag).indicatorID == igcamp;
        if sum(gcampInd) > 0
            bl = mean(antiAnalysis.av(visualTrials).tag(itag).cycTC{2}...
                (basewin_0,gcampInd),1);
            y = antiAnalysis.av(visualTrials).tag(itag).cycTC{2}(tcStartFrame:end,gcampInd) - bl;
            yerr = ste(antiAnalysis.av(visualTrials).tag(itag).cycTC{2}...
                (tcStartFrame:end,gcampInd),2);
            h = shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,indicatorColors(igcamp,:));
            L(igcamp) = h.mainLine;
        end
        nc(igcamp) = sum(gcampInd);
    end
    figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    hline(0,'k:')
    vline(respWinTT,'k--');
    title(antiAnalysis.av(visualTrials).tag(itag).name)
    legend(L(nc>0), strcat(indicatorName(nc > 0)',...
        cellfun(@(x) [': n=' num2str(x)],num2cell(nc(nc > 0))','unif',0)),...
        'location','northeast')
end
figure(cycTC1_indicator_fig);
print([fnout '_anti_cyc1TC_indicator'],'-dpdf','-fillpage')
figure(cycTC2_indicator_fig);
print([fnout '_anti_cyc2TC_indicator'],'-dpdf','-fillpage')

%% stim size analysis
    %heatmaps
for itag = 1:ntag
    %heatmaps - naive only
    ind = cellInfo.tag(itag).lateWinResp | ...
        cellInfo.tag(itag).lateWinSupp | ...
        cellInfo.tag(itag).firstResp | ...
        cellInfo.tag(itag).lateCycResp;
    bxInd = cellInfo.tag(itag).behav;
    figure
    suptitle({['Naive ' antiAnalysis.av(visualTrials).tag(itag).name];...
        'Anticipation: Any Anticipation-Modulated Cells';...
        'Visual Trials'})
    colormap(brewermap([],'*RdBu'));
    n = cell(3,nsize+1);
    n(:,1) = {'Size';'Late Suppressed';'Late Resp'};    
    n(1,:) = cat(2,{[]},num2cell(stimSize));
    for isz = 1:nsize
        szInd = ind & ~bxInd & cellInfo.tag(itag).stimSize == stimSize(isz);
        subplot(2,nsize,isz)
        lateWinTC = antiAnalysis.av(visualTrials).tag(itag).longTC(...
            :,szInd);
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
        title(sprintf('Stim size %s deg',num2str(stimSize(isz))))
        n{2,isz+1} = sum(cellInfo.tag(itag).lateWinSupp & szInd)./length(szInd);
        n{3,isz+1} = sum((cellInfo.tag(itag).lateWinResp | ...
            cellInfo.tag(itag).lateCycResp) & szInd)./length(szInd);
    end
    print([fnout '_anti_heatmaps_' antiAnalysis.av(visualTrials).tag(itag).name '_size'],'-dpdf','-fillpage')
    
%     figure
    f=figure;
    suptitle({['Naive ' antiAnalysis.av(visualTrials).tag(itag).name];...
        'Fraction of Late Modulated Cells'})
    uit = uitable(f);
    uit.Data = n;
    print([fnout '_anti_modCellsTable_' antiAnalysis.av(visualTrials).tag(itag).name '_size'],'-dpdf','-fillpage')
    
end


% fraction of suppressed neurons 
fractionFirstStimSupp = nan(nsize,ntag);
fractionTotalSupp = nan(nsize,ntag);
for itag = 1:ntag
    for isize = 1:nsize
        ind_supp = cellInfo.tag(itag).behav == 0 & ...
            cellInfo.tag(itag).stimSize == stimSize(isize) & ...
             cellInfo.tag(itag).lateWinSupp;
        ind_resp = cellInfo.tag(itag).behav == 0 & ...
            cellInfo.tag(itag).stimSize == stimSize(isize) & ...
             cellInfo.tag(itag).firstResp;
        fractionFirstStimSupp(isize,itag) = sum(ind_supp & ind_resp)./sum(ind_resp);
        fractionTotalSupp(isize,itag) = sum(ind_supp)./length(ind_supp);
    end
end

figure
subplot 121
bar(fractionFirstStimSupp'.*100,'grouped')
legend(cellfun(@num2str,num2cell(stimSize),'unif',0))
figXAxis([],'',[0 ntag+1],1:ntag,tagName_unique)
figYAxis([],'Supp/First Resp (%)',[0 25])
figAxForm
subplot 122
bar(fractionTotalSupp'.*100,'grouped')
L = legend(cellfun(@num2str,num2cell(stimSize),'unif',0));
title(L,'Stim Size')
figXAxis([],'',[0 ntag+1],1:ntag,tagName_unique)
figYAxis([],'Supp/Total (%)',[0 25])
figAxForm
print([fnout '_fractionSupp_size'],'-dpdf')
%% adaptation of each cell type, behaving and naive 
figure
suptitle('Visual Trials Only')
nc = nan(ntag,2);
for itag = 1:ntag
    ind = cellInfo.tag(itag).firstResp & cellInfo.tag(itag).minRespCells;
    bxInd = cellInfo.tag(itag).behav;
    nc(itag,1) = sum(ind & bxInd);
    nc(itag,2) = sum(ind & ~bxInd);
    
    cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),...
        antiAnalysis.av(visualTrials).tag(itag).cycTC,'unif',0);
    adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);
    
    subplot 221
    hold on
    y = cellfun(@(x) mean(x(ind & bxInd)),cycResp);
    yerr = cellfun(@(x) ste(x(ind & bxInd),2),cycResp);
    errorbar(1:nCycles,y,yerr,'.-')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    title('Behaving')
    
    subplot 222
    hold on
    y = cellfun(@(x) mean(x(ind & ~bxInd)),cycResp);
    yerr = cellfun(@(x) ste(x(ind & ~bxInd),2),cycResp);
    errorbar(1:nCycles,y,yerr,'.-')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',cycTCLim)
    figAxForm
    title('Naive')
    
    subplot 223
    hold on
    y = cellfun(@(x) mean(x(ind & bxInd)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind & bxInd),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.-')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',adaptRespLim)
    figAxForm
    title('Behaving')
    
    subplot 224
    hold on
    y = cellfun(@(x) mean(x(ind & ~bxInd)),adaptResp);
    yerr = cellfun(@(x) ste(x(ind & ~bxInd),2),adaptResp);
    errorbar(1:nCycles,y,yerr,'.-')
    figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
    figYAxis([],'dF/F',adaptRespLim)
    figAxForm
    title('Naive')
end
subplot 221
legend(strcat(tagName_unique,...
    cellfun(@(x) [': n=' num2str(x)],num2cell(nc(:,1))','unif',0)),...
    'location','northeast')
subplot 222
legend(strcat(tagName_unique,...
    cellfun(@(x) [': n=' num2str(x)],num2cell(nc(:,2))','unif',0)),...
    'location','northeast')
print([fnout 'adaptCellTypes'],'-dpdf','-fillpage')

%% attention modulation of each cell type, behaving mice with attention and naive mice

setFigParams4Print('portrait')
for itag = 1:ntag
    siFig = figure;
    n = zeros(1,2);
    clear si_mat;
    hp_mat = [];
    si_mat = [];
    suptitle(sprintf('%s Cells',antiAnalysis.av(iav).tag(itag).name))
    for ibx = 1:2        
        if ibx == 1
            bxInd = cellInfo.tag(itag).behav & cellInfo.tag(itag).attention;
            miceInd = cellfun(@(x) any(strcmp(x,antiAnalysis.av(1).tag(itag).name)),...
                {mice.tagTypes}) & cell2mat({mice.attention}) == 1;
            if sum(miceInd) == 0
                continue
            end
            figure
            suptitle(sprintf('Behavior Mice with Attention, n=%s',...
                num2str(sum(miceInd))));
        else
            bxInd = ~cellInfo.tag(itag).behav;
            miceInd = cellfun(@(x) any(strcmp(x,antiAnalysis.av(1).tag(itag).name)),...
                {mice.tagTypes}) & cell2mat({mice.behav}) == 0;
            figure
            suptitle(sprintf('Naive Mice, n=%s',...
                num2str(sum(miceInd))));
        end
        subplot 311    
        ind = cellInfo.tag(itag).lateWinResp | ...
            cellInfo.tag(itag).firstResp | ...
            cellInfo.tag(itag).lateCycResp;
        for iav = 1:2
            y = mean(antiAnalysis.av(iav).tag(itag).longTC...
                (tcStartFrame:end,ind & bxInd),2);
            yerr = ste(antiAnalysis.av(iav).tag(itag).longTC...
                (tcStartFrame:end,ind & bxInd),2);
            hold on
            shadedErrorBar_chooseColor(tt_longTC,y,yerr,cueColor{iav});
        end
        figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
        figYAxis([],'dF/F',longTCLim)  
        vline(lateWinTT,'k--')
        hline(0,'k:')
        figAxForm([],0)
        title(sprintf('Anti Responsive %s Cells (%s/%s)',...
            antiAnalysis.av(iav).tag(itag).name,...
            num2str(sum(ind & bxInd)),...
            num2str(sum(bxInd))))

        subplot 323
        x = mean(antiAnalysis.av(visualTrials).tag(itag).longTC(...
            lateWinFr,ind & bxInd),1);
        xerr = ste(x,2);    
        y = mean(antiAnalysis.av(auditoryTrials).tag(itag).longTC(...
            lateWinFr,ind & bxInd),1);
        yerr = ste(y,2);
        plot(x,y,'.')
        hold on
        errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.','MarkerSize',20)
        plot(winScatLim,winScatLim,'k--')
        [~,p] = ttest(x,y);
        figXAxis([],'Visual (dF/F)',winScatLim)
        figYAxis([],'Auditory (dF/F)',winScatLim)
        figAxForm
        title(sprintf('Late Window, All Resp. Cells (%s/%s), p = %s',...
            num2str(sum(ind & bxInd)),...
            num2str(sum(bxInd)),...
            num2str(round(p,2,'significant'))))

        ind = cellInfo.tag(itag).lateCycResp & bxInd;
        subplot 325
        for iav = 1:2
            y = mean(antiAnalysis.av(iav).tag(itag).lateCycTC...
                (tcStartFrame:end,ind),2);
            yerr = ste(antiAnalysis.av(iav).tag(itag).lateCycTC...
                (tcStartFrame:end,ind),2);
            hold on
            shadedErrorBar_chooseColor(tt_cycTC,y,yerr,cueColor{iav});
        end
        figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
        figYAxis([],'dF/F',cycTCLim)  
        hline(0,'k:')
        vline(respWinTT,'k--')
        figAxForm
        title(sprintf('Late Cycle Resp. Cells (%s/%s)',...
            num2str(sum(ind)),...
            num2str(sum(bxInd))))

        subplot 326
        x = mean(antiAnalysis.av(visualTrials).tag(itag).lateCycTC...
                (respwin,ind),1);
        xerr = ste(x,2);
        y = mean(antiAnalysis.av(auditoryTrials).tag(itag).lateCycTC...
                (respwin,ind),1);
        yerr = ste(y,2);
        plot(x,y,'.')
        hold on
        errorbar(mean(x),mean(y),yerr,yerr,xerr,xerr,'.','MarkerSize',20)
        plot(respScatLim,respScatLim,'k--')
        figXAxis([],'Visual (dF/F)',respScatLim)
        figYAxis([],'Auditory (dF/F)',respScatLim)
        figAxForm
        [~,p] = ttest(x,y);
        title(sprintf('Late Cycle Resp. Cells (%s/%s),p=%s',...
            num2str(sum(ind)),...
            num2str(sum(bxInd)),...
            num2str(round(p,2,'significant'))))
        if ibx == 1
            print([fnout 'attn_behav_' antiAnalysis.av(1).tag(itag).name],'-dpdf','-fillpage')
        else 
            print([fnout 'attn_naive_' antiAnalysis.av(1).tag(itag).name],'-dpdf','-fillpage')
        end
        
        figure(siFig)
        si_resp = cellInfo.tag(itag).lateCycSI(...
            cellInfo.tag(itag).lateCycResp & bxInd);
        hp = cdfplot(si_resp);
        n(ibx) = length(si_resp);
        hold on 
        figXAxis([],'Selectivity Index',siLim)
        figYAxis([],'Fraction of Cells',[0 1])
        figAxForm
        vline(0,'k-')
        h = vline(mean(si_resp),'-');
        h.Color = hp.Color;
        hp_mat = [hp_mat,hp]; 
        si_mat = cat(2,si_mat,{si_resp});
        if length(hp_mat) == 2
            legend(hp_mat,{['Behav., n=' num2str(n(1))],...
                ['Naive, n=',num2str(n(2))]},'location','northwestoutside')
            [~,p] = kstest2(si_mat{1},si_mat{2});
            title(sprintf('p = %s', num2str(round(p,2,'significant'))))
        else
            if ibx == 1
                title('Behavior')
            else
                title('Naive')
            end
        end
        print([fnout 'attn_SI_' antiAnalysis.av(1).tag(itag).name],'-dpdf')
    end
end

%% passive vs. naive analysis

figure
for itag = 1:ntag
    ind = cellInfo.tag(itag).passive==1;
    subplot(1,ntag,itag)
    plot(cellInfo.tag(itag).lateCycSI(ind),cellInfo.tag(itag).lateCycSI_pass(ind),'.','MarkerSize',10)
    figXAxis([],'Behaving V-A Selectivity',siLim)
    figYAxis([],'Passive V-A Selectivity',siLim)
    figAxForm
    hold on
    hline(0,'k:')
    vline(0,'k:')
end

figure
for itag = 1:ntag
    ind = cellInfo.tag(itag).passive & (cellInfo.tag(itag).lateCycResp | cellInfo.tag(itag).targetResp);
    subplot(1,ntag,itag)
    plot(cellInfo.tag(itag).aurocTD(ind),cellInfo.tag(itag).aurocTD_pass(ind),'.','MarkerSize',10)
    figXAxis([],'Behaving auROC (T vs. D)',[0 1])
    figYAxis([],'Passive auROC (T vs. D)',[0 1])  
    figAxForm
    hold on
    hline(.5,'k:')
    vline(0.5,'k:')  
    plot([0 1],[0 1],'k--')
end

%% decoding analysis

figure
suptitle('open circles = naive')
iav = 1;
subplot 321
hold on
pctCorrect_stim = nan(nexp,ntag);
for itag = 1:ntag
    for iexp = 1:nexp
        if ~isfield(decodeAnalysis(iexp).av(iav).tag(itag),'pctCorrectTarget')
            continue
        elseif isempty(decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectTarget)
            continue
        else
            y = decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectTarget;
            pctCorrect_stim(iexp,itag) = y;
            h = plot(itag,y,'o','MarkerSize',5);
            h.MarkerEdgeColor = 'k';
            h.LineWidth = 1;
            if decodeAnalysis(iexp).behavior
                h.MarkerFaceColor = [0 0 0];
            else
                h.MarkerFaceColor = [1 1 1];
            end
        end
    end
end
y = nanmean(pctCorrect_stim(isBehav,:),1);
yerr = ste(pctCorrect_stim(isBehav,:),1);
h = errorbar(1:ntag,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = h.MarkerEdgeColor;
y = nanmean(pctCorrect_stim(~isBehav,:),1);
yerr = ste(pctCorrect_stim(~isBehav,:),1);
h = errorbar(1:ntag,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Cell Type',[0 ntag+1],1:ntag,tagName_unique)
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5, 'k:')
cellTypeGroup = repmat(1:ntag,[nexp,1]);
behavGroup = repmat(isBehav',[1,ntag]);
[p,~,stats] = anovan(pctCorrect_stim(:),{cellTypeGroup(:),behavGroup(:)},...
    'model','interaction','display','off');
title(sprintf('Stimulus Model, p-cell:%s,p-bx:%s',num2str(p(1)),...
    num2str(p(2))))

subplot 322
hold on
ind = ~isnan(pctCorrect_stim(:,1)) & ~isnan(pctCorrect_stim(:,3)) & isBehav';
h = plot(pctCorrect_stim(ind,1),pctCorrect_stim(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [1 0.5 0];
h.MarkerEdgeColor = [1 0.5 0];
ind = ~isnan(pctCorrect_stim(:,1)) & ~isnan(pctCorrect_stim(:,3)) & ~isBehav';
h = plot(pctCorrect_stim(ind,1),pctCorrect_stim(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [1 1 1];
h.MarkerEdgeColor = [1 0.5 0];
ind = ~isnan(pctCorrect_stim(:,2)) & ~isnan(pctCorrect_stim(:,3)) & isBehav';
h = plot(pctCorrect_stim(ind,2),pctCorrect_stim(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [0 0.5 1];
h.MarkerEdgeColor = [0 0.5 1];
ind = ~isnan(pctCorrect_stim(:,2)) & ~isnan(pctCorrect_stim(:,3)) & ~isBehav';
h = plot(pctCorrect_stim(ind,2),pctCorrect_stim(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [1 1 1];
h.MarkerEdgeColor = [0 0.5 1];
figXAxis([],'Interneuron Frac. Correct',[0.4 1])
figYAxis([],'Pyr. Neuron Frac. Correct',[0.4 1])
figAxForm
plot([0.4 1],[0.4 1],'k-')
hline(0.5,'k:')
vline(0.5,'k:')
ind = (~isnan(pctCorrect_stim(:,1))|~isnan(pctCorrect_stim(:,2))) &...
    ~isnan(pctCorrect_stim(:,3));
[~,p] = ttest(nansum(pctCorrect_stim(ind,1:2),2),pctCorrect_stim(ind,3));
title(sprintf('org=%s, bl=%s, p=%s',tagName_unique{1},tagName_unique{2},num2str(p)))

subplot 323
hold on
pctCorrect_ch = nan(nexp,ntag);
for itag = 1:ntag
    for iexp = 1:nexp
        if ~isfield(decodeAnalysis(iexp).av(iav).tag(itag),'pctCorrectDetect')
            continue
        elseif ~decodeAnalysis(iexp).behavior
            continue
        elseif isempty(decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectDetect)
            continue
        else
            y = decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectDetect;
            pctCorrect_ch(iexp,itag) = y;
            h = plot(itag,y,'o','MarkerSize',5);
            h.MarkerEdgeColor = 'k';
            h.LineWidth = 1;
            h.MarkerFaceColor = [0 0 0];
        end
    end
end
y = nanmean(pctCorrect_ch,1);
yerr = ste(pctCorrect_ch,1);
h = errorbar(1:ntag,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = h.MarkerEdgeColor;
figXAxis([],'Cell Type',[0 ntag+1],1:ntag,tagName_unique)
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5, 'k:')
[p,~,stats] = anova1(pctCorrect_ch,[],'off');
title(sprintf('Choice Model, p-cell:%s',num2str(p(1))))

subplot 324
hold on
ind = ~isnan(pctCorrect_ch(:,1)) & ~isnan(pctCorrect_ch(:,3)) & isBehav';
h = plot(pctCorrect_ch(ind,1),pctCorrect_ch(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [1 0.5 0];
h.MarkerEdgeColor = [1 0.5 0];
ind = ~isnan(pctCorrect_ch(:,2)) & ~isnan(pctCorrect_ch(:,3)) & isBehav';
h = plot(pctCorrect_ch(ind,2),pctCorrect_ch(ind,3),'o','LineWidth',1);
h.MarkerFaceColor = [0 0.5 1];
h.MarkerEdgeColor = [0 0.5 1];
figXAxis([],'Interneuron Frac. Correct',[0.4 1])
figYAxis([],'Pyr. Neuron Frac. Correct',[0.4 1])
figAxForm
plot([0.4 1],[0.4 1],'k-')
hline(0.5,'k:')
vline(0.5,'k:')
ind = (~isnan(pctCorrect_ch(:,1))|~isnan(pctCorrect_ch(:,2))) &...
    ~isnan(pctCorrect_ch(:,3));
[~,p] = ttest(nansum(pctCorrect_ch(ind,1:2),2),pctCorrect_ch(ind,3));
title(sprintf('org=%s, bl=%s, p=%s',tagName_unique{1},tagName_unique{2},num2str(p)))

tagColors = [1 0.5 0;...    %PV
            0 0.5 1;...     %SOM
            0 0 0];         %PYR
subplot 325
hold on
for itag = 1:ntag-1
    ind = ~isnan(pctCorrect_ch(:,itag)) & ~isnan(pctCorrect_stim(:,itag)) & isBehav';
    h = plot(1:2,[pctCorrect_stim(ind,itag),pctCorrect_ch(ind,itag)],'o-','LineWidth',1);
    for i = 1:sum(ind)
        h(i).MarkerFaceColor = tagColors(itag,:);
        h(i).MarkerEdgeColor = tagColors(itag,:);
        h(i).Color = tagColors(itag,:);
    end
end
figXAxis([],'Model',[0 3],1:2,{'Stim.','Ch.'})
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5,'k:')

print([fnout 'decodingModel_compareCellType'],'-dpdf','-fillpage')

%% just SOM cells 

figure
suptitle('open circles = naive')
iav = 1;
subplot 321
hold on
pctCorrect_stim = nan(nexp,2);
for itag = 2:ntag
    for iexp = 1:nexp
        if ~isfield(decodeAnalysis(iexp).av(iav).tag(itag),'pctCorrectTarget')
            continue
        elseif isempty(decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectTarget)
            continue
        else
            y = decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectTarget;
            pctCorrect_stim(iexp,itag-1) = y;
            h = plot(itag-1,y,'o','MarkerSize',5);
            h.MarkerEdgeColor = 'k';
            h.LineWidth = 1;
            if decodeAnalysis(iexp).behavior
                h.MarkerFaceColor = [0 0 0];
            else
                h.MarkerFaceColor = [1 1 1];
            end
        end
    end
end
y = nanmean(pctCorrect_stim(isBehav,:),1);
yerr = ste(pctCorrect_stim(isBehav,:),1);
h = errorbar(1:2,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = h.MarkerEdgeColor;
y = nanmean(pctCorrect_stim(~isBehav,:),1);
yerr = ste(pctCorrect_stim(~isBehav,:),1);
h = errorbar(1:2,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = [1 1 1];
figXAxis([],'Cell Type',[0 3],1:2,tagName_unique(2:3))
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5, 'k:')
cellTypeGroup = repmat(1:2,[nexp,1]);
behavGroup = repmat(isBehav',[1,2]);
[p,~,stats] = anovan(pctCorrect_stim(:),{cellTypeGroup(:),behavGroup(:)},...
    'model','interaction','display','off');
title(sprintf('Stimulus Model, p-cell:%s,p-bx:%s',num2str(p(1)),...
    num2str(p(2))))

subplot 322
hold on
ind = ~isnan(pctCorrect_stim(:,1)) & ~isnan(pctCorrect_stim(:,2)) & isBehav';
h = plot(pctCorrect_stim(ind,1),pctCorrect_stim(ind,2),'o','LineWidth',1);
h.MarkerFaceColor = [0 0.5 1];
h.MarkerEdgeColor = [0 0.5 1];
ind = ~isnan(pctCorrect_stim(:,1)) & ~isnan(pctCorrect_stim(:,2)) & ~isBehav';
h = plot(pctCorrect_stim(ind,1),pctCorrect_stim(ind,2),'o','LineWidth',1);
h.MarkerFaceColor = [1 1 1];
h.MarkerEdgeColor = [0 0.5 1];
figXAxis([],'SOM Neuron Frac. Correct',[0.4 1])
figYAxis([],'Pyr. Neuron Frac. Correct',[0.4 1])
figAxForm
plot([0.4 1],[0.4 1],'k-')
hline(0.5,'k:')
vline(0.5,'k:')
ind = ~isnan(pctCorrect_stim(:,1)) &...
    ~isnan(pctCorrect_stim(:,2));
[~,p] = ttest(pctCorrect_stim(ind,1),pctCorrect_stim(ind,2));
title(sprintf('org=%s, bl=%s, p=%s',tagName_unique{1},tagName_unique{2},...
    num2str(round(p,2,'significant'))))

subplot 323
hold on
pctCorrect_ch = nan(nexp,2);
for itag = 2:ntag
    for iexp = 1:nexp
        if ~isfield(decodeAnalysis(iexp).av(iav).tag(itag),'pctCorrectDetect')
            continue
        elseif ~decodeAnalysis(iexp).behavior
            continue
        elseif isempty(decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectDetect)
            continue
        else
            y = decodeAnalysis(iexp).av(iav).tag(itag).pctCorrectDetect;
            pctCorrect_ch(iexp,itag-1) = y;
            h = plot(itag-1,y,'o','MarkerSize',5);
            h.MarkerEdgeColor = 'k';
            h.LineWidth = 1;
            h.MarkerFaceColor = [0 0 0];
        end
    end
end
y = nanmean(pctCorrect_ch,1);
yerr = ste(pctCorrect_ch,1);
h = errorbar(1:2,y,yerr,'o','MarkerSize',7);
h.LineWidth = 1;
h.MarkerFaceColor = h.MarkerEdgeColor;
figXAxis([],'Cell Type',[0 3],1:2,tagName_unique(2:3))
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5, 'k:')
[p,~,stats] = anova1(pctCorrect_ch,[],'off');
title(sprintf('Choice Model, p-cell:%s',num2str(p(1))))

subplot 324
hold on
ind = ~isnan(pctCorrect_ch(:,1)) & ~isnan(pctCorrect_ch(:,2)) & isBehav';
h = plot(pctCorrect_ch(ind,1),pctCorrect_ch(ind,2),'o','LineWidth',1);
h.MarkerFaceColor = [0 0.5 1];
h.MarkerEdgeColor = [0 0.5 1];
figXAxis([],'Interneuron Frac. Correct',[0.4 1])
figYAxis([],'Pyr. Neuron Frac. Correct',[0.4 1])
figAxForm
plot([0.4 1],[0.4 1],'k-')
hline(0.5,'k:')
vline(0.5,'k:')
[~,p] = ttest(pctCorrect_ch(ind,1),pctCorrect_ch(ind,2));
title(sprintf('org=%s, bl=%s, p=%s',tagName_unique{1},tagName_unique{2},num2str(p)))

tagColors = [1 0.5 0;...    %PV
            0 0.5 1;...     %SOM
            0 0 0];         %PYR
subplot 325
hold on
for itag = 2
    ind = ~isnan(pctCorrect_ch(:,itag)) & ~isnan(pctCorrect_stim(:,itag)) & isBehav';
    h = plot(1:2,[pctCorrect_stim(ind,itag),pctCorrect_ch(ind,itag)],'-','LineWidth',1);
    for i = 1:sum(ind)
        h(i).MarkerFaceColor = tagColors(itag,:);
        h(i).MarkerEdgeColor = tagColors(itag,:);
        h(i).Color = tagColors(itag,:);
    end
    
end
figXAxis([],'Model',[0 3],1:2,{'Stim.','Ch.'})
figYAxis([],'Frac. Correct',[0 1])
figAxForm
hline(0.5,'k:')

print([fnout 'decodingModel_compareCellType_SOM'],'-dpdf','-fillpage')

% %% plotting params
% respTCLim = [-0.005 0.05];
% cycTCLim = [-0.01 0.12];
% cycTCLim_minRespCells = [-0.005 0.025];
% scatLim_win = [-0.2 0.6];
% scatLim_cyc = [-0.035 0.085];
% hmLim = [-0.1 0.1];
% exCellTCLim = [-0.02 0.15];
% oriRespLim = [-0.05 0.15];
% siLim = [-10 10];
% siOriLim = [-3 3];
% oriBarLim_win = [0 0.08];
% oriBarLim_resp = [0 0.04];
% oriLim_taskResp = [-0.005 0.035];
% oriNLim = [0 120];
% oriTCLim = [-0.005 0.08];
% targetTCLim = [-0.015 0.08];
% outTCLim = [-0.005 0.04];
% firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];
% suppTCLim = [-0.05 0.005];
% suppScatLim_win = [-0.2 0.1];
% suppScatLim_cyc = [-0.015 0.015];
% 
% tcStartFrame = 26;
% cycTCEndTimeMs = 350;
% cycTCEndFr = 45;
% ttLabel_long = 0:500:2500;
% ttLabel_cyc = -200:100:cycTCEndTimeMs;
% ttLabel_target = -1000:250:900;
% preTargetStimLabel = -700:350:0;
% nFr_long = size(antiAnalysis.tag(1).longTC,1);
% tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
% ttLabelFr_rew = ((ttLabel_target./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr)+1);
% ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
%     ((nBaselineFr+nVisDelayFr_target)+1);
% 
% nFr_cyc = size(antiAnalysis.tag(itag).cycTC{1,1},1);
% tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);
% 
% lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
%     .*(1000/frameRateHz);
% respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
%     .*(1000/frameRateHz);
% respWinTT_target = (...
%     [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
%     .*(1000/frameRateHz);
% baseWinTT = (...
%     [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
%     .*(1000/frameRateHz);
% 
% movWinLabelFr = 30:(30+nMovWin-1);
% movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% 
% weightLim = [-3 4.2];
% binnedWeightLim = [-0.4 0.4];
% weightLimSum = [-0.8 0.8];
% siLimSum = [-0.5 2.5];
% 
% % lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);
% 
% %% auROC T vs D
% lateStimAuROC = cell(1,2);
% lateStimAuROCTest = cell(1,2);
% lateStimAuROC_passive = cell(1,2);
% lateStimAuROCTest_passive = cell(1,2);
% for iexp = 1:nexp
%     for itag = 1:2
%         if isempty(tarDataExpt(iexp).tag(itag).targetTC)
%             continue
%         end
%         tarVisResp = cellfun(@(x) ...
%             squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
%             cat(2,tarDataExpt(iexp).tag(itag).targetTC,...
%             {cat(3,tarDataExpt(iexp).tag(itag).targetTC{1},...
%             tarDataExpt(iexp).tag(itag).targetTC{2})}),'unif',0);
%         lateCycResp = [];
%         for icyc = 5:length(cycTC)
%             lateCycResp = cat(2,lateCycResp,...
%                 squeeze(mean(antiDataExpt(iexp).tag(itag).cycTC{icyc}(respwin,:,:),1)...
%                 - mean(antiDataExpt(iexp).tag(itag).cycTC{icyc}(basewin_0,:,:),1)));
%         end
%         
%         nc = size(lateCycResp,1);
%         auroc_target = nan(nc,1);
%         auroc_target_test = nan(nc,1);
%         auroc = nan(nc,3);
%         auroc_test = nan(nc,3);
%         for icell = 1:nc
%             auroc_target(icell) = roc_gh(...
%                 tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
%             auroc_target_test(icell) = ranksum(...
%                 tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) <0.05;
%             
%             auroc(icell,:) = cellfun(@(x) ...
%                 roc_gh(lateCycResp(icell,:),x(icell,:)),tarVisResp);
%             auroc_test(icell,:) = cellfun(@(x) ...
%                 ranksum(lateCycResp(icell,:),x(icell,:)),tarVisResp) < 0.05;
%         end
%         
%         auroc_specific = nan(nc,1);
%         auroc_specific_test = nan(nc,1);
%         for icell = 1:nc
%             if cellInfo.tag(itag).targetRespCells(icell) == 0
%                 auroc_specific(icell) = auroc(icell,3);
%                 auroc_specific_test(icell) = auroc_test(icell,3);
%             elseif auroc_target_test(icell) == 0
%                 auroc_specific(icell) = auroc(icell,3);
%                 auroc_specific_test(icell) = auroc_test(icell,3);
%             elseif auroc_target_test(icell) == 1 && auroc_target(icell) < 0.5
%                 auroc_specific(icell) = auroc(icell,1);
%                 auroc_specific_test(icell) = auroc_test(icell,1);
%             elseif auroc_target_test(icell) == 1 && auroc_target(icell) > 0.5
%                 auroc_specific(icell) = auroc(icell,2);
%                 auroc_specific_test(icell) = auroc_test(icell,2);
%             end
%         end
%         
%         
%         lateStimAuROC{itag} = cat(1,lateStimAuROC{itag},auroc_specific);
%         lateStimAuROCTest{itag} = cat(1,lateStimAuROCTest{itag},auroc_specific_test);
%         
%         if bxExpt
%             if isempty(antiDataExpt_passive(iexp).tag)
%                 lateStimAuROC_passive{itag} = cat(1,...
%                     lateStimAuROC_passive{itag},nan(nc,1));
%                 lateStimAuROCTest_passive{itag} = cat(1,...
%                     lateStimAuROCTest_passive{itag},nan(nc,1));
%             else
%                 tarVisResp = cellfun(@(x) ...
%                     squeeze(mean(x(respwin_target,:,:),1) - mean(x(basewin_0_target,:,:),1)),...
%                     cat(2,tarDataExpt_passive(iexp).tag(itag).targetTC,...
%                     {cat(3,tarDataExpt_passive(iexp).tag(itag).targetTC{1},...
%                     tarDataExpt_passive(iexp).tag(itag).targetTC{2})}),'unif',0);
%                 lateCycResp = [];
%                 for icyc = 5:length(cycTC)
%                     lateCycResp = cat(2,lateCycResp,...
%                         squeeze(mean(antiDataExpt_passive(iexp).tag(itag).cycTC{icyc}(respwin,:,:),1)...
%                         - mean(antiDataExpt_passive(iexp).tag(itag).cycTC{icyc}(basewin_0,:,:),1)));
%                 end
%                 nc = size(lateCycResp,1);
%                 auroc_target = nan(nc,1);
%                 auroc_target_test = nan(nc,1);
%                 auroc = nan(nc,3);
%                 auroc_test = nan(nc,3);
%                 for icell = 1:nc
%                     auroc_target(icell) = roc_gh(...
%                         tarVisResp{1}(icell,:),tarVisResp{2}(icell,:));
%                     auroc_target_test(icell) = ranksum(...
%                         tarVisResp{1}(icell,:),tarVisResp{2}(icell,:)) <0.05;
% 
%                     auroc(icell,:) = cellfun(@(x) ...
%                         roc_gh(lateCycResp(icell,:),x(icell,:)),tarVisResp);
%                     auroc_test(icell,:) = cellfun(@(x) ...
%                         ranksum(lateCycResp(icell,:),x(icell,:)),tarVisResp) < 0.05;
%                 end
% 
%                 auroc_specific = nan(nc,1);
%                 auroc_specific_test = nan(nc,1);
%                 for icell = 1:nc
%                     if cellInfo.tag(itag).targetRespCells(icell) == 0
%                         auroc_specific(icell) = auroc(icell,3);
%                         auroc_specific_test(icell) = auroc_test(icell,3);
%                     elseif auroc_target_test(icell) == 0
%                         auroc_specific(icell) = auroc(icell,3);
%                         auroc_specific_test(icell) = auroc_test(icell,3);
%                     elseif auroc_target_test(icell) == 1 && auroc_target(icell) < 0.5
%                         auroc_specific(icell) = auroc(icell,1);
%                         auroc_specific_test(icell) = auroc_test(icell,1);
%                     elseif auroc_target_test(icell) == 1 && auroc_target(icell) > 0.5
%                         auroc_specific(icell) = auroc(icell,2);
%                         auroc_specific_test(icell) = auroc_test(icell,2);
%                     end
%                 end
%                 lateStimAuROC_passive{itag} = cat(1,...
%                     lateStimAuROC_passive{itag},auroc_specific);
%                 lateStimAuROCTest_passive{itag} = cat(1,...
%                     lateStimAuROCTest_passive{itag},auroc_specific_test);
%             end
%         end
%     end
% end
% 
% figure
% suptitle('Target or Distractor Responsive')
% for itag = 1:2
%     ind = ~isnan(lateStimAuROC_passive{itag}) & ...
%         (cellInfo.tag(itag).lateCycRespCells | cellInfo.tag(itag).targetRespCells);
%     subplot(3,2,itag)    
%     h1 = cdfplot(lateStimAuROC{itag}(ind));
%     hold on
%     h2 = cdfplot(lateStimAuROC_passive{itag}(ind));
%     figXAxis([],'auROC (Target vs. Distractor)',[0 1])
%     figYAxis([],'Fraction of Cells',[0 1])
%     figAxForm
%     h3 = vline(mean(lateStimAuROC{itag}(ind)),'-');
%     h3.Color = h1.Color;
%     h3 = vline(mean(lateStimAuROC_passive{itag}(ind)),'-');
%     h3.Color = h2.Color;
%     vline(0.5,'k:')
%         
%     n = sum(ind);
%     [~,p] = kstest2(lateStimAuROC{itag}(ind),lateStimAuROC_passive{itag}(ind));
%     title(sprintf('%s Cells, n=%s, p=%s',tagName_unique{itag},num2str(n),...
%         num2str(round(p,2,'significant'))))
%     legend([h1,h2],{'Behav.','Pass'},'location','northeastoutside')
%     
%     subplot(3,2,2+itag)
%     plot(lateStimAuROC{itag}(ind),lateStimAuROC_passive{itag}(ind),...
%         '.','MarkerSize',10)
%     hold on
%     plot([0 1],[0 1],'k--')
%     figXAxis([],'Behaving',[0 1])
%     figYAxis([],'Passive Viewing',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     vline(0.5,'k:') 
%     title(sprintf('%s Cells',tagName_unique{itag}))   
%     
%     nTaskTune = [sum(ind & ~lateStimAuROCTest{itag});...
%         sum(ind & lateStimAuROCTest{itag} & lateStimAuROC{itag} < 0.5);...
%         sum(ind & lateStimAuROCTest{itag} & lateStimAuROC{itag} > 0.5)];
%     nTaskTune_pass = [sum(ind & lateStimAuROCTest_passive{itag}==0);...
%         sum(ind & lateStimAuROCTest_passive{itag}==1 & lateStimAuROC_passive{itag} < 0.5);...
%         sum(ind & lateStimAuROCTest_passive{itag}==1 & lateStimAuROC_passive{itag} > 0.5)];
%     subplot(3,2,4+itag)
%     h = bar([nTaskTune,nTaskTune_pass]','stacked');    
%     figXAxis([],'',[0 3],1:2,{'Behav.','Pass.'})
%     figYAxis([],'N Cells',[])
%     figAxForm
%     title(sprintf('%s Cells',tagName_unique{itag}))   
%     legend({'No Pref.','Distractor','Target'},'location','northeastoutside')
% end
% print([fnout '_tarVsDistAuROC'],'-dpdf','-fillpage')
% 
% figure
% for itag = 1:2
%     ind = ~isnan(lateStimAuROC_passive{itag}) & ...
%         (cellInfo.tag(itag).lateCycRespCells | cellInfo.tag(itag).targetRespCells);
%     taskTuning = nan(length(ind),1);
%     taskTuning(lateStimAuROCTest{itag}==0) = 3;
%     taskTuning(lateStimAuROCTest{itag}==1 & lateStimAuROC{itag} < 0.5) = 1;
%     taskTuning(lateStimAuROCTest{itag}==1 & lateStimAuROC{itag} > 0.5) = 2;
%     taskTuning_pass = nan(length(ind),1);
%     taskTuning_pass(lateStimAuROCTest_passive{itag}==0) = 3;
%     taskTuning_pass(lateStimAuROCTest_passive{itag}==1 & lateStimAuROC_passive{itag} < 0.5) = 1;
%     taskTuning_pass(lateStimAuROCTest_passive{itag}==1 & lateStimAuROC_passive{itag} > 0.5) = 2;
%     noSwitchInd = (taskTuning-taskTuning_pass) == 0;
%     
%     h = [];
%     start = 1;
%     subplot(1,2,itag)
%     h{1} = plot(lateStimAuROC{itag}(ind & noSwitchInd),lateStimAuROC_passive{itag}(ind & noSwitchInd),...
%         '.','MarkerSize',10);
%     h{1}.Color = [0.5 0.5 0.5];
%     plotInd = false(1,4);
%     plotInd(1) = true;
%     hold on
%     if sum(ind & ~noSwitchInd & taskTuning_pass == 3) > 0
%         h{2} = plot(lateStimAuROC{itag}(ind & ~noSwitchInd & taskTuning_pass == 3),...
%             lateStimAuROC_passive{itag}(ind & ~noSwitchInd & taskTuning_pass == 3),...
%             '.','MarkerSize',10);
%         plotInd(2) = true;
%     end
%     if sum(ind & ~noSwitchInd & taskTuning == 3) > 0
%         h{3} = plot(lateStimAuROC{itag}(ind & ~noSwitchInd & taskTuning == 3),...
%             lateStimAuROC_passive{itag}(ind & ~noSwitchInd & taskTuning== 3),...
%             '.','MarkerSize',10);
%         plotInd(3) = true;
%     end
%     if sum(ind & ~noSwitchInd & taskTuning ~= 3 & taskTuning_pass ~= 3) > 0
%         h{4} = plot(lateStimAuROC{itag}(ind & ~noSwitchInd & taskTuning ~= 3 & taskTuning_pass ~= 3),...
%             lateStimAuROC_passive{itag}(ind & ~noSwitchInd & taskTuning ~= 3 & taskTuning_pass ~= 3),...
%             '.','MarkerSize',10);
%         plotInd(4) = true;
%     end
%     plot([0 1],[0 1],'k--')
%     figXAxis([],'Behaving',[0 1])
%     figYAxis([],'Passive Viewing',[0 1])
%     figAxForm
%     hline(0.5,'k:')
%     vline(0.5,'k:') 
%     title(sprintf('%s Cells',tagName_unique{itag}))  
%     legendLabel = {'No Change','T/D to NP','NP to T/D','T/D to T/D'};
%     legend(legendLabel(plotInd),'location','northeastoutside')
%     
%     
% %     for icell = 1:sum(ind)
% %        h=plot(1:2,[taskTuning(icell),taskTuning_pass(icell)],'.-','MarkerSize',10);
% %        h.Color = [0 0 0];
% %        hold on
% %     end
% %     figXAxis([],'',[0 3],1:2,{'Behav.','Pass.'})
% %     figYAxis([],'Task Stim Pref.',[0 4],1:3,{'D','T','NP'})
% %     figAxForm
% %     title(sprintf('%s Cells',tagName_unique{itag}))   
% end
% print([fnout '_tarVsDistAuROC_behav2pass'],'-dpdf','-fillpage')
% 
% %%
% exExpt = 4;
% respCells = respCellsExpt(exExpt).tag(2).targetRespCells;
% exCell = randsample(find(respCells),1);
% d = antiDataExpt(exExpt).tag(2);
% % dp = antiDataExpt_passive(exExpt).tag(2);
% 
% figure
% subplot 211
% y = squeeze(d.longTC(tcStartFrame:end,exCell,:));
% for i = 1:size(y,2)
%     hold on
%     h=plot(tt_longTC,y(:,i),'-');
%     h.Color = [.5 .5 .5];
% end
% h=plot(tt_longTC,mean(y,2),'-');
% h.LineWidth = 2;
% figXAxis([],'Time From Start (ms)',[tt_longTC(1) tt_longTC(end)])
% figYAxis([],'dF/F',[])
% figAxForm([],0)
% 
% respCells = respCellsExpt(exExpt).tag(1).lateCycRespCells;
% exCell = randsample(find(respCells),1);
% d = antiDataExpt(exExpt).tag(1);
% % dp = antiDataExpt_passive(exExpt).tag(1);
% 
% subplot 212
% y = squeeze(d.longTC(tcStartFrame:end,exCell,:));
% for i = 1:size(y,2)
%     hold on
%     h=plot(tt_longTC,y(:,i),'-');
%     h.Color = [.5 .5 .5];
% end
% h=plot(tt_longTC,mean(y,2),'-');
% h.LineWidth = 2;
% figXAxis([],'Time From Start (ms)',[tt_longTC(1) tt_longTC(end)])
% figYAxis([],'dF/F',[])
% figAxForm([],0)
% %%
% tarAdaptFig = figure;
% for itag = 1:2
%     ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells;
%     cycResp = cellfun(@(x) mean(x(respwin,:),1) - mean(x(basewin_0,:),1),antiAnalysis.tag(itag).cycTC,'unif',0);
% 
%     adaptResp = cellfun(@(x) x./cycResp{1},cycResp,'unif',0);
% 
%     setFigParams4Print('landscape')
%     figure
%     if bxExpt
%         suptitle(sprintf('behav, %s, Min. & First Stim Resp. Neurons',...
%             cellInfo.tag(itag).name))
%     else
%         suptitle(sprintf('naive, %s, Min. & First Stim Resp. Neurons',...
%             cellInfo.tag(itag).name))
%     end
%     for icyc = 1:nCycles
%         subplot(2,nCycles,icyc)
%         bl = mean(antiAnalysis.tag(itag).cycTC{icyc}(basewin_0,ind),1);
%         y = antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind) - bl;
%         yerr = ste(antiAnalysis.tag(itag).cycTC{icyc}(26:end,ind),2);
%         shadedErrorBar_chooseColor(tt_cycTC,mean(y,2),yerr,[0 0 0]);
%         figXAxis([],'Time from Stim (ms)',[tt_cycTC(1) 350])
%         figYAxis([],'dF/F',cycTCLim)
%         figAxForm
%         hline(0,'k:')
%         hold on
%         vline(respWinTT,'k--');
%         if icyc == 1
%             title(sprintf('Stim #%s (%s/%s)',num2str(icyc),num2str(sum(ind)),...
%                 num2str(length(ind))))
%         else
%             title(sprintf('Stim #%s',num2str(icyc)))
%         end
%     end
% 
%     subplot 223
%     y = cellfun(@(x) mean(x(ind)),cycResp);
%     yerr = cellfun(@(x) ste(x(ind),2),cycResp);
%     errorbar(1:nCycles,y,yerr,'.')
%     figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
%     figYAxis([],'dF/F',cycTCLim)
%     figAxForm
%     title(sprintf('First Stim Resp. Cells (%s/%s)',num2str(sum(ind)),num2str(length(ind))))
%     subplot 224
%     y = cellfun(@(x) mean(x(ind)),adaptResp);
%     yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
%     errorbar(1:nCycles,y,yerr,'.')
%     figXAxis([],'Stim #',[0 nCycles+1],1:nCycles,1:nCycles)
%     figYAxis([],'Norm dF/F',[0 1.5])
%     figAxForm
%     title(sprintf('All Resp. Cells, n=%s', num2str(sum(ind))))
% 
%     print([fnout 'adaptation_anti_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
%     
%     ind = cellInfo.tag(itag).firstRespCells & cellInfo.tag(itag).minRespCells & ...
%         cellInfo.tag(itag).targetRespCells;
%     tarResp = cellfun(@(x) ...
%         mean(x(respwin_target,:),1) - mean(x(basewin_0_target,:),1),...
%         cat(2, antiAnalysis.tag(itag).cycTC{1}, ...
%         tarAnalysis.tag(itag).tc),'unif',0);
% 
%     adaptResp = cellfun(@(x) x./tarResp{1},tarResp,'unif',0);
%     
%     figure(tarAdaptFig)
%     subplot(1,2,itag)
%     y = cellfun(@(x) mean(x(ind)),adaptResp);
%     yerr = cellfun(@(x) ste(x(ind),2),adaptResp);
%     errorbar(1:3,y,yerr,'.')
%     figXAxis([],'Stim (deg)',[0 4],1:3,[0,tar2analyze])
%     figYAxis([],'Norm dF/F',[0 4])
%     figAxForm
%     title(sprintf('Target & First Stim. Resp. %s Cells, n=%s',...
%         antiAnalysis.tag(itag).name,num2str(sum(ind))))
%     print([fnout 'adaptation_target_' cellInfo.tag(itag).name],'-dpdf','-fillpage')
% end
% 
% %%
% 
% hmLim = [-0.3 0.3];
% figure
% suptitle('Anticipation: Any Task-Modulated Cells')
% colormap(brewermap([],'*RdBu'));
% for itag = 1:2
%     subplot(2,2,itag)
% %     ind = logical(cellInfo.tag(itag).lateWinRespCells);
% %     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     ind = cellInfo.tag(itag).lateWinRespCells | ...
%         cellInfo.tag(itag).lateWinSuppCells | ...
%         cellInfo.tag(itag).firstRespCells | ...
%         cellInfo.tag(itag).lateCycRespCells | ...
%         cellInfo.tag(itag).targetRespCells;
%     lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
%     lateWinResp = mean(lateWinTC(lateWinFr,:),1);
%     [~,lateWinSortInd] = sort(lateWinResp);
%     hm = flipud(lateWinTC(:,lateWinSortInd)');
%     imagesc(hm(:,tcStartFrame:end))
%     hold on
%     figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
%     figYAxis([],'Cell #',[])
%     figAxForm
%     colorbar
%     caxis(hmLim)
%     title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     if itag == 2
%         subplot(2,2,3)
%         ind = cellInfo.tag(itag).isFlex & (cellInfo.tag(itag).lateWinRespCells | ...
%             cellInfo.tag(itag).lateWinSuppCells | ...
%             cellInfo.tag(itag).firstRespCells | ...
%             cellInfo.tag(itag).lateCycRespCells);
%         lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
%         lateWinResp = mean(lateWinTC(lateWinFr,:),1);
%         [~,lateWinSortInd] = sort(lateWinResp);
%         hm = flipud(lateWinTC(:,lateWinSortInd)');
%         imagesc(hm(:,tcStartFrame:end))
%         hold on
%         figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
%         figYAxis([],'Cell #',[])
%         figAxForm
%         colorbar
%         caxis(hmLim)
%         title(sprintf('%s Cells,flex-GCaMP6s',antiAnalysis.tag(itag).name))
%         subplot(2,2,4)
%         ind = ~cellInfo.tag(itag).isFlex & (cellInfo.tag(itag).lateWinRespCells | ...
%             cellInfo.tag(itag).lateWinSuppCells | ...
%             cellInfo.tag(itag).firstRespCells | ...
%             cellInfo.tag(itag).lateCycRespCells);
%         lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
%         lateWinResp = mean(lateWinTC(lateWinFr,:),1);
%         [~,lateWinSortInd] = sort(lateWinResp);
%         hm = flipud(lateWinTC(:,lateWinSortInd)');
%         imagesc(hm(:,tcStartFrame:end))
%         hold on
%         figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
%         figYAxis([],'Cell #',[])
%         figAxForm
%         colorbar
%         caxis(hmLim)
%         title(sprintf('%s Cells,flex-tdTom',antiAnalysis.tag(itag).name))
%     end
% end
% 
% print([fnout '_anti_heatmapsAllCells'],'-dpdf','-fillpage')
% 
% 
% figure
% suptitle('Target Resp Cells')
% colormap(brewermap([],'*RdBu'));
% for itag = 1:2
%     for it = 1:2
%         if it == 1
%             offset = 0;
%         else
%             offset = 2;
%         end
%         subplot(2,2,itag+offset)
%         ind = logical(cellInfo.tag(itag).targetRespCells);
%     %     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     %     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
%         tc = tarAnalysis.tag(itag).tc{it}(:,ind);
%         if it == 1
%             resp = mean(tc(respwin,:),1);
%             [~,lateWinSortInd] = sort(resp);
%         else
%             tc1 = tarAnalysis.tag(itag).tc{it}(:,ind);
%             resp = mean(tc1(respwin,:),1);
%             [~,lateWinSortInd] = sort(resp);
%         end
%             
%         hm = flipud(tc(:,lateWinSortInd)');
%         imagesc(hm(:,tcStartFrame:end))
%         imagesc(hm(:,:))
%         hold on
%         figXAxis([],'Time from Target (ms)',[],ttLabelFr_target,ttLabel_target)
%         figYAxis([],'Cell #',[])
%         figAxForm
%         colorbar
%         caxis(hmLim)
%         title(sprintf('%s Cells, Target=%s',antiAnalysis.tag(itag).name,...
%             num2str(tar2analyze(it))))
%     end
% end
% 
% print([fnout '_tar_heatmapsTarCells'],'-dpdf','-fillpage')
% 
% %% decoding analysis
% decodeAnalysis = struct;
% for iexp = 1:nexp
%     rng(0)
%     if sum(decodeDataExpt(iexp).cellTag) == 0
%         continue
%     end
%     cellInd_notag = decodeDataExpt(iexp).cellInd & ~decodeDataExpt(iexp).cellTag;
%     cellInd_tag = decodeDataExpt(iexp).cellInd & decodeDataExpt(iexp).cellTag;
%     cellInd = false(length(cellInd_notag),1);
%     if sum(cellInd_notag) > maxCellN
%         ind = find(cellInd_notag);
%         cellSampleID = randsample(ind,maxCellN);
%         cellInd(cellSampleID) = true;
%     else
%         cellInd = cellInd_notag;
%     end
%     maxTagCellN = round(0.2*sum(cellInd));
%     if sum(cellInd_tag) > maxTagCellN
%         ind = find(cellInd_tag);
%         cellSampleID = randsample(ind,maxTagCellN);
%         cellInd(cellSampleID) = true;
%     else
%         cellInd(cellInd_tag) = true;
%     end
%     fprintf('Expt %s: %s\n',num2str(iexp),decodeDataExpt(iexp).exptName)
%     decodeAnalysis(iexp).exptName = decodeDataExpt(iexp).exptName;
%     decodeAnalysis(iexp).nCellsSelected = sum(cellInd);
%     decodeAnalysis(iexp).cellInd = cellInd;
%     decodeAnalysis(iexp).cellTag = decodeDataExpt(iexp).cellTag;
%     
%     respAllCells = zscore(decodeDataExpt(iexp).resp)';
%     trOut = decodeDataExpt(iexp).trOut;
%     
%     fprintf('Select Trials...\n')
%     trStimID = discretize(decodeDataExpt(iexp).trStim,oriBins);
%     nStimPerBin = histcounts(trStimID);
%     minBinN = min(nStimPerBin(nStimPerBin >= minTrN_mdl));
%     respStimSort = cell(1,nStimBins);
%     trOutStimSort = cell(1,nStimBins);
%     for istim = 1:nStimBins
%         ind = find(trStimID == istim);
%         if length(ind) >= minTrN_mdl
%             if istim == 1
%                 matchTrialsInd = [];
%                 if sum(nStimPerBin >= minTrN_mdl) == 2
%                     n = minBinN;
%                 elseif minBinN == nStimPerBin(istim)
%                     error('not enough FA/CR trials')
%                 else
%                     n = (nStimBins-1).*minBinN;
%                     if n > length(ind)
%                         error('not enough FA/CR trials')
%                     end
%                 end
%                 indSample = randsample(ind,n);
%                 matchTrialsInd = cat(2,matchTrialsInd,indSample);
%             else
%                 indSample = randsample(ind,minBinN);
%                 matchTrialsInd = cat(2,matchTrialsInd,indSample);
%             end
%             respStimSort{istim} = respAllCells(indSample,cellInd);
%             trOutStimSort{istim} = trOut(indSample);
%         end
%     end
%     nMatchedTrials = cumsum(cellfun(@length,trOutStimSort));
%     for istim = 1:nStimBins
%         if istim == 1
%             stimSortInd = cell(1,nStimBins);
%             stimSortInd{istim} = 1:nMatchedTrials;
%         else
%             stimSortInd{istim} = ...
%                 (nMatchedTrials(istim-1)+1):nMatchedTrials(istim);
%         end
%     end
%     
%     resp = respAllCells(matchTrialsInd,cellInd);
%     [detectTrInd, targetTrInd] = getStimAndBehaviorYs(...
%         decodeDataExpt(iexp).trOut(matchTrialsInd));
%     
%     fprintf('Run Model...\n')
%     C = eye(size(resp,2));
%     p=1;
%     [~,~,detectGLM] = glmfit(resp*C,detectTrInd,'binomial');
%     [~,~,targetGLM] = glmfit(resp*C,targetTrInd,'binomial');
% 
%     dv_detect = mean(detectTrInd);
%     dv_target = mean(targetTrInd);
% 
%     pctCorrectDetect_train = getPctCorr_trainData(detectGLM,resp,detectTrInd,dv_detect);
%     pctCorrectDetect_ho = getPctCorr_hoData(resp,detectTrInd,dv_detect);
% 
%     pctCorrectTarget_train = getPctCorr_trainData(targetGLM,resp,targetTrInd,dv_target);
%     pctCorrectTarget_ho = getPctCorr_hoData(resp,targetTrInd,dv_target);
%     
%     pctCorrDetect_xStim_train = nan(1,nStimBins);
%     pctCorrDetect_xStim_ho = nan(1,nStimBins);
%     pctCorrTarget_xStim_train = nan(1,nStimBins);
%     pctCorrTarget_xStim_ho = nan(1,nStimBins);
%     for istim = 1:nStimBins
%         if nStimPerBin(istim) >= minTrN_mdl
%             [detectStimInd, targetStimInd] = getStimAndBehaviorYs(...
%                 trOutStimSort{istim});
%             pctCorrDetect_xStim_train(istim) = getPctCorr_trainData(...
%                 detectGLM,respStimSort{istim},detectStimInd,dv_detect);
%             pctCorrDetect_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
%                 resp,detectTrInd,stimSortInd{istim},dv_detect);
%             pctCorrTarget_xStim_train(istim) = getPctCorr_trainData(...
%                 targetGLM,respStimSort{istim},targetStimInd,dv_target);
%             pctCorrTarget_xStim_ho(istim) = getPctCorr_hoData_subGroup(...
%                 resp,targetTrInd,stimSortInd{istim},dv_target);
%         end
%     end
%     
%     
%     fprintf('Add Results to Struct...\n')
%     decodeAnalysis(iexp).dvDetect = dv_detect;
%     decodeAnalysis(iexp).dvTarget = dv_target;
%     decodeAnalysis(iexp).detectWeight = detectGLM.beta(2:end);
%     decodeAnalysis(iexp).targetWeight = targetGLM.beta(2:end);
%     decodeAnalysis(iexp).pctCorrectAllDetect_train = pctCorrectDetect_train;
%     decodeAnalysis(iexp).pctCorrectAllDetect_holdout = pctCorrectDetect_ho;
%     decodeAnalysis(iexp).pctCorrectAllTarget_train = pctCorrectTarget_train;
%     decodeAnalysis(iexp).pctCorrectAllTarget_holdout = pctCorrectTarget_ho;
%     decodeAnalysis(iexp).pctCorrectXStimDetect_train = pctCorrDetect_xStim_train;
%     decodeAnalysis(iexp).pctCorrectXStimDetect_holdout = pctCorrDetect_xStim_ho;
%     decodeAnalysis(iexp).pctCorrectXStimTarget_train = pctCorrTarget_xStim_train;
%     decodeAnalysis(iexp).pctCorrectXStimTarget_holdout = pctCorrTarget_xStim_ho;
% 
% end
% ind = cellfun(@(x) ~isempty(x),{decodeAnalysis.exptName});
% decodeAnalysis = decodeAnalysis(ind);
% 
% %%
% nexp_dc = size(decodeAnalysis,2);
% targetWeight = [];
% detectWeight = [];
% cellTag = [];
% pctCorrect_target = nan(nexp_dc,4);
% pctCorrect_detect = nan(nexp_dc,4);
% for iexp = 1:nexp_dc
%     targetWeight = cat(1,targetWeight,decodeAnalysis(iexp).targetWeight);
%     detectWeight = cat(1,detectWeight,decodeAnalysis(iexp).detectWeight);
%     ind = false(1,sum(decodeAnalysis(iexp).cellInd));
%     cellTag = cat(2,cellTag,...
%         decodeAnalysis(iexp).cellTag(decodeAnalysis(iexp).cellInd));
%     
%     pctCorrect_target(iexp,1:3) = decodeAnalysis(iexp).pctCorrectXStimTarget_holdout;
%     pctCorrect_target(iexp,end) = decodeAnalysis(iexp).pctCorrectAllTarget_holdout;
%     pctCorrect_detect(iexp,1:3) = decodeAnalysis(iexp).pctCorrectXStimDetect_holdout;
%     pctCorrect_detect(iexp,end) = decodeAnalysis(iexp).pctCorrectAllDetect_holdout;
% end
% 
% binEdges = -2:0.2:2;
% figure
% subplot 221
% histogram(targetWeight(~cellTag),binEdges,'Normalization','probability')
% hold on
% histogram(targetWeight(cellTag==1),binEdges,'Normalization','probability')
% title('Stimulus Model Weights')
% figXAxis([],'Weight',[-2 2])
% figYAxis([],'Fraction of Cells',[0 0.5])
% h(1) = vline(mean(targetWeight(~cellTag)),'b');
% h(2) = vline(mean(targetWeight(cellTag==1)),'r');
% legend(h,{antiAnalysis.tag(1).name,antiAnalysis.tag(2).name},'location','northeast')
% figAxForm
% subplot 222
% histogram(detectWeight(~cellTag),binEdges,'Normalization','probability')
% hold on
% histogram(detectWeight(cellTag==1),binEdges,'Normalization','probability')
% title('Choice Model Weights')
% figXAxis([],'Weight',[-2 2])
% figYAxis([],'Fraction of Cells',[0 0.5])
% h(1) = vline(mean(detectWeight(~cellTag)),'b');
% h(2) = vline(mean(detectWeight(cellTag==1)),'r');
% legend(h,{antiAnalysis.tag(1).name,antiAnalysis.tag(2).name},'location','northeast')
% figAxForm
% stimLabel = {'0';'HT';'ET';'All'};
% subplot 223
% y = pctCorrect_target;
% yerr = ste(pctCorrect_target,1);
% plot(1:4,y,'k-');
% hold on
% errorbar(1:4,mean(y),yerr,'.')
% figXAxis([],'Stim',[0 5], 1:4, stimLabel)
% figYAxis([],'Fraction Correct',[0 1])
% hline(0.5,'k--')
% figAxForm
% title('Stimulus Model Performance')
% subplot 224
% y = pctCorrect_detect;
% yerr = ste(pctCorrect_detect,1);
% plot(1:4,y,'k-');
% hold on
% errorbar(1:4,mean(y),yerr,'.')
% figXAxis([],'Stim',[0 5], 1:4, stimLabel)
% figYAxis([],'Fraction Correct',[0 1])
% hline(0.5,'k--')
% figAxForm
% title('Choice Model Performance')
% 
% print([fnout '_decodeModelPerf&Weights'],'-dpdf','-fillpage')
% 
% %% reward analysis
% binEdges = 0:0.1:1;
% hmLim = [-0.2 0.2];
% tcLim = [-0.08 0.08];
% rewHM = figure;
% suptitle('Hard Target, Reward: Target Resp. Cells')
% colormap(brewermap([],'*RdBu'));
% rewAnal = figure;
% suptitle('Hard Target, Reward: Target Resp. Cells')
% colormap(brewermap([],'*RdBu'));
% for itag = 1:2
% %     ind = cellInfo.tag(itag).lateWinRespCells | ...
% %         cellInfo.tag(itag).lateWinSuppCells | ...
% %         cellInfo.tag(itag).firstRespCells | ...
% %         cellInfo.tag(itag).lateCycRespCells | ...
% %         cellInfo.tag(itag).targetRespCells;
%     ind = cellInfo.tag(itag).targetRespCells==1;
%     figure(rewHM)
%     subplot(2,2,itag)
%     rewTC = tarAnalysis.tag(itag).rewTC{1,1}(:,ind);
%     rewResp = mean(rewTC(rewwin,:),1);
%     [~,rewWinSortInd] = sort(rewResp);
%     hm = flipud(rewTC(:,rewWinSortInd)');
%     imagesc(hm)
%     hold on
%     figXAxis([],'Time from Reward (ms)',[],ttLabelFr_rew,ttLabel_target)
%     figYAxis([],'Cell #',[])
%     figAxForm
%     colorbar
%     caxis(hmLim)
%     title(sprintf('Hits, %s Cells',antiAnalysis.tag(itag).name))
%     subplot(2,2,itag+2)
%     rewTC = tarAnalysis.tag(itag).rewTC{1,2}(:,ind);
%     hm = flipud(rewTC(:,rewWinSortInd)');
%     imagesc(hm)
%     hold on
%     figXAxis([],'Time from Reward (ms)',[],ttLabelFr_rew,ttLabel_target)
%     figYAxis([],'Cell #',[])
%     figAxForm
%     colorbar
%     caxis(hmLim)
%     title(sprintf('Misses, %s Cells',antiAnalysis.tag(itag).name))
%     
%     figure(rewAnal)
%     subplot(2,2,itag)
%     y = tarAnalysis.tag(itag).rewTC{1,1}(:,ind);
%     yerr = ste(y,2);
%     h=shadedErrorBar_chooseColor(tt_targetTC,mean(y,2),yerr,[0 0 0]);
%     hold on
%     y = tarAnalysis.tag(itag).rewTC{1,2}(:,ind);
%     yerr = ste(y,2);
%     if ~all(isnan(yerr(:)))
%         shadedErrorBar_chooseColor(tt_targetTC,mean(y,2),yerr,[.75 0 0]);
%     end
%     figXAxis([],'Time from Reward',[tt_targetTC(1) tt_targetTC(end)],ttLabel_target,ttLabel_target)
%     figYAxis([],'dF/F',tcLim)
%     hline(0,'k--')
%     figAxForm
%     title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     
%     auroc = cellInfo.tag(itag).aurocReward;
%     subplot(2,2,3)
%     hold on
% %     if ~all(isnan(auroc(ind)))
% %         histogram(auroc(ind),binEdges,'Normalization','probability');
% %     end
% %     figXAxis([],'Reward auROC (Miss vs. Hit)',[0 1])
% %     figYAxis([],'Fraction of Cells,',[0 0.5])
% %     figAxForm
%     if ~all(isnan(auroc(ind)))
%         cdfplot(auroc(ind));
%     end
%     figXAxis([],'Reward auROC (Miss vs. Hit)',[0 1])
%     figYAxis([],'Fraction of Cells,',[0 1])
%     figAxForm
%     legend({antiAnalysis.tag(itag).name,antiAnalysis.tag(itag).name},...
%         'location','northwest')    
% end
% figure(rewHM)
% print([fnout '_reward_heatmaps'],'-dpdf','-fillpage')
% figure(rewAnal)
% print([fnout '_reward_TC'],'-dpdf','-fillpage')
% %% cluster analysis
% allCells = [];
% allCells_pass = [];
% cellTag = [];
% for itag = 1:2
%     ind = cellInfo.tag(itag).lateWinRespCells | ...
%         cellInfo.tag(itag).lateWinSuppCells | ...
%         cellInfo.tag(itag).firstRespCells | ...
%         cellInfo.tag(itag).lateCycRespCells;
%     lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
%     lateWinResp = mean(lateWinTC(lateWinFr,:),1);
%     [~,lateWinSortInd] = sort(lateWinResp);
%     lateWinTC_pass = antiAnalysis_passive.tag(itag).longTC(:,ind);
%     allCells = cat(2,allCells,lateWinTC(:,lateWinSortInd));
%     allCells_pass = cat(2,allCells_pass,lateWinTC_pass(:,lateWinSortInd));
%     cellTag = cat(2,cellTag,repmat({cellInfo.tag(itag).name},1,sum(ind)));
% end
% 
% nc = size(allCells,2);
% allCells_norm = allCells./max(allCells,[],1);
% allCells_pass_norm = allCells_pass./max(allCells,[],1);
% 
% % allCells_pca = pca(allCells_norm);
% 
% nCluster = 3;
% [clustID, clustCentroid] = kmeans(allCells_norm',nCluster,'MaxIter',1000000,'Replicates',1000);
% [~,clustID_pass] = pdist2(clustCentroid,allCells_pass_norm','euclidean','smallest',1);
% 
% figure
% subplot 221
% histogram(clustID,0:(nCluster+1),'Normalization','probability')
% figXAxis([],'Cluster #',[0 (nCluster+2)])
% figYAxis([],'Fraction of Cells',[0 1])
% figAxForm
% subplot 222
% histogram(clustID_pass,0:(nCluster+1),'Normalization','probability')
% figXAxis([],'Cluster #',[0 (nCluster+2)])
% figYAxis([],'Fraction of Cells',[0 1])
% figAxForm
% subplot 223
% for ic = 1:nCluster
%     cInd = clustID == ic;
%     hold on
%     plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
% end
% tagInd = strcmp(cellTag,'SOM+');
% plot(ones(1,sum(tagInd)).*(nCluster+1),find(tagInd),'k.','MarkerSize',10)
% figXAxis([],'Cluster #',[0 (nCluster+2)],1:(nCluster+1),...
%     cat(2,cellfun(@num2str,num2cell(1:nCluster),'unif',0),{'SOM+'}))
% figYAxis([],'Fraction of Cells',[0 nc])
% figAxForm
% subplot 224
% for ic = 1:nCluster
%     cInd = clustID_pass == ic;
%     hold on
%     plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
% end
% figXAxis([],'Cluster #',[0 (nCluster+1)])
% figYAxis([],'Fraction of Cells',[0 nc])
% figAxForm
% print([fnout 'clusterAnalysisBxPass_clusterID'],'-dpdf','-fillpage')
% 
% clusterTC = cell(1,nCluster);
% clusterTC_pass = cell(1,nCluster);
% for ic = 1:nCluster
%     cInd = clustID == ic;
%     clusterTC{ic} = allCells(:,cInd);
%     cInd = clustID_pass == ic;
%     clusterTC_pass{ic} = allCells_pass(:,cInd);    
% end
% 
% figure
% for ic = 1:nCluster
%     subplot(2,nCluster,ic)
%     y = mean(clusterTC{ic}(tcStartFrame:end,:),2);
%     yerr = ste(clusterTC{ic}(tcStartFrame:end,:),2);
%     if all(~isnan(yerr))
%         shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
%     else
%         plot(tt_longTC,y,'k')
%     end
%     figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
%     figYAxis([],'dF/F',[-0.1 0.1])
%     figAxForm
%     title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC{ic},2))))
%     
%     subplot(2,nCluster,ic+nCluster)
%     y = mean(clusterTC_pass{ic}(tcStartFrame:end,:),2);
%     yerr = ste(clusterTC_pass{ic}(tcStartFrame:end,:),2);
%     if all(~isnan(yerr))
%         shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
%     else
%         plot(tt_longTC,y,'k')
%     end
%     figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
%     figYAxis([],'dF/F',[-0.1 0.1])
%     figAxForm
%     title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC_pass{ic},2))))
% end
% print([fnout 'clusterAnalysisBxPass_clusterTC'],'-dpdf','-fillpage')
% 
% switchClustID = clustID'-clustID_pass ~= 0;
% 
% figure
% subplot 121
% y = sum(switchClustID & ~tagInd)./sum(~tagInd);
% plot(1,y,'.','MarkerSize',20)
% hold on
% y = sum(switchClustID & tagInd)./sum(tagInd);
% plot(2,y,'.','MarkerSize',20);
% figXAxis([],'',[0 3],1:2,{'SOM-','SOM+'})
% figYAxis([],'Fraction Switch Cluster',[0 0.2])
% figAxForm
% subplot 122
% hold on
% for ic = 1:nCluster
%     cInd = clustID' == ic & switchClustID;
%     n = histcounts(clustID_pass(cInd),1:nCluster+1);
%     scatter(ic.*(ones(1,nCluster)),1:nCluster,(n*5)+1)
% end
% figXAxis([],'Behavior Cluster #',[0 nCluster+1])
% figYAxis([],'Passive Cluster #',[0 nCluster+1])
% figAxForm
% title(sprintf('Cluster Switch Cells, n=%s',num2str(sum(switchClustID))))
% print([fnout 'clusterAnalysisBxPass_nClusterSwitch'],'-dpdf','-fillpage')
% 
% %%
%     figure
%     suptitle('Target Resp.')
% for itag = 1:2
%     ind = logical(cellInfo.tag(itag).targetRespCells);
%     tarResp_bx = mean(tarAnalysis.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis.tag(itag).easyTC(basewin_0,:),1);
%     tarResp_pass = mean(tarAnalysis_passive.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis_passive.tag(itag).easyTC(basewin_0,:),1);
% 
%     tarRespLim = [-0.01 0.12];
% %     for icyc = 1:nCycles
%         subplot(1,2,itag)
%         x = tarResp_bx(ind);
%         y = tarResp_pass(ind);
%         plot(x,y,'.','MarkerSize',10)
%         hold on
%         plot(tarRespLim,tarRespLim,'k--')
%         figXAxis([],'Behavior',tarRespLim)
%         figYAxis([],'Passive',tarRespLim)
%         figAxForm
%         title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
% %     end
% end
% 
% figure
%     suptitle('Late Resp.')
% for itag = 1:2
%     ind = logical(cellInfo.tag(itag).lateCycRespCells);
%     lateResp_bx = mean(antiAnalysis.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis.tag(itag).lateCycTC(basewin_0,:),1);
%     lateResp_pass = mean(antiAnalysis_passive.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis_passive.tag(itag).lateCycTC(basewin_0,:),1);
% 
%     lateRespLim = [-0.01 0.1];
% %     for icyc = 1:nCycles
%         subplot(1,2,itag)
%         x = lateResp_bx(ind);
%         y = lateResp_pass(ind);
%         plot(x,y,'.','MarkerSize',10)
%         hold on
%         plot(lateRespLim,lateRespLim,'k--')
%         figXAxis([],'Behavior',lateRespLim)
%         figYAxis([],'Passive',lateRespLim)
%         figAxForm
%         title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
% %     end
% end
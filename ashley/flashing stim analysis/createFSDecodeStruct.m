function decodeDataExpt = createFSDecodeStruct(datasetStr,cellsOrDendrites,doAuditory,oriBinSize,ampBinSize);
%cellsOrDendrites: 1 == cells; 2 == dendrites
ds = datasetStr;
rc = behavConstsAV;
imgParams_FSHVA
bxParams_FSAV_attnV1ms

eval(ds)

mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if expt(1).folder(1:5) == 'ashle'
    titleStr = ds(6:end);
    if cellsOrDendrites == 1
        load(fullfile(rc.caOutputDir,ds,...
            [mouse_str '_trOutcomeStruct_cells' ds(5:end) '.mat']));
        fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeStruct_cells' ds(5:end) '.mat']); 
    elseif cellsOrDendrites == 2
        load(fullfile(rc.caOutputDir,ds,...
            [mouse_str '_trOutcomeStruct_dend' ds(5:end) '.mat']));
        fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeStruct_dend' ds(5:end) '.mat']);
    end
elseif expt(1).folder(1:5) == 'linds'
    titleStr = ds;
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_trOutcomeStruct_cells_' ds '.mat']));
    fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeStruct_cells_' ds '.mat']); 
end

nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nexp = size(expt,2);
if doAuditory
    nav = 2;
else
    nav = 1;
end

if oriBinSize>0
    orientations = 0:oriBinSize:(180-oriBinSize);
    oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
    nOri = length(orientations);
end

decodeDataExpt = struct;
decodeDataExpt.av(visualTrials).name = 'Visual';
if doAuditory
    decodeDataExpt.av(auditoryTrials).name = 'Auditory';
end
nTrialsPerExpt = nan(1,nexp);

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if imouse == 1 && iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end

        d = mouse(imouse).expt(iexp);
        exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);
        
        nTrialsPerExpt(exptN) = length(d.av(visualTrials).align(1).outcome) + ...
                length(d.av(auditoryTrials).align(1).outcome);

        cycLengthFr = d.info.cycTimeFrames;
        
        maxCycles = max(cat(2,d.av(visualTrials).align(alignStart).nCycles,...
                d.av(auditoryTrials).align(alignStart).nCycles));
            
        decodeDataExpt(exptN).exptArea = mouse(imouse).expt(iexp).area;
        area_ind = find(strcmp(area_list,mouse(imouse).expt(iexp).area));
        base_d = basewin_0{area_ind};
        base_t = basewin_0_target{area_ind};
        resp_d = respwin{area_ind};
        resp_t = respwin_target{area_ind};
            
        for iav = 1:nav
            dd = d.av(iav);
            trOut = [];
            trStim = [];
            trResp = [];
            trTC = [];
%             figure;
            for ialign = 2:4
                ddd = dd.align(ialign);
                tc = ddd.respTC;
                if ialign == 4
%                     subplot(2,1,ialign-2)
%                     plot(mean(mean(tc,2),3))
%                     hold on
%                     vline(base_t,':k')
%                     vline(resp_t)
%                     suptitle(['mouse ' num2str(imouse) '- expt ' num2str(iexp) '- area ' mouse(imouse).expt(iexp).area])
                    tc_temp = circshift(tc,-1,1);
                    tc_temp(end,:,:) = nan;
                    trTC = cat(3,trTC,tc);
                    trResp = cat(2,trResp,squeeze(mean(tc(resp_t,:,:),1)...
                        - mean(tc(base_t,:,:),1)));
                else
%                     if ialign==3
%                     subplot(2,1,ialign-2)
%                     plot(mean(mean(tc,2),3))
%                     hold on
%                     vline(base_d,':k')
%                     vline(resp_d)
%                     end
                    trResp = cat(2,trResp,squeeze(mean(tc(resp_d,:,:),1)...
                        - mean(tc(base_d,:,:),1)));
                    trTC = cat(3,trTC,circshift(tc,-1,1)- mean(tc(base_d,:,:),1));
                end
                if ialign == alignCR
                    trOut = cat(2,trOut,repmat({'cr'},[1,length(ddd.outcome)]));
                else
                    trOut = cat(2,trOut,ddd.outcome);
                end
                if isempty(ddd.ori) & isempty(ddd.amp)
                    trStim = cat(2,trStim,zeros(1,length(ddd.outcome)));
                elseif iav == 1
                    trStim = cat(2,trStim,ddd.ori);
                elseif iav == 2
                    trStim = cat(2,trStim,ddd.amp);
                end
            end
            trOut(strcmp(trOut,'success')) = {'h'};
            trOut(strcmp(trOut,'ignore')) = {'m'};
            trOut(strcmp(trOut,'failure')) = {'fa'};
            
            decodeDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            decodeDataExpt(exptN).mouseName = [mouse(imouse).mouse_name];
            decodeDataExpt(exptN).exptDate = [d.date];
            decodeDataExpt(exptN).av(iav).outcome = trOut;
            decodeDataExpt(exptN).av(iav).stim = trStim;
            decodeDataExpt(exptN).av(iav).resp = trResp;
            decodeDataExpt(exptN).av(iav).tc = trTC(1:(nBaselineFr*2),:,:);

            if iav == 1 & oriBinSize>0
                binEdges = oriBins;
                trStimID = discretize(trStim,binEdges);
            elseif iav == 1
                [oris x trStimID] = unique(trStim);
            elseif iav == 2
                binEdges = ampBins;
            end
            
            rng(0)
            
            if size(d.av(visualTrials).align,2) > 4
                r = squeeze(mean(...
                    d.av(visualTrials).align(5).respTC(respwin_target,:,:),1) - ...
                    mean(d.av(visualTrials).align(5).respTC(basewin_0_target,:,:),1));
                trOut =  d.av(visualTrials).align(5).outcome;
                ind = cellfun(@(x) ~isempty(x),trOut);
                trOut = trOut(ind);
                r = r(:,ind);
                trOut(strcmp(trOut,'FA')) = {'h'};
                trOut(strcmp(trOut,'CR')) = {'m'};
                decodeDataExpt(exptN).av(visualTrials).catchResp = r;
                decodeDataExpt(exptN).av(visualTrials).catchOutcome = trOut;
                decodeDataExpt(exptN).av(visualTrials).catchStim = ...
                    d.av(visualTrials).align(5).ori;
            end
            
            do = d.oriTuning;
            if size(do.oriResp,2) == 8
                decodeDataExpt(exptN).oriResp = do.oriResp(:,1:2:8);
                decodeDataExpt(exptN).oriRespErr = do.oriRespSem(:,1:2:8);
            elseif size(do.oriResp,2) ~= 4
                error('error in orientations used for passivie tuning')
            else
                decodeDataExpt(exptN).oriResp = do.oriResp;
                decodeDataExpt(exptN).oriRespErr = do.oriRespSem;
            end
            decodeDataExpt(exptN).fit = do.oriFit;
            decodeDataExpt(exptN).isTuned = do.oriFitReliability < tuningReliabilityThresh;
            [~,decodeDataExpt(exptN).fitPeak] = max(do.oriFit,[],1);
%             [~,oriPref] = histc(oriTuningExpt(exptN).fitPeak,oriBinEdges);
%             oriPref(oriPref == length(orientations)+1 | oriPref == length(oriBinEdges) == 1) = 1;
%             oriTuningExpt(exptN).oriPref = oriPref;
            decodeDataExpt(exptN).tuningReliability = do.oriFitReliability;
        end
    end
end

save(fnout, 'decodeDataExpt');

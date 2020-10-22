function FSData = createFSTuningDataStruct(ds);

rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms
eval(ds)

mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

if expt(1).folder(1:5) == 'ashle'
titleStr = ds(6:end);
if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_decodeStruct_cells' ds(5:end) '.mat']));
    fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeAnalysis_cells' ds(5:end) '.mat']); 
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,ds,...
        [mouse_str '_decodeStruct_dend' ds(5:end) '.mat']));
    fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_decodeAnalysis_dend' ds(5:end) '.mat']);
end
elseif expt(1).folder(1:5) == 'linds'
titleStr = ds;
load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_decodeStruct_cells_' ds '.mat']));
fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_FSData_cells_' ds '.mat']); 
end

nexp = size(expt,2);
for iexp = 1:nexp
    fprintf([decodeDataExpt(iexp).exptDate '_i' decodeDataExpt(iexp).mouseName '\n'])
    nCells = size(decodeDataExpt(iexp).av(visualTrials).resp,1);
    [oris x trStimID] = unique(decodeDataExpt(iexp).av(visualTrials).stim);
    nstim = length(oris);
    FSData(iexp).sigRespInd = zeros(nCells,1);
    for istim = 1:nstim
        ind = find(trStimID == istim);
        stimResp{istim} = decodeDataExpt(iexp).av(visualTrials).resp(:,ind);
        if istim == 1
            baseResp{istim} = squeeze(mean(decodeDataExpt(iexp).av(visualTrials).tc(basewin_0,:,ind),1));
        else
            baseResp{istim} = squeeze(mean(decodeDataExpt(iexp).av(visualTrials).tc(basewin_0_target,:,ind),1));
        end
        [hResp{istim} pResp{istim}] = ttest(baseResp{istim},stimResp{istim},'tail','left','dim',2,'alpha',0.05./(nstim-1));

        if istim>1
            for iCell = 1:nCells
                FSData(iexp).auROC(iCell,istim) = roc_gh(stimResp{1}(iCell,:),stimResp{istim}(iCell,:));
            end
        end
        FSData(iexp).sigRespInd = FSData(iexp).sigRespInd + hResp{istim};
    end
    FSData(iexp).oris = oris;
end
save(fnout,'FSData')

function oriData = createOriDataStruct(ds, tuningThresh);

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
    fnout = fullfile(rc.caOutputDir,ds,[mouse_str '_tuningData_cells_' ds '.mat']); 
end
nexp = size(expt,2);
for iexp = 1:nexp
    fprintf([decodeDataExpt(iexp).exptDate '_i' decodeDataExpt(iexp).mouseName '\n'])
    for i = 1:nexp
        if strcmp(expt(i).date,decodeDataExpt(iexp).exptDate)
            if strcmp(expt(i).mouse,decodeDataExpt(iexp).mouseName)
                break
            end
        end
    end
    dir_str = ['runs-' expt(i).dirtuning];
    load(fullfile(rc.lindseyAnalysis, [decodeDataExpt(iexp).exptDate '_i' decodeDataExpt(iexp).mouseName], [decodeDataExpt(iexp).exptDate '_i' decodeDataExpt(iexp).mouseName '_' dir_str], [decodeDataExpt(iexp).exptDate '_i' decodeDataExpt(iexp).mouseName '_' dir_str '_oriTuningAndFits.mat']))
    cellInd = fitReliability <= tuningThresh;
    ind = find(cellInd);
    oriData(iexp).OSI = zeros(1,length(ind));
    oriData(iexp).prefOri = zeros(1,length(ind));
    for iCell = 1:length(ind)
        i = ind(iCell);
        [max_val max_ind] = max(vonMisesFitAllCellsAllBoots(:,1,i),[],1);
        oriData(iexp).prefOri(:,iCell) = max_ind-1;
        if max_ind<91
            min_ind = max_ind+90;
        else
            min_ind = max_ind-90;
        end
        min_val = vonMisesFitAllCellsAllBoots(min_ind,1,i);
        if min_val < 0
            min_val = 0;
        end
        oriData(iexp).OSI(:,iCell) = (max_val-min_val)./(max_val+min_val);
    end
    oriData(iexp).nCells = length(cellInd);
    oriData(iexp).nWellfitCells = length(ind);
    oriData(iexp).wellfitCells = cellInd;
end

save(fnout,'oriData')

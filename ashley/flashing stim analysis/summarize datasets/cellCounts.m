clear all
close all
datasetStr = '_V1';
cellsOnly = 0;
% load data
targetResp = 13;
base1Resp = 12;
baseSustResp = 8;
av = behavParamsAV_naive;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end

rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
if cellsOnly
load(fullfile(rc.caOutputDir,dataGroup, 'cells only',[mouse_str '_CaSummary' datasetStr '.mat']));
fnout = fullfile(rc.caOutputDir, dataGroup, 'cells only', [titleStr '_' mouse_str]); %% maybe lose mouse_str
else
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' mouse_str]); %% maybe lose mouse_str
end

%% get vectors of cell indices for all task-responsive groups

tarCellsInd = [];
b1CellsInd = [];
bSCellsInd = [];
cc = 0;

tarCellsInd_inv = [];
b1CellsInd_inv = [];
bSCellsInd_inv = [];
cc_inv = 0;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        tar = mouse(imouse).expt(iexp).cells(targetResp).ind;
        b1 = mouse(imouse).expt(iexp).cells(base1Resp).ind;
        bS = mouse(imouse).expt(iexp).cells(baseSustResp).ind;
        tarCellsInd = cat(1,tarCellsInd,tar+cc);
        b1CellsInd = cat(1,b1CellsInd,b1+cc);
        bSCellsInd = cat(1,bSCellsInd,bS+cc);
        cc = mouse(imouse).expt(iexp).info.nCells+cc;
        
        if mouse(imouse).expt(iexp).info.isCatch
        tar = mouse(imouse).expt(iexp).cells(targetResp).ind;
        b1 = mouse(imouse).expt(iexp).cells(base1Resp).ind;
        bS = mouse(imouse).expt(iexp).cells(baseSustResp).ind;
        tarCellsInd_inv = cat(1,tarCellsInd_inv,tar+cc_inv);
        b1CellsInd_inv = cat(1,b1,b1CellsInd_inv,b1+cc_inv);
        bSCellsInd_inv = cat(1,bSCellsInd_inv,bS+cc_inv);
        cc_inv = mouse(imouse).expt(iexp).info.nCells+cc_inv;
        end
    end
end

allCells = 1:cc;
allCells_inv = 1:cc_inv;

bAllCells = unique(cat(1,b1CellsInd,bSCellsInd));
overlapCells = intersect(tarCellsInd,bAllCells);
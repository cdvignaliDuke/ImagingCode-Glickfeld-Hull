%% dF/F by trial
tStartInd = double(1:nON+nOFF:size(data_reg,3));
dataF = mean(reshape(data_reg(:,:,linspaceNDim(tStartInd,tStartInd+(nOFF-1),nOFF)),size(data_reg,1),size(data_reg,2),nTrials,nOFF),4);

dataDF = zeros(size(data_reg));
dFoverF = zeros(size(data_reg));
for itrial = 1:nTrials
    dataDF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)) = bsxfun(@minus,data_reg(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)),dataF(:,:,itrial));
    dFoverF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)) = bsxfun(@rdivide,dataDF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)),dataF(:,:,itrial));
end

writetiff(dFoverF,'dFoverFmov');

%%
if VSR == 2

%% for direction tuning exp
    % sort like-trials in visstimret experiment

tDirection = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(tDirection);
Oris = Dirs(1:length(Dirs)/2);


dirInd = [];
tStartDirInd = [];
dFoverF_dirsort = [];
dFoverF_dirmean = zeros(size(dFoverF,1),size(dFoverF,2),length(Dirs),nON+nOFF);

for i = 1:length(Dirs)
    dirInd{i} = find(tDirection == Dirs(i));
    tStartDirInd{i} = linspaceNDim(tStartInd(dirInd{i}),tStartInd(dirInd{i})+(nON+nOFF-1),nON+nOFF);
    tempind = tStartDirInd{i};
    dFoverF_dirsort{i} = reshape(dFoverF(:,:,reshape(tempind,1,numel(tempind))),size(dFoverF,1),size(dFoverF,2),size(tempind,1),size(tempind,2));
    dFoverF_dirmean(:,:,i,:) = squeeze(mean(dFoverF_dirsort{i},3));
end

dFoverF_orimean = zeros(size(dFoverF,1),size(dFoverF,2),length(Oris),nON+nOFF);
for i = 1:length(Oris)
    dFoverF_orimean(:,:,i,:) = squeeze(mean(cat(3,dFoverF_dirsort{i},dFoverF_dirsort{i+length(Oris)}),3));
end

%% create first response stack for each dir, save as tif

firstRespStack = dFoverF_dirmean(:,:,:,nOFF+2);
firstRespStackOri = dFoverF_orimean(:,:,:,nOFF+2);
dirRespStack = mean(dFoverF_dirmean(:,:,:,nOFF+1:end),4);
oriRespStack = mean(dFoverF_orimean(:,:,:,nOFF+1:end),4);
try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end
writetiff(firstRespStack,'firstRespStack');
writetiff(firstRespStackOri,'firstRespStackOri');
writetiff(oriRespStack,'oriRespStack');
writetiff(dirRespStack,'dirRespStack');

%% for retinotopy exp
elseif VSR == 1
    
tAz = double(cell2mat(input.tGratingAzimuthDeg));
[nAz azPos] = histc(tAz,unique(tAz));
tEl = double(cell2mat(input.tGratingElevationDeg));
[nEl elPos] = histc(tEl,unique(tEl));
if any(tEl == 0) | any(tAz == 0)
    if any(tEl == 0)
        tEl2 = tEl+1;
    else
        tEl2 = tEl;
    end
    if any(tAz == 0)
        tAz2 = tAz+1;
    else
        tAz2 = tAz;
    end
    pos = cart2pol(tAz2,tEl2);
else
    pos = cart2pol(tAz,tEl);
end
[n posN] = histc(pos,unique(pos));
Rets = unique(posN);



retInd = [];
tStartRetInd = [];
dFoverF_retsort = [];
dFoverF_retmean = zeros(size(dFoverF,1),size(dFoverF,2),length(Rets),nON+nOFF);

for i = 1:length(Rets)
    retInd{i} = find(posN == Rets(i));
    tStartRetInd{i} = linspaceNDim(tStartInd(retInd{i}),tStartInd(retInd{i})+(nON+nOFF-1),nON+nOFF);
    tempind = tStartRetInd{i};
    dFoverF_retsort{i} = reshape(dFoverF(:,:,reshape(tempind,1,numel(tempind))),size(dFoverF,1),size(dFoverF,2),size(tempind,1),size(tempind,2));
    dFoverF_retmean(:,:,i,:) = squeeze(mean(dFoverF_retsort{i},3));
end 

%% writetiffs
firstRespStack_ret = dFoverF_retmean(:,:,:,nOFF+2);
retRespStack = mean(dFoverF_retmean(:,:,:,nOFF+1:end),4);
try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end
writetiff(firstRespStack_ret,'firstRespStack');
writetiff(retRespStack,'retRespStack');

%% for moving dots directions
elseif VSR == 3
    
   tDirection = cell2mat(input.tDotDirectionDeg);
Dirs = unique(tDirection);


dirInd = [];
tStartDirInd = [];
dFoverF_dirsort = [];
dFoverF_dirmean = zeros(size(dFoverF,1),size(dFoverF,2),length(Dirs),nON+nOFF);

for i = 1:length(Dirs)
    dirInd{i} = find(tDirection == Dirs(i));
    tStartDirInd{i} = linspaceNDim(tStartInd(dirInd{i}),tStartInd(dirInd{i})+(nON+nOFF-1),nON+nOFF);
    tempind = tStartDirInd{i};
    dFoverF_dirsort{i} = reshape(dFoverF(:,:,reshape(tempind,1,numel(tempind))),size(dFoverF,1),size(dFoverF,2),size(tempind,1),size(tempind,2));
    dFoverF_dirmean(:,:,i,:) = squeeze(mean(dFoverF_dirsort{i},3));
end


%% create first response stack for each dir, save as tif

firstRespStack = dFoverF_dirmean(:,:,:,nOFF+2);
dirRespStack = mean(dFoverF_dirmean(:,:,:,nOFF+1:end),4);
try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end
writetiff(firstRespStack,'firstRespStack');
writetiff(dirRespStack,'dirRespStack');
 
end
    



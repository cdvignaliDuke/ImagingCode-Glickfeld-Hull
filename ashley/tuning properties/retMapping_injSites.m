% experiment specs
expt(1).SubNum = '617';
expt(1).mouse = 'AW17';
expt(1).date = '160511';
expt(1).dateFolder = '160511';
expt(1).img_loc  = {'V1';'AL'};
expt(1).visstim  = {'drifting';'drifting'};
expt(1).time_mat = ['1334'; '1345'];
expt(1).runs = ['001'; '002'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';

%%
SubNum = expt(1).SubNum;
mouse = expt(1).mouse;
date = expt(1).date;

sites = unique(expt(1).img_loc);
visstim = unique(expt(1).visstim);

down = 5;



%%
for irun = 1:expt(1).nrun
    disp(irun)
    %% load the data
    time = expt(1).time_mat(irun,:);
    ImgFolder = expt(1).runs(irun,:);

    fName = [ImgFolder '_000_000'];

    mworks = ['data-' 'i' SubNum '-' date '-' time]; 
    load (['Y:\home\andrew\Behavior\Data\' mworks]);

    CD = ['Z:\data\' mouse '\two-photon imaging\' expt(1).dateFolder '\' ImgFolder];
    cd(CD);
    imgMatFile = [fName '.mat'];
    load(imgMatFile);
    nframes = info.config.frames;
    tic
    data = sbxread(fName,0,nframes);
    toc
    data = squeeze(data);    
%%
    data_down = stackGroupProject(data,down);
    clear data
    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    data = data_sub;
    clear data_sub
    % register
    data_avg = mean(data(:,:,100:110),3);
    figure; imagesq(data_avg); colormap(gray)
    [out data_reg] = stackRegister(data, data_avg);
    clear data
    
    %%
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
    
    %%
    writetiff(data_reg, 'Retinotopy');
    clear data_reg
end
%% if data already registered
for irun = 1:expt(1).nrun
    time = expt(1).time_mat(irun,:);
    ImgFolder = expt(1).runs(irun,:);
    mworks = ['data-' 'i' SubNum '-' date '-' time]; 
    load (['Y:\home\andrew\Behavior\Data\' mworks]);
    data_reg = readtiff(fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder, 'Retinotopy.tif'));
    %%    
    nON = double(input.nScansOn/down);
    nOFF = double(input.nScansOff/down);
    nStim = double(input.gratingElevationStepN.*input.gratingAzimuthStepN);
    nRep = floor(size(data_reg,3)/((nON+nOFF)*nStim));
    data_reg = data_reg(:,:,1:(((nON+nOFF)*nStim)*nRep));
    nTrials = (nStim.*nRep);

    %% create dF/F stack, mask_cell, and timecourses

    nOFF_ind = zeros(1,(nOFF*nStim*nRep));
    start = 1;
    for iStim = 1:(nRep*nStim)
        nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
        start = start+nOFF;
    end

    nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
    nON_avg = mean(data_reg(:,:,nON_ind),3);
    nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

    % dF average F
    dF_data = bsxfun(@minus,data_reg, nOFF_avg);
    dFoverF = bsxfun(@rdivide,dF_data,nOFF_avg);

    max_dF = max(dF_data,[],3);
    maxDFoverF = max(dFoverF,[],3);
    figure; imagesq(maxDFoverF); colormap(gray)
    
    %save max DF/F
    writetiff(maxDFoverF, 'maxDFoverF');

    bwout = imCellEditInteractive(maxDFoverF);
    mask_cell = bwlabel(bwout);

    data_TC = stackGetTimeCourses(data_reg,mask_cell);
    figure; tcOffsetPlot(data_TC)
    
    %% subtract neuropil from timecourses
    buf = 4;
    np = 6;
    nCells = size(data_TC,2);
    neuropil = imCellNeuropil(mask_cell,buf,np);

    npTC = zeros(size(data_TC));
    for i = 1:nCells
        tempNPmask = squeeze(neuropil(:,:,i));
        if sum(sum(tempNPmask)) > 0
        npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
        end
    end

    data_TC_mavg = tsmovavg(data_TC,'s',10,1);
    npTC_mavg = tsmovavg(npTC,'s',10,1);
    
    ii= 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for i = 1:100
        x(i,:) = skewness(data_TC_mavg-tcRemoveDC(npTC_mavg*ii(i)));
    end
    [max_skew ind] =  max(x,[],1);
    np_w = 0.01*ind;
    data_TCsub = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);
    data_TC = data_TCsub;
    clear data_reg
    %% save
    fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fileSave);
    save('mask&TCDir.mat','neuropil','mask_cell','data_TC');
    %% 
    VSsize = input.gratingDiameterDeg;
    VSdirectionDeg = input.gratingDirectionDeg;
    tAz = double(cell2mat(input.tGratingAzimuthDeg));
    Az = unique(tAz);
    [nAz azPos] = histc(tAz,Az);
    tEl = double(cell2mat(input.tGratingElevationDeg));
    El = unique(tEl);
    [nEl elPos] = histc(tEl,El);
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
    Rets = unique(pos);

    AzPos = NaN(1,nStim);
    ElPos = NaN(1,nStim);
    for i = 1:nStim
        AzPos(i) = unique(tAz(pos == Rets(i)));
        ElPos(i) = unique(tEl(pos == Rets(i)));
    end

%% dF/F by trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);
stimON_ind = nOFF+1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end


%% sort data by trial type
dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanRetResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
errbarRets = zeros(size(data_TC,2),nStim);
for i = 1:nStim
    trials = find(pos(:,1:nTrials) == Rets(i));
    dFoverF_meanRetResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
    errbarRets(:,i) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end

figure;
for i = 1:nStim
    plot(dFoverF_meanRetResp(:,3,i));
    hold on
end

%% plot tuning curves
dFoverF_meanOFFRetResp = (squeeze(mean(dFoverF_meanRetResp(1:10,:,:),1)));

RetRespPerCell = (squeeze(mean(dFoverF_meanRetResp(11:end,:,:),1)));

% az x el matrix per cell
resp = squeeze(mean(dFoverFCellsTrials(11:end,:,:),1));
resp_AzElMat = zeros(size(resp,1),length(El),length(Az));
errResp_AzElMat = zeros(size(resp,1),length(El),length(Az));

for iaz = 1:length(Az)
    azind = intersect(find(tAz == Az(iaz)),1:nTrials);
    for iel = 1:length(El)
        elind = intersect(find(tEl == El(iel)),1:nTrials);
        resp_AzElMat(:,iel,iaz) = mean(resp(:,intersect(azind,elind)),2);
        errResp_AzElMat = squeeze(std(resp(:,intersect(azind,elind)),[],2)/sqrt(double(length(intersect(azind,elind)))));
    end
end

%%
figure;
imagesc(squeeze(mean(resp_AzElMat,1)));
set(gca,'xTick', [1:length(Az)])
set(gca,'xTickLabel', cellfun(@num2str,num2cell(Az),'UniformOutput',false))
set(gca,'yTick', [1:length(El)])
set(gca,'yTickLabel', cellfun(@num2str,num2cell(El),'UniformOutput',false))
xlabel('Az')
ylabel('El')
colorbar

title({'Retinotopy';strjoin([expt(1).img_loc(irun) '; ' expt(1).visstim(irun) ' gratings']) })

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

% get position used for behavior


print(fullfile('Z:\analysis\',mouse,'two-photon imaging',['RetPreferences' strjoin(['-' expt(1).img_loc(irun) '_' expt(1).visstim(irun) ' gratings'])] ), '-dpdf')
end

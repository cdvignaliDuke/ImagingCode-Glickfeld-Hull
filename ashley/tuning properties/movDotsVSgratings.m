% load dataset info
awMovingDotsDirectionTuning_V1

iexp = 1;
runs = expt(iexp).runs;
nrun = expt(iexp).nrun;
runs_mat = runs(1,:);
for i = 2:nrun
    runs_mat = [runs_mat '-' runs(i,:)];
end
imgsiz = expt(iexp).imgsiz;

%% register each dataset
for i = 1:nrun

    if strcmp(expt(iexp).regImg,expt(iexp).runs(i,:)) == 0;

        imgMatFile = [expt(iexp).runs(i,:) '_000_000.mat'];
        cd(fullfile('Z:\data\',expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, expt(iexp).runs(i,:)));
        load(imgMatFile);

        nframes = info.config.frames;
        tic
        data = sbxread([expt(iexp).runs(i,:) '_000_000'],0,nframes);
        toc
        data = squeeze(data);
        
        mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).time_mat(i,:)]; 
        load (fullfile('Y:\home\andrew\Behavior\Data',mworks));
        
        down = 10;
        nON = double(input.nScansOn)./down;
        nOFF = double(input.nScansOff)./down;
        nStim = double(input.gratingDirectionStepN);

        %average signals in time
        data_down = stackGroupProject(data,down);
        clear data

        %remove negative data (by addition)
        data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
        clear data_down

        % load image for registration
        load(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, expt(iexp).regImg,'regImg.mat'))

        [out data_reg] = stackRegister(data_sub, data_avg);
        clear data_sub
        
        cd(fullfile('Z:\analysis\',expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, expt(iexp).runs(i,:)))
        if strcmp('DirectionTuning_V1',expt(iexp).data_reg_name{i}) == 1
            writetiff(data_reg, 'DirectionTuning_V1');
        else
            writetiff(data_reg, 'DirectionTuning_movDots');
        end
        
        nRep = size(data_reg,3)./((nON+nOFF)*nStim);
        nTrials = (nStim.*nRep);        
        
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
        writetiff(maxDFoverF, 'maxDFoverF');
    end
    clear data_reg input
end
%% create a max dF/F stack from all experiments done, use max of that to find cells to analyze
    maxDFstack = zeros(imgsiz(1),imgsiz(2),nrun);
for i = 1:nrun
    maxDFstack(:,:,i) = readtiff(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, expt(iexp).runs(i,:),'maxDFoverF.tif'));
end

maxDFoverF = max(maxDFstack,[],3);
figure; imagesq(maxDFoverF); colormap(gray)

bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);

%% get timecourses for each experiment with this cell mask, subtract neuropil

nCells = length(unique(mask_cell))-1;
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);

for i = 1:nrun
    data_reg = readtiff(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, expt(iexp).runs(i,:),[expt(iexp).data_reg_name{i} '.tif']));
    data_TC_temp = stackGetTimeCourses(data_reg,mask_cell);
    
    npTC = zeros(size(data_TC_temp));
    for iNP = 1:size(data_TC_temp,2)
        tempNPmask = squeeze(neuropil(:,:,iNP));
        if sum(sum(tempNPmask)) > 0
        npTC(:,iNP) = stackGetTimeCourses(data_reg,tempNPmask);
        end
    end
        %get weights by maximizing skew
    ii= 0.01:0.01:1;
    x = zeros(length(ii), nCells);
    for iNP = 1:100
        x(iNP,:) = skewness(data_TC_temp-tcRemoveDC(npTC*ii(iNP)));
    end
    [max_skew ind] =  max(x,[],1);
    skew(buf,:) = max_skew;
    np_w = 0.01*ind;
    npSubTC = data_TC_temp-bsxfun(@times,tcRemoveDC(npTC),np_w);

    data_TC{i} = npSubTC;
    clear data_reg npSubTC
end

%% save data
try
    save(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, runs_mat,'dataTCstruct.mat'),'data_TC','maxDFoverF','maxDFstack','mask_cell','neuropil');
catch
    cd(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date));
    mkdir(runs_mat);
    save(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, runs_mat,'dataTCstruct.mat'),'data_TC','maxDFoverF','maxDFstack','mask_cell','neuropil');
end

%% load saved data
% load dataset info
awMovingDotsDirectionTuning_V1

iexp = 1;
runs = expt(iexp).runs;
nrun = expt(iexp).nrun;
runs_mat = runs(1,:);
for i = 2:nrun
    runs_mat = [runs_mat '-' runs(i,:)];
end
imgsiz = expt(iexp).imgsiz;

load(fullfile('Z:\Analysis\', expt(iexp).mouse, expt(iexp).folder, expt(iexp).date, runs_mat,'dataTCstruct.mat'));

%% dF/F
down = expt(iexp).downSampleRate;
nON = expt(iexp).framesONframesOFF(1)/down;
nOFF = expt(iexp).framesONframesOFF(2)/down;

dFoverF = {};
nStim = {};
nRep = {};
nTrials = {};
for i = 1:nrun
    
    mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).time_mat(i,:)]; 
    load (fullfile('Y:\home\andrew\Behavior\Data', mworks));
    
    stimOFF_ind{i} = 1:nOFF+nON:size(data_TC{i},1);
    stimON_ind{i} = nOFF+1:nOFF+nON:size(data_TC{i},1);
    nStim{i} = double(input.dotDirectionStepN);
    nRep{i} = size(data_TC{i},1)./((nON+nOFF)*nStim{i});
    nTrials{i} = (nStim{i}.*nRep{i});
    if strcmp(expt(iexp).dataset_visstim{i},'moving dots') == 1
        DirectionDeg{i} = cell2mat(input.tDotDirectionDeg);
        nStim{i} = double(input.dotDirectionStepN);
    elseif strcmp(expt(iexp).dataset_visstim{i},'gratings') == 1
        DirectionDeg{i} = double(cell2mat(input.tGratingDirectionDeg));
        nStim{i} = double(input.gratingDirectionStepN);
    end
    nRep{i} = size(data_TC{i},1)./((nON+nOFF)*nStim{i});
    nTrials{i} = (nStim{i}.*nRep{i});
    Dirs{i} = unique(DirectionDeg{i});
    
    dF = zeros(size(data_TC{i}));
    tempDFoverF = zeros(size(data_TC{i}));
    for itrial = 1:nTrials{i}
        indAll = stimOFF_ind{i}(itrial):stimOFF_ind{i}(itrial)+(nOFF+nON-1);
        indF = stimOFF_ind{i}(itrial)+5:stimOFF_ind{i}(itrial)+(nOFF-1);
        dF(indAll,:) = bsxfun(@minus,data_TC{i}(indAll,:),mean(data_TC{i}(indF,:),1));
        tempDFoverF(indAll,:) = bsxfun(@rdivide,dF(indAll,:),mean(data_TC{i}(indF,:),1));
    end
    
    dFoverF{i} = tempDFoverF;
    clear input dF tempDFoverF
end
%% get tuning curves for 4 directions for each dataset

dFoverF_trialChop = {};
dFoverF_meanRespDir = {};

for i = 1:nrun
    dFoverF_trialChop{i} = zeros(10+nON,size(dFoverF{i},2),nTrials{i});
    for itrial = 1:nTrials{i}
        dFoverF_trialChop{i}(:,:,itrial) = dFoverF{i}(stimON_ind{i}(itrial)-10:stimON_ind{i}(itrial)+(nON-1),:);
    end
end

for i = 1:nrun
    dFoverF_meanRespDir{i} = zeros(size(dFoverF_trialChop{i},1),size(dFoverF_trialChop{i},2),nStim{i});
    for istim = 1:nStim{i}
        trials = find(DirectionDeg{i}(:,1:nTrials{i}) == Dirs{i}(istim));
        dFoverF_meanRespDir{i}(:,:,istim) = mean(dFoverF_trialChop{i}(:,:,trials),3);
    end
    dirRespCurves{i} = squeeze(mean(dFoverF_meanRespDir{i}(11:end,:,:),1));
    dirRespCurvesErr{i} = squeeze(std(dFoverF_meanRespDir{i}(11:end,:,:),[],1));
end

%% plot tuning curves overlayed for all experiments (all cells)

col_mat = {'k', 'g' , 'b' , 'r'};


% for a few randomly chosen cells
nCells = size(dirRespCurves{1},1);
cellSet = randperm(nCells,100);
figure;
for icell = 1:length(cellSet);
    subplot(3,3,icell)
    for i = 1:nrun
        errorbar(Dirs{i},dirRespCurves{i}(cellSet(icell),:),dirRespCurvesErr{i}(cellSet(icell),:),'Color',col_mat{i})
        hold on
    end
    if icell == 8;
        legend(expt(iexp).dataset_exp,'Location','SouthEast')
    end
    title(['cell ' num2str(cellSet(icell))])
    xlim([-10 280])
    ylim([-.2 1])
end

figure;
for icell = 1:length(cellSet)
    subplot(10,10,icell)
    plot(Dirs{3},dirRespCurves{3}(cellSet(icell),:))
    title(['cell ' num2str(cellSet(icell))])
    xlim([-10 280])
    ylim([-.2 .5])
end

cMap = colormap(jet(nStim{3}));
figure;
for icell = 1:length(cellSet)
    subplot(10,10,icell)
    for istim = 1:nStim{3}
        plot(dFoverF_meanRespDir{3}(:,icell,istim),'color',cMap(istim,:))
        hold on
        vline(10,'k')
        title(['cell ' num2str(cellSet(icell))])
        hold on
    end
end
%% find DSI (mov dots and gratings) and OSI


%% save figures


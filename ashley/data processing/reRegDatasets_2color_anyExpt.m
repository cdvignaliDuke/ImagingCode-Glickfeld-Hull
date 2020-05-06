clear 
clear global
close all
ds = 'retWithPupilMirror'; %dataset info
rc = behavConstsAV; %directories
eval(ds)
slct_expt = 1; %which expt from ds to analyze
doPreviousReg = false;
findNewRegImg = false;
doGreenOnly = true;
nRedFrames = nan;

%% enter in imaging run info
exptType = 'retinotopy';
if strcmp(exptType,'retinotopy')
    runNames = {'noMirror','withMirror'};
    runs = {expt(slct_expt).retNoMirror{1},expt(slct_expt).retWithMirror{1}};
    run_times = {expt(slct_expt).retNoMirror{2},expt(slct_expt).retWithMirror{2}};
    nrun = length(runs);

    maskDataInd = 1;
end
%%

% analysis path
mouse = expt(slct_expt).mouse;
sn = mouse;
expDate = expt(slct_expt).date;
if strcmp(rc.name, 'ashle')
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    if ~exist(fullfile(fn),'dir')
        mkdir(fn)
    end
    if ~exist(fullfile(fn,runFolder),'dir')
        mkdir(fn,'data processing')
    end
    fnout = fullfile(fn,'data processing');
else
    error('you do not belong here')
end

% load and register data
data_g_allRun = cell(1,nrun);
mw_allRun = cell(1,nrun);
for irun = 1:nrun
    runFolder = runs{irun};
    expTime = run_times{irun};
    fName = [runFolder '_000_000'];
    if doGreenOnly
        data_g_temp = loadsbx_choosepmt(...
            1,mouse,expDate,runFolder,fName);
        mw_allRun{irun} = loadMworksFile(sn,expDate,expTime);
        if irun == 1
            if findNewRegImg
                regImgStartFrame = compareRegImg_2color(rc,ds,sn,expDate,data_g_allRun{irun});
            else
                regImgStartFrame = expt(slct_expt).regImgStartFrame;
            end
            regImg = mean(data_g_temp(:,:,regImgStartFrame:(regImgStartFrame+99)),3);
        end
        
        data_g_temp = data_g_temp - min(data_g_temp(:));
        [outs,data_g_allRun{irun}] = stackRegister(data_g_temp,regImg);
        if ~exist(fullfile(fn,runFolder),'dir')
            mkdir(fn,runFolder)
        end
        save(fullfile(fn,runFolder,'regOuts&Img'),'outs','regImg')        
    else
        error('Need to add in code for channel 2')
    end    
end
clear data_g_temp

%% get mask
nBaselineFr = params.nBaselineMs./1000.*params.frameRate;
if doGreenOnly
    if strcmp(exptType,'retinotopy')
        mw = mw_allRun{maskDataInd};
        Az = celleqel2mat_padded(mw.tGratingAzimuthDeg);
        El = celleqel2mat_padded(mw.tGratingElevationDeg);
        Azs = unique(Az);
        Els = unique(El);
        if min(Els,[],2)<0
            Els = fliplr(Els);
        end
        nStim = length(Azs).*length(Els);
        on = mw.nScansOn;
        off = mw.nScansOff;

        nTrials = length(Az);
        nFr = size(data_g_allRun{maskDataInd},3);
        if nFr > ((on+off)*nTrials)
            nFr = ((on+off)*nTrials);
            d = data_g_allRun{maskDataInd}(:,:,1:nFr);
        else
            d = data_g_allRun{maskDataInd};
        end
        d = double(reshape(d,[size(d,1),size(d,2),on+off,nTrials]));
        f = mean(d(:,:,(off-nBaselineFr):off,:),3);
        dff = (d-f)./f;
        clear d f

        stims = nan(2,nStim);
        maskImg = nan(size(dff,1),size(dff,2),nStim);
        start = 1;
        for iEl = 1:length(Els)
            ind1 = find(El == Els(iEl));
            for iAz = 1:length(Azs)
                stims(:,start) = [Els(iEl), Azs(iAz)];
                ind2 = find(Az == Azs(iAz));
                ind = intersect(ind1,ind2);
                maskImg(:,:,start) = mean(mean(dff(:,:,(off+2):(off+on),ind),4),3);
                start = start+1;
            end
        end    
        clear mw dff
    end
    figure;colormap gray;imagesc(mean(maskImg,3));
    writetiff(maskImg,fullfile(fnout,'maskImages'))
% save some F images as tifs
    nImg = 100;
    nFr = cellfun(@(x) size(x,3),data_g_allRun);
    for irun = 1:nrun
        frInterval = 1:floor((nFr(irun)-100)/nImg):nFr(irun)-100;
        F = nan(size(data_g_allRun{irun},1),size(data_g_allRun{irun},2),nImg);
        for iimg = 1:nImg
            F(:,:,iimg) = mean(data_g_allRun{irun}(:,:,frInterval(iimg):frInterval(iimg)+99),3);
        end
        writetiff(F,fullfile(fnout,['FImages_' runNames{irun}]))
    end
end

%% select cells
if doGreenOnly
    mask_cell = maskFromMultiMaxDFFStack(maskImg);
    close all
    figure; setFigParams4Print('portrait')
    imagesc(mask_cell);
    title({sprintf('%s cells',...
        num2str(length(unique(mask_cell(:)))-1));[mouse '-' expDate]})
        print(fullfile(fnout,'final_mask_cells'),'-dpdf')
        save(fullfile(fnout,'final_mask_cells.mat'),'mask_cell');
end
%% extract and save tcs
buf = 4;
np = 6;
if doGreenOnly    
    data_g_tc_allRun = cellfun(@(x) stackGetTimeCourses(x,mask_cell),...
        data_g_allRun,'unif',0);
    data_g_tc_subnp_allRun = cellfun(@(x,y) ...
        getWeightedNeuropilTimeCourse(x,y,mask_cell,buf,np),...
        data_g_allRun,data_g_tc_allRun,'unif',0);
    for irun = 1:nrun
        d = data_g_tc_allRun{irun};
        d_np = data_g_tc_subnp_allRun{irun};
        save(fullfile(fn,runs{irun},'timecourses_cells.mat'),...
            'd_np','d','buf','np')
    end
end

clear all
close all
ds = 'retWithPupilMirror';
eval(ds)
%%
rc = behavConstsAV;
fnout = fullfile(rc.ashleyAnalysis,'Expt Summaries',ds);

%% analysis params
nbaselinefr = round((params.nBaselineMs/1000)*params.frameRate);
% tt_stimOnS_label = -0.75:0.25:0.75;
% tcStartTimeS = -1;

nexp = size(expt,2);
runNames = {'noMirror','withMirror'};
%% get RF fits for each experiment
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    runs = {expt(iexp).retNoMirror{1},expt(iexp).retWithMirror{1}};
    run_times = {expt(iexp).retNoMirror{2},expt(iexp).retWithMirror{2}};
    nrun = length(runs);
    fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    
    load(fullfile(fn,runs{1},'timecourses_cells.mat'));
    mw = loadMworksFile(mouse,expDate,run_times{1});
    tc = d_np;
    clear d d_np
    [nFr,nCells] = size(tc);
    
    on = mw.nScansOn;
    off = mw.nScansOff;
    tAz = celleqel2mat_padded(mw.tGratingAzimuthDeg);
    tEl = celleqel2mat_padded(mw.tGratingElevationDeg);
    azs = unique(tAz);
    els = unique(tEl);
    nAz = length(azs);
    nEl = length(els);
    if min(els)<0
        els = fliplr(els);
    end
    nStim = nAz*nEl;
    nTrials = length(tAz);
    if nFr > nTrials*(on+off)
        nFr = nTrials*(on+off);
        tc = tc(1:nFr,:);
    end
    
    tc_tr = reshape(tc,[on+off,nTrials,nCells]);    
    f = mean(tc_tr((off-nbaselinefr):off,:,:),1);
    dff = (tc_tr - f)./f;
    dff_resp = squeeze(mean(dff((off+1):(on+off),:,:),1));
    dff_base = squeeze(mean(dff(round(off/2):off,:,:),1));
    
    stims = nan(2,nStim);
    stimsTrialInd = cell(1,nStim);
    dff_stim = nan(on+off,nStim,nCells);
    istim = 1;
    for iel = 1:nEl
        for iaz = 1:nAz
            ind = tEl == els(iel) & tAz == azs(iaz);
            stims(:,istim) = [els(iel) azs(iaz)];
            dff_stim(:,istim,:) = mean(dff(:,ind,:),2);
            stimsTrialInd{istim} = find(ind);
            istim = istim+1;
        end
    end
    
    tuning = squeeze(mean(dff_stim((off+1):(on+off),:,:),1));
    retinotopy = reshape(tuning,[nEl,nAz,nCells]);
    
    figure; colormap gray
    for icell = 1:25
        subplot(5,5,icell)
        imagesc(retinotopy(:,:,icell))
        figAxForm
        title(sprintf('Cell %s',num2str(icell)))
    end
    suptitle([mouse '-' expDate])
    
    [AzAz, ElEl] = meshgrid(azs,els);
    grid2.AzAz = AzAz;
    grid2.ElEl = ElEl;    
    dAz = median(diff(azs));
    dEl = median(diff(els));
    Az_vec00 = azs(1):(dAz/10):azs(end);
    El_vec00 = els(1):(dEl/10):els(end);
    [AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
    grid2.AzAz00 = AzAz00;
    grid2.ElEl00 = ElEl00;
%     Nshuf = 50;

    for icell = 1:nCells
        b = retinotopy(:,:,icell);
        data = b';
        PLOTIT_FIT = 1;
        SAVEALLDATA = 0;
        Fit_2Dellipse_LG_Ret_AW
    end
end
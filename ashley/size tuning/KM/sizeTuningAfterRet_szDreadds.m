%% size tuning script by kevin
% modifying newScript.m, singleChannelTCScript.m
% modifying sizeTuningAfterRet.m
% reads in Size Tuning Curves from single channel imaging data

% here we extract size tuning, using the neuron mask from retOnly run

%% get path names
clear all;clc;

ds = 'szTuning_dreadds_V1';
iexp = 1;
rc = behavConstsAV;
eval(ds)

%%
mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
runs = expt(iexp).sizeTuningFolder;
expTime = expt(iexp).sizeTuningTime;
nrun = length(runs);
frame_rate = params.frameRate;

retFolder = cell2mat(expt(iexp).retinotopyFolder);
fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

fprintf(['2P imaging size tuning analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' expDate '\nExperiments:\n'])
for irun=1:nrun
    fprintf([runs{irun} ' - ' expTime{irun,:} '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    %load imaging data
    dataFolder = runs{irun};
    fprintf(['\nLoading run ' num2str(irun) '...'])
    fName = [dataFolder '_000_000'];
    switch expt(iexp).saveLoc
        case 'ashley'
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName);
        case 'kevin'
            fdir = ['\\crash.dhe.duke.edu\data\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,[],fdir);
        otherwise
            error('identify data source')
    end

    % load behavior data
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' mouse '-' expDate '-' expTime{irun} '.mat'];
    load(fName);
    ntrials = size(input.tGratingDirectionDeg,2);
    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    nframes = size(data_temp,3);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    temp(irun).run = mat2cell(irun.*ones(1,ntrials),1,ntrials);
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete')
end
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
for i=1:length(input.tGratingContrast) % replace int64(con=1) with double
    if ~(class(input.tGratingContrast{i})=="double")
        input.tGratingContrast{i} = double(input.tGratingContrast{i});
    end
end
clear data_temp
clear temp

%% Check stability
regIntv = 3000;
nep = floor(size(data,3)./regIntv);
fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

% plot 500 frame means at each register interval
[n, n2] = subplotn(nep);
figure(1);clf;
colormap(gray)
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end
%% Register data

if nrun>1 
    dataFolder = [runs{1} '_' runs{nrun}];
end
        
fprintf('\nBegin registering...\n')   

% load from ref folder:
% reg_shifts.mat (out, data_avg)
% mask_cell.mat ()
% trialData.mat ()
% then create new directory and save analysis
fprintf(['Reference specified: ' retFolder '\n'])

% load from folder specified by ref_str
fprintf(['Loading from folder: ' retFolder '\n'])
load(fullfile(fnout,retFolder,[mouse '_' expDate 'ret_reg_shifts.mat']))
load(fullfile(fnout,retFolder,[mouse '_' expDate '_mask_cell.mat']))

% register
fprintf('stackRegister with reference image\n')
[out, data_reg] = stackRegister(data,data_avg);

% save
fprintf('Registration complete, now saving...\n')
mkdir(fullfile(fnout,dataFolder))
save(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']), 'out', 'data_avg','meanrng')
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']), 'input')

clear data % depending on memory

%% test stability
% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

% figure 3 shows imagesq (scaled?) of the registered data 1-10000, then prints to FOV_avg.mat
fprintf('Examine FOV\n')
figure(3);
if nframes>10000
    imagesq(mean(data_reg(:,:,1:10000),3));
else
    imagesq(mean(data_reg,3));
end
truesize;
% print to file
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOV_avg.pdf']))


%% find activated cells

% check useFilt (only for dFoF calculation and segmentation, timecourse extraction not filtered)
useFilt = 0;

% calculate dF/F
fprintf('\nBegin image analysis...\n')

% max by trial type
if isfield(input, 'nScansOn')
    % nScansOn -> passive vis stim ret
    fprintf('nScansOn method - get dF/F\n')
    
    % load defining variables
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    sz = size(data_reg);
    
    fprintf('Calculating dF/F...\n')
    
    % split into trials by reshape
    data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_fOn = mean(data_tr(:,:,nOff:(nOn+nOff),:),3);
    data_dfof = squeeze((data_fOn-data_f)./data_f);
    
    % previous code would do data_dfof for every frame
    % data_dfof = (double(data_tr)-data_f)./data_f;
    
    clear data_tr data_f data_fOn % data_fLate
    
    if useFilt
        fprintf('Filtering images and time averaging...\n')
        switch useFilt
            case 1
                % filter with a 20x20 gaussian sigma=0.7
                fprintf('useFilt=1: Gaussian smoothing filter\n')
                myfilter = fspecial('gaussian',[20 20], 0.7);
                data_dfof = imfilter(data_dfof,myfilter);
            case 2
                % filter with median filter (default 3x3)
                fprintf('useFilt=2: Median filter\n')
                data_dfof = medfilt2(data_dfof, [3 3]);
            case 3
                % filter with wiener adaptive filter (default 3x3)
                fprintf('useFilt=3: Wiener adaptive filter\n')
                data_dfof = wiener2(data_dfof,[5 5]);
        end
    end
    
    fprintf('done\n')
end

% with dF/F, average by each stimulus, depending on experiment
if input.doRetStim
    % doRetStim -> retinotopy
    % requires data_dfof from above nScansOn method
    fprintf('input.doRetStim method - varying Az+El position\n')
    fprintf('ERROR: this script only does size tuning not retinotopy')
    
elseif input.doSizeStim
    % doSizeStim -> size tuning
    % requires data_dfof from above nScansOn method
    fprintf('input.doSizeStim method - varying grating diameter\n')
    
    % store vectors for grating diameter
    % what is celleqel2mat_padded?
    sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
    szs = unique(sz_mat);
    
    nStim = length(szs);
    fprintf([num2str(nStim) ' unique size stimuli\n'])
    
    fprintf('Averaging dF/F for each size...\n')
    data_dfof_avg = zeros(sz(1), sz(2), nStim);
    for isz = 1:nStim
        ind = find(sz_mat == szs(isz));
        data_dfof_avg(:,:,isz) = mean(data_dfof(:,:,ind),3);
    end
    
    % plot dF/F for all stimuli
    figure(5);clf;
    [n, n2] = subplotn(nStim);
    for i = 1:nStim
        subplot(n,n2,i);
        imagesc(data_dfof_avg(:,:,i));
        clim([0 max(data_dfof_avg(:))])
        title(['Size: ' num2str(szs(i)) ' deg'])
    end
    % print to pdf
    print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOV_resp_Size.pdf']), '-dpdf')
    
elseif input.doTFStim && ~input.doMatrix
    % doTFStim + !doMatrix -> temporal frequency?
    % requires data_dfof from above nScansOn method
    fprintf('input.doTFStim method - varying grating TF+SF\n')
    fprintf('ERROR: this script only does size tuning not TF tuning')
    
elseif input.doDirStim
    % doDirStim -> directional stimulation tuning curve
    % requires data_dfof from above nScansOn method
    fprintf('input.doDirStim method - varying grating direction\n')
    fprintf('ERROR: this script only does size tuning not direction tuning')
    
end

% take max across stimuli
fprintf('Final step: take max across stimuli\n')
data_dfof_max = max(data_dfof_avg,[],3);
figure(9);clf;
if input.doDirStim
    imagesc(max(data_dfof_avg_ori,[],3))
    clim([0 max(data_dfof_avg_ori(:))])
else
    imagesc(data_dfof_max)
    clim([0 max(data_dfof_max(:))])
end
title('Maximum dF/F across all stimuli')

% save stimActFOV.mat containing: data_dfof_max, data_dfof_avg, nStim
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg', 'nStim')


%% load cell masks from retinotopy runs

fprintf(['Loading masks from retinotopy runs: ' retFolder '\n'])

% loads 'mask_cell', 'mask_np'
load(fullfile(fnout, retFolder, [mouse '_' expDate '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

% translate if necessary (should not be necessary after register with ref)
% mask_cell = imtranslate(mask_cell, [0 2]); % [x(+right) y(+down)]
% mask_np = imtranslate(mask_np, [0 2]);

% load ret fit data, in order to select only goodfit_ind cells
% fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
% load(fn_out);

% [mask_cell, nCells] = bwlabel(mask_cell); % bwlabel labels all individual cells
nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

figure(11);clf;
[n, n2] = subplotn(nStim);
for i = 1:nStim
    subplot(n,n2,i);
    shade_img = imShade(data_dfof_avg(:,:,i), mask_cell);
    imagesc(shade_img)
    if input.doSizeStim
        title([num2str(szs(i)) ' deg'])
    elseif input.doRetStim
        title(num2str(Stims(i,:)))
    end
    clim([0 max(data_dfof_avg(:))])
    colormap(gray)
end
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_retSz_maskOverlap.pdf']),'-dpdf')


%% Get time courses, including neuropil subtraction

fprintf('\nBegin time course extraction...\n')
down = 5; % downsample rate
data_reg_down  = stackGroupProject(data_reg,down);
fprintf(['Downsampling at M=' num2str(down) '\n'])

fprintf('Extracting cell signal...\n')
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
fprintf([num2str(nCells) ' total cells extracted\n'])

fprintf('Extracting neuropil signal...\n')
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
    np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
    fprintf(['Cell #' num2str(i) ' / ' num2str(nCells) '\n'])
end

fprintf('Subtract neuropil signal, maximizing skewness\n')
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
fprintf(['Maximum skews: ' num2str(max_skew) '\n'])
fprintf(['at inds: ' num2str(ind) '\n'])

np_w = 0.01*ind;
fprintf(['np_w = ' num2str(np_w) '\n'])

npSub_tc = data_tc-(tcRemoveDC(np_tc).*np_w);
clear data_reg_down

fprintf('Neuropil subtraction complete, saving data...\n')

save(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']), 'input')

fprintf('Calculating dF/F...\n')
%get dF/F
nCells = size(npSub_tc,2);
tc_mat = zeros(nOn+nOff, nCells, ntrials);
% reshape into trials, could do this in one line with reshape?
for itrial = 1:ntrials
    tc_mat(:,:,itrial) = npSub_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)),:);
end
tc_f = mean(tc_mat(nOff/2:nOff,:,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
clear tc_mat tc_f

fprintf('Time course extraction complete.\n')


%% plot tuning
load(fullfile(fnout, retFolder, [mouse '_' expDate '_lbub_fits.mat']));

tuning_mat = zeros(nStim, nrun, 2, nCells);
tRun = [];
for irun = 1:nrun
    tRun = [tRun irun.*ones(1,input.trialsSinceReset(1,irun))];
end

Ind_struct = [];
if nCells<36
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt = (1-nOff:nOn)*(1000./frame_rate);
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if ~sum(iCell==goodfit_ind) % only plot goodfit cells
        continue
    end
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind_temp = find(sz_mat == szs(iStim));
        plot(tt', squeeze(mean(tc_dfof(:,iCell,ind_temp),3)))
        title(iCell)
        for irun = 1:nrun
            ind = intersect(find(tRun == irun),ind_temp);
            hold on
            tuning_mat(iStim,irun,1,iCell) = mean(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
            tuning_mat(iStim,irun,2,iCell) = std(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
            Ind_struct(iStim).all_trials = ind;
        end
    end
    ylim([-0.05 0.25])
    vline(nOff)
    start = start + 1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if ~sum(iCell==goodfit_ind) % only plot goodfit cells
        continue
    end
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for irun = 1:nrun
        ret_mat = tuning_mat(:,irun,1,iCell);
        ret_mat = ret_mat';
        plot(szs,ret_mat)
        hold on 
    end
    title(num2str(chop(max(mean(tuning_mat(:,:,1,iCell),2),[],1),2)))
    start = start+1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning.mat']), 'tc_dfof', 'tuning_mat', 'szs', 'Ind_struct')

%% present selected cells
% 
% % input stimulus location based on experimental choice
% stimEl = double(input.gratingElevationDeg);
% stimAz = double(input.gratingAzimuthDeg);
% fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])
% 
% % load ret fits - loads 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind'
% fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
% load(fn_out);
% cellAz = lbub_fits(:,4,4);
% cellEl = lbub_fits(:,5,4);
% fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
% 
% fprintf('Calculating cell RF distances to stimulus...\n')
% cellDists = sqrt((cellAz-stimAz).^2+(cellEl-stimEl).^2);
% 
% % compare distance to select cells
% cutOffRadius = 10;
% fprintf(['Isolating cells with RF centers within ' num2str(cutOffRadius) ' degrees\n'])
% 
% cutCells = find(cellDists < cutOffRadius)';
% goodCutCells = intersect(cutCells,goodfit_ind);
% nCutCells = length(goodCutCells);
% fprintf([num2str(nCutCells) ' cells selected\n'])
% 
% nSize = length(szs);
% conTrials = cell2mat(input.tGratingContrast);
% cons = unique(conTrials);
% nCon = length(cons);
% conInds = cell(nCon,1);
% for i = 1:nCon
%     conInds{i} = find(conTrials == cons(i));
% end
% 
% sizeTune = zeros(nSize,nCon,nCutCells);
% sizeSEM = sizeTune;
% [n, n2] = subplotn(min([36 nCutCells]));
% 
% figure;
% start = 1;
% f = 1;
% suptitle(['Size Tuning Curves for cells within ' num2str(cutOffRadius) ' deg of stim'])
% for i = 1:nCutCells
%     iCell = goodCutCells(i);
%     
%     if start >36
%         set(gcf, 'Position', [0 0 800 1000]);
%         print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning_in' num2str(cutOffRadius) 'deg' num2str(f) '.pdf']), '-dpdf')
%         start = 1;
%         f= f+1;
%         figure;
%         suptitle(['Size Tuning Curves for cells within ' num2str(cutOffRadius) ' deg of stim'])
%     end
%     % this first section does all trials overlaid for each stim
% %     figure;
% %     for iCond = 1:nStim
% %         subplot(4,2,iCond)
% %         ind_all = Ind_struct(iCond).all_trials;
% %         plot(squeeze(tc_dfof(:,iCell,ind_all)))
% %         title(['Size: ' num2str(szs(iCond)) ' deg'])
% %         ylim([-0.05 0.4])
% %     end
% %     suptitle(['Cell #: ' num2str(iCell)])
%     % this second section shows all trials averaged
% %     figure;
% %     for iCond = 1:nStim
% %         subplot(4,2,iCond)
% %         ind_all = Ind_struct(iCond).all_trials;
% %         plot(squeeze(mean(tc_dfof(:,iCell,ind_all),3)))
% %         title(['Size: ' num2str(szs(iCond)) ' deg'])
% %         ylim([-0.05 0.4])
% %     end
% %     suptitle(['Cell #: ' num2str(iCell)])
% 
%     % this third section shows size tuning curves
%     
%     subplot(n,n2,start)
%     for iCon = 1:nCon
%         for iSize = 1:nSize
%             ind = intersect(Ind_struct(iSize).all_trials,conInds{iCon});
%             %stimOff = mean(mean(tc_dfof((nOff/2):nOff,iCell,ind_all),3),1);
%             stimOn = mean(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind),3),1);
%             sd = std(mean(tc_dfof((nOff+1):(nOff+nOn),iCell,ind),1));
%             sizeTune(iSize,iCon,i) = stimOn;
%             sizeSEM(iSize,iCon,i) = sd/sqrt(length(ind));
%         end
%         errorbar(szs,sizeTune(:,iCon,i),sizeSEM(:,iCon,i))
%         hold on
%     end
%     ylim([-1.2*(max(max(sizeSEM(:,:,i)))) 1.2*(max(max(sizeTune(:,:,i)))+max(max(sizeSEM(:,:,i))))])
%     title(['Cell #: ' num2str(iCell)])
%     hold off
%     
%     start=start+1;
% end
% % Construct a Legend with the data from the sub-plots
% hL = legend(num2str(cons'));
% % Programatically move the Legend
% newPosition = [0.5 0.07 0.2 0.2];
% newUnits = 'normalized';
% set(hL, 'Position', newPosition, 'Units', newUnits);
% 
% 
% set(gcf, 'Position', [0 0 800 1000]);
% print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning_in' num2str(cutOffRadius) 'deg' num2str(f) '.pdf']), '-dpdf')
% 
% fprintf('\ndone\n')
% return
% 
% 

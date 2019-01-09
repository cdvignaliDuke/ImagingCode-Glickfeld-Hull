%% Size Detection Behavior script by kevin
% modified from retOnly.m

% reads in Tuning Curves from single channel imaging data

%% get path names
clear all;clc;

mouse = 'i560';
date = '180212';
ImgFolder = char('001');
time = char('1356');
doFromRef = 1;
ref = char('002');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

fprintf(['2P imaging size detection task analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder(irun,:) ' - ' time(irun,:) '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    % load 2p imaging data
    CD = ['H:\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['H:\home\valerie\Data\2p\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['H:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    % load behavior/experimental data
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    % read in frames with sbxread
    nframes = min(input.counterValues{end}(end),info.config.frames);
    fprintf(['\nReading run ' num2str(irun) ' - ' num2str(nframes) ' frames \n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    toc
    
    temp(irun) = input;
    
    % store value of number of trials
    ntrials = size(temp(irun).tGratingContrast,2);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% First examine behavior

% define binary vectors for success, fail, ignore
tSuccess = strcmp(input.trialOutcomeCell,'success');
tFail = strcmp(input.trialOutcomeCell,'incorrect');
tIgnore = strcmp(input.trialOutcomeCell,'ignore');

% look at ignores to choose trial cutoff range
figure(1);clf; plot(tIgnore); title('Ignore trials')

% store vectors for trial side, size, contrast
lTrial = cell2mat(input.tLeftTrial);
rTrial = 1-lTrial;
siz_mat = cell2mat(input.tGratingDiameterDeg);
con_mat = cell2mat(input.tGratingContrast);
sizs = unique(siz_mat);
cons = unique(con_mat);
nSiz = length(sizs);
nCon = length(cons);

% for each stimulus pair, calculate percent success, define Ind_struct
nStim = nSiz*nCon;
fprintf([num2str(nStim) ' unique size+contrast stimuli\n'])

pctSuccBoth = zeros(1,nStim);
pctSuccLeft = zeros(1,nStim);
pctSuccRight = zeros(1,nStim);
pctFailBoth = zeros(1,nStim);
pctFailLeft = zeros(1,nStim);
pctFailRight = zeros(1,nStim);

start=1;
Stims = zeros(nStim,2);
Ind_struct = [];

indR = find(rTrial);
for iSiz = 1:nSiz
    indS = find(siz_mat == sizs(iSiz)); % find indices of each size
    for iCon = 1:nCon
        indC = find(con_mat == cons(iCon)); % find indices of each con
        ind = intersect(indS,indC); % choose common siz and con indices
        pctSuccBoth(start) = sum(tSuccess(ind))/sum(tSuccess(ind)+tFail(ind)+tIgnore(ind));
        pctFailBoth(start) = sum(tFail(ind))/sum(tSuccess(ind)+tFail(ind)+tIgnore(ind));
        Ind_struct(start).all_trials = ind;
        
        indL = setdiff(ind,indR); % only left trials
        pctSuccLeft(start) = sum(tSuccess(indL))/sum(tSuccess(indL)+tFail(indL)+tIgnore(indL));
        pctFailLeft(start) = sum(tFail(indL))/sum(tSuccess(indL)+tFail(indL)+tIgnore(indL));
        Ind_struct(start).left_trials = indL;
        
        ind = intersect(ind,indR); % only right trials
        pctSuccRight(start) = sum(tSuccess(ind))/sum(tSuccess(ind)+tFail(ind)+tIgnore(ind));
        pctFailRight(start) = sum(tFail(ind))/sum(tSuccess(ind)+tFail(ind)+tIgnore(ind));
        Ind_struct(start).right_trials = ind;
        
        % stores combination to Stims and iterates start
        Stims(start,:) = [sizs(iSiz) cons(iCon)];
        start = start+1;
    end
end

% trial numbers per stim
tNumBoth = arrayfun(@(x) numel(x.all_trials), Ind_struct);
tNumLeft = arrayfun(@(x) numel(x.left_trials), Ind_struct);
tNumRight = arrayfun(@(x) numel(x.right_trials), Ind_struct);

% size x con map, with trial numbers as both(right)
% currently NOT excluding ignore trials
trialMapBoth = reshape(tNumBoth,nCon,nSiz);
trialMapLeft = reshape(tNumLeft,nCon,nSiz);
trialMapRight = reshape(tNumRight,nCon,nSiz);
figure(2);clf
% this for loop contains the trial map plot, just wanted to fold it up
for foldplot=1
    subplot(2,2,[1 2])
    imagesc(trialMapBoth)
    title('Trial counts: both(right)')
    xlabel('Size (deg)')
    ylabel('Contrast')
    xticks(1:nSiz)
    yticks(1:nCon)
    set(gca,'XTickLabel',sizs);
    set(gca,'YTickLabel',cons);
    for i=1:nCon
        for j=1:nSiz
            text(j,i,[num2str(trialMapBoth(i,j)) '(' num2str(trialMapRight(i,j)) ')'],'HorizontalAlignment','center')
        end
    end
    subplot(2,2,3)
    imagesc(trialMapLeft)
    title('Trial counts: left')
    xlabel('Size (deg)')
    ylabel('Contrast')
    xticks(1:nSiz)
    yticks(1:nCon)
    set(gca,'XTickLabel',sizs);
    set(gca,'YTickLabel',cons);
    for i=1:nCon
        for j=1:nSiz
            text(j,i,num2str(trialMapLeft(i,j)),'HorizontalAlignment','center')
        end
    end
    subplot(2,2,4)
    imagesc(trialMapRight)
    title('Trial counts: right')
    xlabel('Size (deg)')
    ylabel('Contrast')
    xticks(1:nSiz)
    yticks(1:nCon)
    set(gca,'XTickLabel',sizs);
    set(gca,'YTickLabel',cons);
    for i=1:nCon
        for j=1:nSiz
            text(j,i,num2str(trialMapRight(i,j)),'HorizontalAlignment','center')
        end
    end
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Behav_trialCounts.pdf']), '-dpdf')

% plot pct success/fail/ignore with varying contrast, size
figure(3);clf;
subplot(2,2,[1 2])
bar([pctSuccBoth', pctFailBoth', 1-pctSuccBoth'-pctFailBoth'], 'stacked')
title('Both sides behavioral outcomes')
ylabel('Percent trials')
legend('Success','Fail','Ignore','Location','best')
subplot(2,2,3)
bar([pctSuccLeft', pctFailLeft', 1-pctSuccLeft'-pctFailLeft'], 'stacked')
title('Left side behavioral outcomes')
ylabel('Percent trials')
legend('Success','Fail','Ignore','Location','best')
subplot(2,2,4)
bar([pctSuccRight', pctFailRight', 1-pctSuccRight'-pctFailRight'], 'stacked')
title('Right side behavioral outcomes')
ylabel('Percent trials')
legend('Success','Fail','Ignore','Location','best')

% plot psychometric (just percent success normalized with only failures)
pctSuccBothNorm = pctSuccBoth./(pctSuccBoth+pctFailBoth);
pctSuccLeftNorm = pctSuccLeft./(pctSuccLeft+pctFailLeft);
pctSuccRightNorm = pctSuccRight./(pctSuccRight+pctFailRight);
figure(4);clf;
Markers = {'o-','+-','*-','x-','s-','d-','^','v','h'};
subplot(2,1,1)
for iSiz = 1:nSiz
    indS = find(Stims(:,1)==sizs(iSiz));
    indS = intersect(indS,1:160); % cutoff at trial 160
    x = Stims(indS,2); % use contrast as x
    hold on
    plot(x,pctSuccBothNorm(indS),Markers{iSiz})
    hold off
end
title('Both sides combined')
xlabel('Contrast')
ylabel('%correct')
ylim([0 1])
legend(num2str(sizs'),'Location','best')
subplot(2,1,2)
for iSiz = 1:nSiz
    indS = find(Stims(:,1)==sizs(iSiz));
    indS = intersect(indS,1:160); % cutoff at trial 160
    x = Stims(indS,2); % use contrast as x
    hold on
    plot([-fliplr(x') x'],[fliplr(pctSuccLeftNorm(indS)) pctSuccRightNorm(indS)],Markers{iSiz})
    hold off
end
title('Left + Right')
xlabel('Contrast')
ylabel('%correct')
ylim([0 1])
legend(num2str(sizs'),'Location','best')
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Behav_psychometric.pdf']), '-dpdf')

%% Choose register interval
regIntv = 15000;
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

chooseInt = 1; %nep/2

sz = uint64(size(data));

fprintf('\nBegin registering...\n')
if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']),'file')
    % checks if analysis already present
    % load reg_shifts.mat (out, data_avg) and save the current input file
    fprintf('Found previous analysis! Loading...\n')
    
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    
    % register
    fprintf('stackRegister_MA, using shifts from previous registration\n')
    % uses previous registration shifts (out) to re-register data quickly
    [outs, data_reg]=stackRegister_MA(data,[],[],double(out));
    fprintf('Previous registration loaded...\n')
    
    % save new input
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    
elseif doFromRef
    % if doFromRef specified, along with ref (ref needs to exist, no error catch)
    % load from ref folder:
    % reg_shifts.mat (out, data_avg)
    % mask_cell.mat ()
    % trialData.mat ()
    % then create new directory and save analysis
    fprintf(['Reference specified: ' ref '\n'])
    
    ref_str = ['runs-' ref];
    
    % use multiple refs?
    %     ref_str = ['runs-' ref(1,:)];
    %     if size(ref,1)>1
    %         ref_str = [ref_str '-' ref(end,:)];
    %     end
    
    % load from folder specified by ref_str
    fprintf(['Loading from folder: ' ref_str '\n'])
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    
    % register
    fprintf('stackRegister with reference\n')
    [out, data_reg] = stackRegister(data,data_avg);
    
    % save
    fprintf('Registration complete, now saving...\n')
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    % else means no previous analysis present
    % use data_avg selected above (could move here?)
    % then create new directory and save analysis
    fprintf('Creating new analysis!\n')
    
    meanrng = regIntv*(chooseInt)+(1:500);
    data_avg = mean(data(:,:,meanrng),3);
    fprintf(['\nRegister frame averaged from ' num2str(meanrng(1)) ' - ' num2str(meanrng(end)) '\n'])
    
    % register
    fprintf('stackRegister\n')
    [out, data_reg] = stackRegister(data,data_avg);
    
    % save
    fprintf('Registration complete, now saving...\n')
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.mat']))

%% find activated cells
% calculate dF/F
fprintf('\nBegin image analysis...\n')

useFilt = 0;
max_trial = 300; % choose by examining tIgnore plot (behavior fig 1)

% load defining variables
nITI = input.nFramesIti;
nDelay = input.nFramesDelay;
nOn = input.nFramesOn;
nOff = input.nFramesTooFast;
ntrials = size(input.tGratingContrast,2);

fprintf(['Truncating to ' num2str(max_trial) '/' num2str(ntrials) ' trials\n'])
% truncate lTrial, rTrial, siz_mat, con_mat to max_trial
lTrial = lTrial(1:max_trial);
rTrial = rTrial(1:max_trial);
siz_mat = siz_mat(1:max_trial);
con_mat = con_mat(1:max_trial);

data_tr = zeros(sz(1),sz(2),nITI+nDelay+nOn+nOff, max_trial);
% split into trials sequentially
fprintf('Splitting into trials...\n')
for iTr = 1:max_trial
    cOn = input.cStimOn{iTr};
    data_tr(:,:,:,iTr) = data_reg(:,:,cOn-nITI-nDelay+1:cOn+nOn+nOff);
end

fprintf('Calculating dF/F...\n')
data_f = mean(data_tr(:,:,nITI-30:nITI,:),3);
data_fOn = mean(data_tr(:,:,nITI+nDelay+(1:nOn),:),3);
data_dfof = squeeze((data_fOn-data_f)./data_f);
clear data_tr data_f data_fOn

if useFilt
    % filter with a 20x20 gaussian sigma=0.7
    % then average over time from nOff:(nOff+nOn), all time when stim off?
    % squeeze to 3d matrix (2Dimage x nStim)
    fprintf('Filtering images and time averaging...\n')
    myfilter = fspecial('gaussian',[20 20], 0.7);
    %data_dfof_unfilt = data_dfof;
    data_dfof = imfilter(data_dfof,myfilter);
end

fprintf('done\n')

% with dF/F, average by each trial type (side+size+contrast)
% requires data_dfof from above
fprintf('Average for each condition (side+size+contrast)\n')
fprintf([num2str(nStim) ' unique size+contrast stimuli\n'])

fprintf('Averaging dF/F for each size+con condition on right trials...\n')
data_dfof_avg = zeros(sz(1), sz(2), nStim);

for iStim = 1:nStim
    ind = Ind_struct(iStim).right_trials; %right trial inds for this stim
    ind = intersect(ind,1:max_trial); %truncate to max_trial
    ind = setdiff(ind,find(tIgnore)); %remove ignore trials
    data_dfof_avg(:,:,iStim) = mean(data_dfof(:,:,ind),3);
end

% plot dF/F for all stimuli
figure(4);clf;
for i = 1:nStim
    subplot(nSiz,nCon,i);
    imagesc(data_dfof_avg(:,:,i));
    clim([0 max(data_dfof_avg(:))])
    title(['Size: ' num2str(Stims(i,1)) ', Con: ' num2str(Stims(i,2))])
end
% print to pdf
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp_Ret.pdf']), '-dpdf')

% take max across stimuli
fprintf('Final step: take max across stimuli\n')
data_dfof_max = max(data_dfof_avg,[],3);
figure(5);clf;
imagesc(data_dfof_max)
clim([0 max(data_dfof_max(:))])
title('Maximum dF/F across all stimuli')

% save stimActFOV.mat containing: data_dfof_max, data_dfof_avg, nStim
set(gcf, 'Position', [0 0 800 1000]);
save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg', 'nStim')

% examine total image average over all trials
fprintf('Examine image average dF/F by trial\n')
tot_avg = squeeze(mean(mean(data_dfof,1),2));
figure(6);clf;
plot(1:max_trial,tot_avg)
title('Total image average dF/F vs trial')
xlabel('Trial')
ylabel('dF/F')

%% cell segmentation
% here use GUI to select cells

fprintf('\nBegin cell segmentation...')
mask_all = zeros(sz(1), sz(2));
mask_exp = mask_all;
mask_data = data_dfof_avg;

% start with max projection
fprintf('\nFirst Step: Max Projection\n')
mask_data_temp = data_dfof_max;
mask_data_temp(mask_exp >= 1) = 0;
bwout = imCellEditInteractive(mask_data_temp);
mask_all = mask_all+bwout;
mask_exp = imCellBuffer(mask_all,3)+mask_all;

% by each stim
for iStim = 1:size(mask_data,3)
    fprintf(['Stim ' num2str(iStim) ' / ' num2str(size(mask_data,3)) '\n'])
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(mask_exp >= 1) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell = bwlabel(mask_all); % bwlabel labels all individual cells

%imagesc(data_dfof_max)

% repeat for max projection
fprintf('\nFinal Step: Max Projection\n')
mask_data_temp = data_dfof_max;
mask_data_temp(mask_exp >= 1) = 0;
bwout = imCellEditInteractive(mask_data_temp);
mask_all = mask_all+bwout;

[mask_cell, nCells] = bwlabel(mask_all); % bwlabel labels all individual cells
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

figure(11);clf;
[n, n2] = subplotn(nStim);
for i = 1:nStim
    subplot(n,n2,i);
    shade_img = imShade(data_dfof_avg(:,:,i), mask_all);
    imagesc(shade_img)
    title(num2str(Stims(i,:)))
    clim([0 max(data_dfof_avg(:))])
    colormap(gray)
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')
fprintf('Neuropil mask generated\n')

% clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir

%% Alternatively, load cell masks from other run
% only run instead of cell segmentation

RetImgFolder = char('002');
nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);
fprintf(['Loading masks from retinotopy runs: ' ret_str '\n'])

% loads 'mask_cell', 'mask_np'
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

% save to this folder now
save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')
fprintf('Saved loaded masks to current folder\n')

fprintf('Examine cell masks overlaid onto dF/F...\n')
figure(11);clf;
[n, n2] = subplotn(nStim);
for i = 1:nStim
    subplot(n,n2,i);
    shade_img = imShade(data_dfof_avg(:,:,i), mask_cell);
    imagesc(shade_img)
    title(num2str(Stims(i,:)))
    clim([0 max(data_dfof_avg(:))])
    colormap(gray)
end

% dfof_max
figure(12);clf;
shade_img = imShade(data_dfof_max, mask_cell);
imagesc(shade_img)
title('Max projection')
clim([0 max(data_dfof_max(:))])
colormap(gray)

% dfof average all trials
figure(13);clf;
shade_img = imShade(mean(data_dfof,3), mask_cell);
imagesc(shade_img)
title('Max projection')
clim([0 max(data_dfof_max(:))])
colormap(gray)
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVmaskoverlay.pdf']), '-dpdf')

%% Get time courses, including neuropil subtraction

fprintf('\nBegin time course extraction...\n')
down = 5; % downsample rate
fprintf(['Downsampling at M=' num2str(down) '\n'])
data_reg_down  = stackGroupProject(data_reg,down);

fprintf('Extracting cell signal...\n')
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
fprintf([num2str(nCells) ' total cells extracted\n'])

fprintf('Extracting neuropil signal...\n')
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(size(data_reg_down,3), nCells);
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
npSig_tc = (tcRemoveDC(np_tc).*np_w);
clear data_reg_down

fprintf('Neuropil subtraction complete, saving data...\n')

save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')

% convert npSub_tc to trials
fprintf('Calculating dF/F time courses...\n')
%get dF/F
nCells = size(npSub_tc,2);
tc_mat = zeros(nITI+nDelay+nOn+nOff, nCells, max_trial);
fulltc_mat = zeros(nITI+nDelay+nOn+nOff, 4, max_trial); % 1 for overall image, 2 for cells average, 3+4 for neuropil average
full_tc = squeeze(mean(mean(data_reg,1),2));
% reshape into trials
% could do this in one line with reshape and take out chunk?
% split into trials sequentially
framesOffset = int64(-nITI-nDelay:nOn+nOff-1); %gives offset indices around cOn
for iTr = 1:max_trial
    cOn = input.cStimOn{iTr};
    tc_mat(:,:,iTr) = npSub_tc(cOn+framesOffset,:);
    fulltc_mat(:,1,iTr) = full_tc(cOn+framesOffset);
    fulltc_mat(:,2,iTr) = mean(tc_mat(:,:,iTr),2);
    fulltc_mat(:,3,iTr) = mean(npSig_tc(cOn+framesOffset,:),2);
    fulltc_mat(:,4,iTr) = mean(np_tc(cOn+framesOffset,:),2);
end
tc_f = mean(tc_mat(nITI-29:nITI,:,:),1); % use last 30 frames as baseline F
fulltc_f = mean(fulltc_mat(nITI-29:nITI,:,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
fulltc_dfof = (fulltc_mat - fulltc_f) ./ fulltc_f;
clear tc_mat tc_f fulltc_mat fulltc_f

%data_fOn = mean(data_tr(:,:,nITI+nDelay+(1:nOn),:),3);

% save tc_mat??
fprintf('Time course extraction complete.\n')

%% heat map of cell responses for all trials

% examine heat map of cell responses over all trials
% this allows examination of long-term effects over trials (ex: eye goop)
% responses are average of all stim on time (nITI+nDelay+(1:nOn)) for a trial
fprintf('Examine cell responses across all trials\n')
trialRespMap = squeeze(mean(tc_dfof(nITI+nDelay+(1:nOn),:,:),1));
maxR = max(trialRespMap(:));
trialRespMap = [trialRespMap;siz_mat/max(siz_mat)*maxR;con_mat/max(con_mat)*maxR]; %add rows of normalized siz and con
figure;
imagesc(trialRespMap);
title('Avg dF/F resp by cell,trial')
xlabel('Trial')
ylabel('Cell')
colorbar

fprintf('Re-order by stimulus\n')
trialRespMapSort = sortrows(trialRespMap',[nCells+1 nCells+2])';
figure;
imagesc(trialRespMapSort);
title('Avg dF/F resp by cell,stimulus')
xlabel('Stim (El,Az)')
ylabel('Cell')
colorbar

%% calculate tuning mat and plot

window = 1:nOn; % this will be the window to define dF/F response (offset by nITI + nDelay)

tt= (-nITI-nDelay:nOn+nOff-1)*(1000./frame_rate);

fprintf('\nPlotting timecourses and measuring stimOn response\n')
tuning_mat = zeros(nStim, 2, nCells); % conditions x 2(mean+stdev) x cells
fulltuning_mat = zeros(nStim, 2, 4);
if nCells<36
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind = intersect(Ind_struct(iStim).all_trials,1:max_trial);
        plot(tt', squeeze(mean(tc_dfof(:,iCell,ind),3)))
        hold on
        tuning_mat(iStim,1,iCell) = mean(mean(tc_dfof(nITI+nDelay+(window),iCell,ind),1),3);
        tuning_mat(iStim,2,iCell) = std(mean(tc_dfof(nITI+nDelay+(window),iCell,ind),1),[],3)./sqrt(length(ind));
    end
    ylim([-0.05 0.25])
    vline(0)
    start = start + 1;
end
for iCell = 1:4
    for iStim = 1:nStim
        ind = intersect(Ind_struct(iStim).all_trials,1:max_trial);
        fulltuning_mat(iStim,1,iCell) = mean(mean(fulltc_dfof(nITI+nDelay+(window),iCell,ind),1),3);
        fulltuning_mat(iStim,2,iCell) = std(mean(fulltc_dfof(nITI+nDelay+(window),iCell,ind),1),[],3)./sqrt(length(ind));
    end
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

fprintf('Plotting size tuning curves\n')
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        set(gcf, 'Position', [0 0 800 1000]);
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    ret_mat = reshape(tuning_mat(:,1,iCell), [nCon nSiz]);
    ret_mat = ret_mat';
    plot(sizs,ret_mat)
    colormap gray
    %clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
    %clim([0 chop(max(tuning_mat(:,1,iCell),[],1),2)])
    title(num2str(chop(max(tuning_mat(:,1,iCell),[],1),2)))
    start = start +1;
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']), 'tc_dfof', 'tuning_mat', 'Stims', 'Ind_struct')

% plot tc and ret_mat for full and cell avg
figure;
subplot(2,2,1)
ret_mat = reshape(fulltuning_mat(:,1,1), [nCon nSiz]);
ret_mat = ret_mat';
plot(sizs,ret_mat)
colormap gray
title(['Full image, max:' num2str(chop(max(fulltuning_mat(:,1,1),[],1),2))])
subplot(2,2,2)
ret_mat = reshape(fulltuning_mat(:,1,2), [nCon nSiz]);
ret_mat = ret_mat';
plot(sizs,ret_mat)
colormap gray
title(['All cells average, max:' num2str(chop(max(fulltuning_mat(:,1,2),[],1),2))])
subplot(2,2,3)
ret_mat = reshape(fulltuning_mat(:,1,3), [nCon nSiz]);
ret_mat = ret_mat';
plot(sizs,ret_mat)
colormap gray
title(['Neuropil-DC, max:' num2str(chop(max(fulltuning_mat(:,1,3),[],1),2))])
subplot(2,2,4)
ret_mat = reshape(fulltuning_mat(:,1,4), [nCon nSiz]);
ret_mat = ret_mat';
plot(sizs,ret_mat)
colormap gray
title(['Neuropil, max:' num2str(chop(max(fulltuning_mat(:,1,4),[],1),2))])

% look at neuropil response
% remove top and bottom trials at each position
% look at standard deviation
% try 15-20 degree step in rough retinotopy

%% begin analysis + fitting

%% start by re-loading size tuning and retinotopy results
% check max_trial

fprintf('\nBegin results analysis\n')
% re-load size tuning data: tc_dfof, tuning_mat, Stims, Ind_struct, input,
% data_dfof_max, data_dfof_avg, nStim
fprintf('Re-loading size tuning results...\n')
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']), 'tc_dfof', 'tuning_mat', 'Stims', 'Ind_struct')
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg', 'nStim')

% re-load behavioral data
fprintf('Re-analyzing behavioral data...\n')
% define binary vectors for success, fail, ignore
tSuccess = strcmp(input.trialOutcomeCell,'success');
tFail = strcmp(input.trialOutcomeCell,'incorrect');
tIgnore = strcmp(input.trialOutcomeCell,'ignore');

% look at ignores to choose trial cutoff range
figure(1);clf; plot(tIgnore); title('Ignore trials')
max_trial = 320; % choose by examining tIgnore plot (behavior fig 1)

% store vectors for trial side, size, contrast
lTrial = cell2mat(input.tLeftTrial);
rTrial = 1-lTrial;
siz_mat = cell2mat(input.tGratingDiameterDeg);
con_mat = cell2mat(input.tGratingContrast);
sizs = unique(siz_mat);
cons = unique(con_mat);
nSiz = length(sizs);
nCon = length(cons);

% load defining variables
nITI = input.nFramesIti;
nDelay = input.nFramesDelay;
nOn = input.nFramesOn;
nOff = input.nFramesTooFast;
ntrials = size(input.tGratingContrast,2);

fprintf(['Truncating to ' num2str(max_trial) '/' num2str(ntrials) ' trials\n'])
% truncate lTrial, rTrial, siz_mat, con_mat to max_trial
lTrial = lTrial(1:max_trial);
rTrial = rTrial(1:max_trial);
siz_mat = siz_mat(1:max_trial);
con_mat = con_mat(1:max_trial);

fprintf('\nRe-loading retinotopy results...\n')
% load retinotopy data
RetImgFolder = char('002');
nret = size(RetImgFolder,1);
ret_str = catRunName(RetImgFolder, nret);

% loads 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind'
% lbub_fits contains RF fit data, extract RF centers from here
% goodfit_ind
fprintf(['Loading fits from retinotopy runs: ' ret_str '\n'])
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
load(fn_out);
cellAz = lbub_fits(:,4,4);
cellEl = lbub_fits(:,5,4);
fprintf('Retinotopy fits loaded, found cell receptive field coordinates\n')
nCells = length(goodfit_ind);
fprintf(['# goodfit cells = ' num2str(nCells) '\n'])

% input stimulus location based on experimental choice
stimEl = double(input.gratingElevationDeg); % should be 0
stimAz = double(input.gratingEccentricityDeg); % should be 30
fprintf(['Stimulus at: El ' num2str(stimEl) ', Az ' num2str(stimAz) '\n'])

fprintf('Calculating cell RF distances to stimulus...\n')
cellDists = sqrt((cellAz(goodfit_ind)-stimAz).^2+(cellEl(goodfit_ind)-stimEl).^2);

figure(1);clf;
scatter(cellAz(goodfit_ind),cellEl(goodfit_ind),10,cellDists)
hold on
plot(stimAz,stimEl,'rx')
hold off
title('RF location for goodfit cells')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
nasalInd = find(cellAz(goodfit_ind)<15); % maybe use 10

% compare distance to index cells into distance ranges
fprintf('Binning cells by RF center distance\n')
bIndex = 0*cellDists;
nBin = 3;
for i=1:nCells
    if cellDists(i)<7
        bIndex(i) = 1;
    elseif cellDists(i)<15
        bIndex(i) = 2;
    else
        bIndex(i) = 3;
    end
end
nBinCells = [sum(bIndex==1) sum(bIndex==2) sum(bIndex==3)];
fprintf([num2str(nBinCells) ' cells selected\n'])

% binStrs = ["0-6","6-12","12+"];
% cBins = categorical({'0-6','6-12','12+'});
% cBins = reordercats(cBins,{'0-6','6-12','12+'});
binStrs = ["0-7","7-15","15+"];
cBins = categorical({'0-7','7-15','15+'});
cBins = reordercats(cBins,{'0-7','7-15','15+'});
% binStrs = ["0-10","10-20","20+"];
% cBins = categorical({'0-10','10-20','20+'});
% cBins = reordercats(cBins,{'0-10','10-20','20+'});

% histogram of # cells in each bin
figure(2);clf;
subplot(2,1,1)
histogram(cellDists,10)
title('RF distance histogram')
xlabel('RF-stim distance')
ylabel('# cells')
subplot(2,1,2)
bar(cBins,nBinCells)
title('Cell Counts by RF Distance Bin')
xlabel('Distance Bin')
ylabel('# cells')

fprintf('Data re-loaded, continue with size-tuning analysis\n')

%% now select window to extract response

tt= (-nITI-nDelay:nOn+nOff-1)*(1000./frame_rate);
window = 13:nOn; % this will be the window to define dF/F response (offset by nITI + nDelay)

fprintf('\nExamine average cell dF/F timecourse to select response window\n')
avTC = squeeze(mean(mean(tc_dfof(:,goodfit_ind,:),2),3)); % average all cells,trials
figure(3);
for i = 1:length(goodfit_ind)
%     for j = 1:nStim
%         iSiz = find(Stims(iStim,1) == sizs);
%         iCon = find(Stims(iStim,2) == cons);
%         ind = Ind_struct(iStim).right_trials; %right trial inds for this stim
%         ind = intersect(ind,1:max_trial); %truncate to max_trial
%         ind = setdiff(ind,find(tIgnore)); %remove ignore trials
%         
%         plot(tt',squeeze(mean(tc_dfof(:,goodfit_ind(i),ind),3)))
%         hold on
%     end
    plot(tt', squeeze(mean(tc_dfof(:,goodfit_ind(i),:),3)), 'LineWidth',1)
    hold on
end
plot(tt', avTC, 'LineWidth', 5)
title('Average cell timecourses')
ylim([-0.05 0.25])
vline(0)
vline(tt(nITI+nDelay+window(end)))
vline(tt(nITI+nDelay+window(1)))
h=vfill([tt(nITI+nDelay+window(1)),tt(nITI+nDelay+window(end))],'gray','FaceAlpha',0.5);
uistack(h,'bottom')


%% size tuning analysis
window = 13:nOn; % this will be the window to define dF/F response (offset by nITI + nDelay)

% build size tuning data (with Mean and SEM) for each bin of goodfit cells
sizeTune = cell(nSiz,nCon,nCells);
sizeMean = zeros(nSiz,nCon,nCells);
sizeSEM = sizeMean;
sizeTuneL = sizeTune;
sizeMeanL = sizeMean;
sizeSEML = sizeSEM;
sizeTuneRC = sizeTune;
sizeMeanRC = sizeMean;
sizeSEMRC = sizeSEM;
sizeTuneRI = sizeTune;
sizeMeanRI = sizeMean;
sizeSEMRI = sizeSEM;
sizeTuneLC = sizeTune;
sizeMeanLC = sizeMean;
sizeSEMLC = sizeSEM;
sizeTuneLI = sizeTune;
sizeMeanLI = sizeMean;
sizeSEMLI = sizeSEM;

% calculate size tuning data for all goodfit cells
% sorts left and right side trials
for i = 1:nCells
    iCell = goodfit_ind(i);
    
    for iStim = 1:nStim
        iSiz = find(Stims(iStim,1) == sizs);
        iCon = find(Stims(iStim,2) == cons);
        
        indR = intersect(Ind_struct(iStim).right_trials,1:max_trial); %right trial inds for this stim, and truncate to max trial
        indR = setdiff(indR,find(tIgnore)); %remove ignore trials
        indRC = intersect(indR,find(tSuccess));
        indRI = intersect(indR,find(tFail));
        indL = intersect(Ind_struct(iStim).left_trials,1:max_trial); %right trial inds for this stim, and truncate to max trial
        indL = setdiff(indL,find(tIgnore)); %remove ignore trials
        indLC = intersect(indL,find(tSuccess));
        indLI = intersect(indL,find(tFail));
        
        % take mean dF/F during stimOn at selected indices
        stimOnR = mean(tc_dfof(nITI+nDelay+(window),iCell,indR),1);
        sdR = std(stimOnR);
        stimOnRC = mean(tc_dfof(nITI+nDelay+(window),iCell,indRC),1);
        sdRC = std(stimOnRC);
        stimOnRI = mean(tc_dfof(nITI+nDelay+(window),iCell,indRI),1);
        sdRI = std(stimOnRI);
        stimOnL = mean(tc_dfof(nITI+nDelay+(window),iCell,indL),1);
        sdL = std(stimOnL);
        stimOnLC = mean(tc_dfof(nITI+nDelay+(window),iCell,indLC),1);
        sdLC = std(stimOnLC);
        stimOnLI = mean(tc_dfof(nITI+nDelay+(window),iCell,indLI),1);
        sdLI = std(stimOnLI);
        
        % this cell matrix stores all stimOn for all trials in ind_all
        % dims (szs, cons, cells)
        sizeTune{iSiz, iCon, i} = squeeze(stimOnR);
        % also take mean and SEM
        sizeMean(iSiz,iCon,i) = mean(stimOnR,3);
        sizeSEM(iSiz,iCon,i) = sdR/sqrt(length(indR));
        % for left trials too
        sizeTuneL{iSiz, iCon, i} = squeeze(stimOnL);
        sizeMeanL(iSiz,iCon,i) = mean(stimOnL,3);
        sizeSEML(iSiz,iCon,i) = sdL/sqrt(length(indL));
        sizeTuneRC{iSiz, iCon, i} = squeeze(stimOnRC);
        sizeMeanRC(iSiz,iCon,i) = mean(stimOnRC,3);
        sizeSEMRC(iSiz,iCon,i) = sdRC/sqrt(length(indRC));
        sizeTuneRI{iSiz, iCon, i} = squeeze(stimOnRI);
        sizeMeanRI(iSiz,iCon,i) = mean(stimOnRI,3);
        sizeSEMRI(iSiz,iCon,i) = sdRI/sqrt(length(indRI));
        sizeTuneLC{iSiz, iCon, i} = squeeze(stimOnLC);
        sizeMeanLC(iSiz,iCon,i) = mean(stimOnLC,3);
        sizeSEMLC(iSiz,iCon,i) = sdLC/sqrt(length(indLC));
        sizeTuneLI{iSiz, iCon, i} = squeeze(stimOnLI);
        sizeMeanLI(iSiz,iCon,i) = mean(stimOnLI,3);
        sizeSEMLI(iSiz,iCon,i) = sdLI/sqrt(length(indLI));
    end
end

% vvv commented out for now
% plot all size tuning curves from each bin
for k = 1:nBin
%     binCells = find(bIndex==k);
%     [n, n2] = subplotn(nBinCells(k));
%     figure(k+1);clf;
%     suptitle(['Right Side Size Tuning, ' char(binStrs(k))])
%     for i = 1:nBinCells(k)
%         iCell = binCells(i);
%         subplot(n,n2,i)
%         for iCon = 1:nCon
%             errorbar(sizs,sizeMean(:,iCon,iCell),sizeSEM(:,iCon,iCell))
%             hold on
%         end
%         ylim([-1.2*(max(max(sizeSEM(:,:,iCell)))) 1.2*(max(max(sizeMean(:,:,iCell)))+max(max(sizeSEM(:,:,iCell))))])
%         title(['Cell #: ' num2str(iCell)])
%         hold off
%     end
%     % Construct a Legend with the data from the sub-plots
%     hL = legend(num2str(cons'));
%     % Programatically move the Legend
%     newPosition = [0.5 0.07 0.2 0.2];
%     newUnits = 'normalized';
%     set(hL, 'Position', newPosition, 'Units', newUnits);
end
% repeat for left
for k = 1:nBin
%     binCells = find(bIndex==k);
%     [n, n2] = subplotn(nBinCells(k));
%     figure(k+5);clf;
%     suptitle(['Left Side Size Tuning, ' char(binStrs(k))])
%     for i = 1:nBinCells(k)
%         iCell = binCells(i);
%         subplot(n,n2,i)
%         for iCon = 1:nCon
%             errorbar(sizs,sizeMeanL(:,iCon,iCell),sizeSEML(:,iCon,iCell))
%             hold on
%         end
%         ylim([-1.2*(max(max(sizeSEM(:,:,iCell)))) 1.2*(max(max(sizeMean(:,:,iCell)))+max(max(sizeSEM(:,:,iCell))))])
%         title(['Cell #: ' num2str(iCell)])
%         hold off
%     end
%     % Construct a Legend with the data from the sub-plots
%     hL = legend(num2str(cons'));
%     % Programatically move the Legend
%     newPosition = [0.5 0.07 0.2 0.2];
%     newUnits = 'normalized';
%     set(hL, 'Position', newPosition, 'Units', newUnits);
end

% save responses as sizeTuneData.mat
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
save(filename, 'sizeTune', 'sizeMean', 'sizeSEM', 'cellDists')
fprintf('Saved sizeTuneData.mat\n')

%% examine average size tuning curves within bins, for different conditions
% last make a figure of mean size tuning curve per bin
figure(5);clf;
%suptitle('Average Size Tuning Curve by Bin')
for k = 1:nBin
    subplot(2,nBin,k)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMean(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEM(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(right) RF dist ' char(binStrs(k)) ', n=' num2str(nBinCells(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
    
    subplot(2,nBin,k+nBin)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanL(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEML(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(left) RF dist ' char(binStrs(k)) ', n=' num2str(nBinCells(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
end

% now only correct
figure(6);clf;
%suptitle('Average Size Tuning Curve by Bin')
for k = 1:nBin
    subplot(2,nBin,k)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanRC(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEMRC(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(right, success) RF dist ' char(binStrs(k)) ', n=' num2str(nBinCells(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
    
    subplot(2,nBin,k+nBin)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanLC(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEMLC(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(left, success) RF dist ' char(binStrs(k)) ', n=' num2str(sum(bIndex==k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
end

% now only incorrect
figure(7);clf;
%suptitle('Average Size Tuning Curve by Bin')
for k = 1:nBin
    subplot(2,nBin,k)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanRI(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEMRI(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(right, fail) RF dist ' char(binStrs(k)) ', n=' num2str(nBinCells(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
    
    subplot(2,nBin,k+nBin)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanLI(:,iCon,find(bIndex == k)),3),1/sqrt(nBinCells(k))*geomean(sizeSEMLI(:,iCon,find(bIndex == k)),3))
        hold on
    end
    title(['(left, fail) RF dist ' char(binStrs(k)) ', n=' num2str(nBinCells(k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    ylim([-0.05 0.12])
    legend(num2str(cons'),'Location','best')
end

% average only nasal RF cells, using both correct and incorrect
figure(8);clf;
%suptitle('Average Size Tuning Curve by Bin')
for k = 1:nBin
    subplot(2,nBin,k)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMean(:,iCon,intersect(find(bIndex == k),nasalInd)),3),1/sqrt(nBinCells(k))*geomean(sizeSEM(:,iCon,intersect(find(bIndex == k),nasalInd)),3))
        hold on
    end
    title(['(right) RF dist ' char(binStrs(k)) ', n=' num2str(sum(bIndex(nasalInd)==k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    legend(num2str(cons'),'Location','best')
    
    subplot(2,nBin,k+nBin)
    for iCon=1:nCon
        errorbar(sizs,mean(sizeMeanL(:,iCon,intersect(find(bIndex == k),nasalInd)),3),1/sqrt(nBinCells(k))*geomean(sizeSEML(:,iCon,intersect(find(bIndex == k),nasalInd)),3))
        hold on
    end
    title(['(left) RF dist ' char(binStrs(k)) ', n=' num2str(sum(bIndex(nasalInd)==k))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    ylim([-0.1 0.25])
    legend(num2str(cons'),'Location','best')
end

%% examine individual cells with left vs right
% also include incorrect vs correct
% look just at cell 34
for iCell=1:nCells
    if ~sum(nasalInd==iCell)
        continue
    end
    figure(9);clf;
    %suptitle('Average Size Tuning Curve by Bin')
    ax1=subplot(1,2,2);
    for iCon=1:nCon
        errorbar(sizs,sizeMean(:,iCon,iCell),1/sqrt(nBinCells(k))*sizeSEM(:,iCon,iCell))
        hold on
    end
    title(['(right) Cell# ' num2str(iCell)])
    xlabel('Size (deg)')
    ylabel('dF/F')
    legend(num2str(cons'),'Location','best')
    ax2=subplot(1,2,1);
    for iCon=1:nCon
        errorbar(sizs,sizeMeanL(:,iCon,iCell),1/sqrt(nBinCells(k))*sizeSEML(:,iCon,iCell))
        hold on
    end
    title(['(left) Cell# ' num2str(iCell)])
    xlabel('Size (deg)')
    ylabel('dF/F')
    linkaxes([ax1 ax2],'y')
    legend(num2str(cons'),'Location','best')
    pause
end

%% same but separate correct + incorrect
for iCell=1:nCells
    figure(9);clf;
    %suptitle('Average Size Tuning Curve by Bin')
    ax1=subplot(1,2,2);
    for iCon=1:nCon
        errorbar(sizs,sizeMean(:,iCon,iCell),1/sqrt(nBinCells(k))*sizeSEM(:,iCon,iCell))
        hold on
    end
    title(['(right) Cell# ' num2str(iCell)])
    xlabel('Size (deg)')
    ylabel('dF/F')
    legend(num2str(cons'),'Location','best')
    ax2=subplot(1,2,1);
    for iCon=1:nCon
        errorbar(sizs,sizeMeanL(:,iCon,iCell),1/sqrt(nBinCells(k))*sizeSEML(:,iCon,iCell))
        hold on
    end
    title(['(left) Cell# ' num2str(iCell)])
    xlabel('Size (deg)')
    ylabel('dF/F')
    linkaxes([ax1 ax2],'y')
    legend(num2str(cons'),'Location','best')
    pause
end

%% fitting

%% Fit size tuning curves
chosen = 1:nCells;
%chosen = [];

% First model: single sigmoid (3 fit params)
% use 90% cutoff as measure of preferred size, no suppression index
% Second model: sum of + and - sigmoids (6 fit params)
% use peak of fit as preferred size, and measure suppression index

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))))
%logfit1 = @(coefs,xdata) 2*coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3)))) - coefs(1)
% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))))

szRng = linspace(0.1,max(sizs));

% create a struct to store variables
s = zeros(nCells, nCon);
s3 = zeros(nCells, nCon,3);
s6 = zeros(nCells, nCon,6);
sizeFits = struct('Rsq1',s,'coefs1',s3,'Rsq2',s,'coefs2',s6,'prefSize',s,'suppInd',s,'Fscore',s,'Ftest',s);

cd('K:\Code')
opts = optimset('Display','off');

figure(6);clf;
for i = 1:nCells
    if ~bIndex(i)
        continue
    end
    fprintf(['\nCell# ' num2str(i) ', bIndex ' num2str(bIndex(i))])
    
    for iCon = 1:nCon
        fprintf(['\nContrast: ' num2str(cons(iCon))])
        
        dumdum = zeros(1,8);%[0]; % define zero point for dF/F
        szs0 = zeros(1,8);%[0.1]; % and size
        nPts = 8;
        for iSz = 1:nSiz
            nPts = nPts + length(sizeTune{iSz,iCon,i});
            dumdum = [dumdum sizeTune{iSz,iCon,i}'];
            szs0 = [szs0 sizs(iSz)*ones(1,length(sizeTune{iSz,iCon,i}))];
        end
        
        SStot = sum((dumdum-mean(dumdum)).^2);
        fprintf([' (SStot: ' num2str(SStot) ')'])
        
        % make initial guesses
        maxMean = max(sizeMean(:,iCon,i));
        
        % single sigmoid model
        fprintf('\nEvaluating fit 1...\n')
        guess1 = [0.5*maxMean 0.35 10]; % guesses for 3 params (Ae, ke, xe)
        lb = [0 0 0];
        ub = [inf 1 80]; %[inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        %[fitx,resnorm1] = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        %guess = lsqcurvefit(logfit1,guess,szs0,dumdum,lb,ub,opts);
        [fitx,OF1] = fminsearchbnd(@(c)sigevalpen1(c,szs0,dumdum),guess1,lb,ub,opts);
        resnorm1 = sigeval1(fitx,szs0,dumdum);
        
        fprintf('Fit 1 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)) ', center: ' num2str(fitx(3))])
        r2 = 1 - resnorm1/SStot; % resnorm is SSE
        fprintf(['\nFit 1 R-sq:' num2str(r2) '\n'])
        fprintf(['Fit 1 OF:' num2str(OF1) ', SSE:' num2str(resnorm1) ', Pen:' num2str(OF1-resnorm1) '\n'])
        sizeFits.Rsq1(i, iCon) = r2;
        sizeFits.coefs1(i,iCon,:) = fitx;
        
        fitout1 = logfit1(fitx,szRng);
        maxResp1 = max(fitout1);
        prefSize1 = szRng(find(fitout1>(0.9*maxResp1),1));
        if ~numel(prefSize1); prefSize1 = 1; end;
        suppInd1 = 0;
        fprintf(['Pref size: ' num2str(prefSize1) ' (dF/F: ' num2str(0.9*maxResp1) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd1) ' (no suppression in single sigmoid model)\n'])
        
        % double sigmoid model
        fprintf('\nEvaluating fit 2...\n')
        guess2 = [maxMean 0.4 10 maxMean-mean(dumdum) 0.3 20]; % guesses for (Ae, k1, x1, Ai, k2 x2)
        lb = [0 0 0 0 0 0];
        ub = [inf 1 80 inf 1 80]; %[inf inf inf inf inf inf]; % here limit steepness to 1 (above this is way too steep)
        %[fitx,resnorm2] = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        %guess = lsqcurvefit(logfit2,guess,szs0,dumdum,lb,ub,opts);
        [fitx,OF2] = fminsearchbnd(@(c)sigevalpen2(c,szs0,dumdum),guess2,lb,ub,opts);
        resnorm2 = sigeval2(fitx,szs0,dumdum);
        
        fprintf('Fit 2 readout:\n')
        fprintf(['Ex amp: ' num2str(fitx(1)) ', steepness: ' num2str(fitx(2)+fitx(5)) ', center: ' num2str(fitx(3))])
        fprintf(['\nInh amp: ' num2str(fitx(4)) ', steepness: ' num2str(fitx(5)) ', center: ' num2str(fitx(3)+fitx(6))])
        r2 = 1 - resnorm2/SStot; % resnorm is SSE
        fprintf(['\nFit 2 R-sq: ' num2str(r2) '\n'])
        fprintf(['Fit 2 OF:' num2str(OF2) ', SSE:' num2str(resnorm2) ', Pen:' num2str(OF2-resnorm2) '\n'])
        sizeFits.Rsq2(i,iCon) = r2;
        sizeFits.coefs2(i,iCon,:) = fitx;
        
        fitout2 = logfit2(fitx,szRng);
        maxResp2 = max(fitout2);
        prefSize2 = szRng(find(fitout2==maxResp2,1));
        suppInd2 = 1 - fitout2(end)/maxResp2;
        fprintf(['Pref size: ' num2str(prefSize2) ' (dF/F: ' num2str(maxResp2) ')\n'])
        fprintf(['Suppression Index: ' num2str(suppInd2) '\n'])
        
        % plot both fits, if chosen cell
        if sum(i==chosen)
            subplot(2,nCon,iCon)
            errorbar([0 sizs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout1,'-')
            plot(szRng,logfit1(guess1,szRng),'g--')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq1(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('mean','data','fit')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
            subplot(2,nCon,iCon+nCon)
            errorbar([0 sizs],[0 sizeMean(:,iCon,i)'],[0 sizeSEM(:,iCon,i)'])
            hold on
            plot(szs0,dumdum,'.b')
            plot(szRng,fitout2,'-')
            plot(szRng,logfit2(guess2,szRng),'g--')
            hold off
            ylim([min([-0.5*maxResp1 min(sizeMean(:,iCon,i))]) 1.4*max([maxResp2 max(sizeMean(:,iCon,i))])])
            title(['Cell #' num2str(i) ', Con ' num2str(cons(iCon)) ', R^2=' num2str(sizeFits.Rsq2(i,iCon))])
            xlabel('Stimulus Size')
            ylabel('dF/F')
            legend('mean','data','fit')
            
            h = findobj(gca,'Type','errorbar');
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        
        % now F-test
        fprintf([num2str(nPts) ' points\n'])
        dfS = nPts - 3; % single sigmoid model df = #sizes - 3 model parameters
        dfD = nPts - 6; % double sigmoid model df = #sizes - 6 model parameters
        Fcrit = finv(0.95,dfD,dfS); % measure critical F at a=0.05
        sizeFits.Fscore(i,iCon) = (resnorm2 - resnorm1)/(dfD-dfS) ./ (resnorm2/dfD);
        Ftest = sizeFits.Fscore(i,iCon)>Fcrit;
        sizeFits.Ftest(i,iCon) = Ftest;
        
        % store values from selected test (0=single,1=double)
        if Ftest
            sizeFits.prefSize(i,iCon) = prefSize2;
            sizeFits.suppInd(i,iCon) = suppInd2;
        else
            sizeFits.prefSize(i,iCon) = prefSize1;
            sizeFits.suppInd(i,iCon) = suppInd1;
        end
    end
    
    fprintf(['\nF-tests for cell i=' num2str(i) ', distIndex ' num2str(bIndex(i)) ':\nF-scores: ' num2str(sizeFits.Fscore(i,:)) '\nChoose model 2? ' num2str(sizeFits.Ftest(i,:)) '\n'])
    
    if sum(i==chosen)
        pause
    end
end

% Save fit results as sizeFitResults.mat
% save sizeFitResults.mat, with sizeFits struct, bIndex
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
%save(filename, 'sizeFits', 'bIndex')
fprintf('Saved sizeFitResults_SP.mat\n')

%% Load instead of fitting
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
load(filename)
fprintf('Loaded sizeTuneData.mat\n')

filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
load(filename) %'allFtest', 'allPrefSize', 'allSuppInd') %, 'bIndex')
fprintf('Loaded sizeFitResults_SP.mat\n')

% WARNING: bIndex is calculated based on predefined bin ranges (right now 0-6, 6-12, 12+)
% so better to use cellDists and recalculate bins if they change

%% examine parameters
bin1inds = find(bIndex == 1);
bin2inds = find(bIndex == 2);
bin3inds = find(bIndex == 3);

figure(3);clf;
subplot(1,3,1)
histogram(sizeFits.prefSize(bin1inds,:),50)
title('PrefSize (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(sizeFits.prefSize(bin2inds,:),50)
title('PrefSize (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(sizeFits.prefSize(bin3inds,:),50)
title('PrefSize (RF dist 15+ deg)')

figure(4);clf;
subplot(1,3,1)
histogram(sizeFits.suppInd(bin1inds,:),30)
title('SuppInd (RF dist 0-7 deg)')
subplot(1,3,2)
histogram(sizeFits.suppInd(bin2inds,:),30)
title('SuppInd (RF dist 7-15 deg)')
subplot(1,3,3)
histogram(sizeFits.suppInd(bin3inds,:),30)
title('SuppInd (RF dist 15+ deg)')

% allRsq1 refers to single model, allRsq2 refers to double model
figure(5);clf;
subplot(1,3,1)
histogram(sizeFits.Rsq1(:),30)
ylim([0 20])
title('R^2 for single model')
subplot(1,3,2)
histogram(sizeFits.Rsq2(:),30)
ylim([0 20])
title('R^2 for double model')
subplot(1,3,3)
histogram(sizeFits.Rsq2(:)-sizeFits.Rsq1(:),30)
ylim([0 20])
title('Paired Difference (double-single)')

% insert adjusted Rsq *********

% now examine distance vs parameters, to pick cutoff
% ke: coefs1(2), coefs2(2)+coefs2(5)
% xe: coefs1(3),coefs2(3)
% ki: coefs2(5)
% xi: coefs2(3)+coefs2(6)
x = repmat(cellDists,1,nCon)
figure(6);clf;
subplot(2,3,1)
histogram(sizeFits.coefs1(:,:,2),20)
% plot(x,sizeFits.coefs1(:,:,2),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - single model')
subplot(2,3,2)
histogram(sizeFits.coefs2(:,:,2)+sizeFits.coefs2(:,:,5),20)
% plot(x,sizeFits.coefs2(:,:,2)+sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_e')
title('k_e - double model')
subplot(2,3,3)
histogram(sizeFits.coefs2(:,:,5),20)
% plot(x,sizeFits.coefs2(:,:,5),'o')
% xlabel('cell dist (deg)')
% ylabel('k_i')
title('k_i - double model')
subplot(2,3,4)
histogram(sizeFits.coefs1(:,:,3),20)
% plot(x,sizeFits.coefs1(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - single model')
subplot(2,3,5)
histogram(sizeFits.coefs2(:,:,3),20)
% plot(x,sizeFits.coefs2(:,:,3),'o')
% xlabel('cell dist (deg)')
% ylabel('x_e')
title('x_e - double model')
subplot(2,3,6)
histogram(sizeFits.coefs2(:,:,3)+sizeFits.coefs2(:,:,6),20)
% plot(x,sizeFits.coefs2(:,:,3)+sizeFits.coefs2(:,:,6),'o')
% xlabel('cell dist (deg)')
% ylabel('x_i')
title('x_i - double model')

%% examine overall data after F-test selection
cCons = categorical(cons);

% plot the number of unsuppressed vs suppressed cells
figure(7);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = nBinCells(k);
    subplot(1,3,k)
    
    binFtest = zeros(nCon,2); % one column each for unsupp, supp cells
    if nCutCells
        for iCon = 1:nCon
            binFtest(iCon,1) = nCutCells - sum(sizeFits.Ftest(binCells,iCon));
            binFtest(iCon,2) = sum(sizeFits.Ftest(binCells,iCon));
        end
    end
    bar(cCons,binFtest)
    title(['Unsupp vs Supp cells, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('# cells')
    legend('Unsupp','Supp')
end

figure(8);clf;
for k=1:3
    binCells = find(bIndex==k);
    nCutCells = length(binCells);
    
    % first show average PrefSize at each contrast
    subplot(2,3,k)
    prefAv = mean(sizeFits.prefSize(binCells,:),1);
    prefSEM = std(sizeFits.prefSize(binCells,:),0,1)/sqrt(nCutCells);
    bar(cCons,prefAv)
    hold on
    errorbar(cCons,prefAv,prefSEM,'.')
    hold off
    title(['Preferred Size, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Pref Size')
    ylim([0 80])
    
    % second show average SI at each contrast
    subplot(2,3,k+3)
    SIav = mean(sizeFits.suppInd(binCells,:),1);
    SISEM = std(sizeFits.suppInd(binCells,:),0,1)/sqrt(nCutCells);
    bar(cCons,SIav)
    hold on
    errorbar(cCons,SIav,SISEM,'.')
    hold off
    title(['Suppression Index, RF dist: ' char(binStrs(k)) ' (n=' num2str(nCutCells) ' cells)'])
    xlabel('Contrast')
    ylabel('Average Supp Ind')
    ylim([0 0.6])
end

%% Examine how prefSize changes with contrast
% will need to store fit params for examining a single size

% first calculate slopes, fitting with a 1st degree polynomial (line)
lineFits = zeros(nCells,2);
for i=1:nCells
    lineFits(i,:) = polyfit(cons,(sizeFits.prefSize(i,:)/(sizeFits.prefSize(i,end))),1); %allPrefSize(i,:)
end

% the plot each cell's prefSize vs con curve, within each bin
for k=1:nBin
    figure(8+k);clf;
    [n, n2] = subplotn(nBinCells(k));
    for i = 1:nBinCells(k)
        if k==1
            iCell = bin1inds(i);
        elseif k==2
            iCell = bin2inds(i);
        else
            iCell = bin3inds(i);
        end
        subplot(n,n2,i)
        plot(cons,(sizeFits.prefSize(iCell,:)/(sizeFits.prefSize(iCell,end))),'o') %allPrefSize(iCell,:)
        hold on
        plot(cons,polyval(lineFits(iCell,:),cons),'r-')
        hold off
        title(['Cell #' num2str(iCell) ', m=' num2str(lineFits(iCell,1))])
        xlabel('Contrast')
        ylabel('Pref Sz')
        ylim([0 5])
    end
end

% readout means and stdevs for each bin
for k=1:nBin
    fprintf(['Bin #' num2str(k) ', Mean slopes: ' num2str(mean(lineFits(bIndex==k,1),1)) ', stdev: ' num2str(std(lineFits(bIndex==k,1),[],1)) '\n'])
end

dataset = 'FS_HVA';
eval(dataset);
iexp = 11;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%% load data
tic
data = [];
clear temp
offset = 0;

%load FS data
nrun = size(expt(iexp).runs,1);
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' expt(iexp).date '_i' expt(iexp).mouse '\' expt(iexp).runs(irun,:)];
    cd(CD);
    imgMatFile = [expt(iexp).runs(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-i' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).time_mat(irun,:) '.mat'];
    load(fName);
    
    if isnan(expt(iexp).frames(irun))
        nframes = [input.counterValues{end}(end) info.config.frames];
    else
        nframes = expt(iexp).frames(irun);
    end
    if nframes<input.counterValues{end}(end)
        ntrials = size(input.trialOutcomeCell,2);
        for itrial = ntrials:-1:1
            if input.counterValues{itrial}(end) <= nframes
                break
            end
        end
        input = trialChopper(input,[1 itrial]);
    end
    temp(irun) = input;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                %temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                if ~isempty(temp(irun).cLeverUp{itrial})
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                else
                    temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial};
                end
                if ~isempty(temp(irun).cTargetOn{itrial})
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
                else
                    temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial};
                end
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end 

fs_input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%load direction data
CD = [LG_base '\Data\2P_images\' expt(iexp).date '_i' expt(iexp).mouse '\' expt(iexp).dirtuning];   
cd(CD);
imgMatFile = [expt(iexp).dirtuning '_000_000.mat'];
load(imgMatFile);
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-i' expt(iexp).mouse '-' expt(iexp).date '-' expt(iexp).dirtuning_time '.mat'];
load(fName);
dir_input = input;
dir_input = trialChopper(dir_input,[1 160]);
nframes = info.config.frames;
    
fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
if size(data_temp,1)== 2
    data_temp = data_temp(1,:,:,:);
end

data_temp = squeeze(data_temp);
data = cat(3,data,data_temp);

clear data_temp
clear temp

%% test stability
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

data_avg = mean(data(:,:,40001:40500),3);

%% register
run_str = ['runs-' expt(iexp).runs(1,:) '_' expt(iexp).dirtuning];
[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_input.mat']), 'input')

clear data
%% test stability
figure; plot(out); xlabel('Frames'); ylabel('pixels')

figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% FS stim responses
tCyc = cell2mat(fs_input.tCyclesOn);
cStart = cell2mat(fs_input.cFirstStim);
cStim = cell2mat(fs_input.cStimOn);
cTarget = celleqel2mat_padded(fs_input.cTargetOn);
nTrials = length(tCyc);
sz = size(data_reg);
data_f = zeros(sz(1),sz(2),nTrials);
data_base = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        data_f(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
        data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+6:cStart(itrial)+9),3);
    else
        data_f(:,:,itrial) = nan(sz(1),sz(2));
        data_base(:,:,itrial) = nan(sz(1),sz(2));
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+9 < sz(3)
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+6:cTarget(itrial)+9),3);
        else
            data_targ(:,:,itrial) = nan(sz(1),sz(2));
        end
    else
        data_targ(:,:,itrial) = nan(sz(1),sz(2));
    end
end
data_base_dfof = (data_base-data_f)./data_f;
data_targ_dfof = (data_targ-data_f)./data_f;
targCon = celleqel2mat_padded(fs_input.tGratingContrast);
if fs_input.doRandCon
    baseCon = ones(size(targCon));
else
    baseCon = celleqel2mat_padded(fs_input.tBaseGratingContrast);
end
ind_con = intersect(find(targCon == 1),find(baseCon == 0));
baseDir = celleqel2mat_padded(fs_input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
targetDelta = round(celleqel2mat_padded(fs_input.tGratingDirectionDeg),0);
deltas = unique(targetDelta);
nDelta = length(deltas);
data_dfof_dir = zeros(sz(1),sz(2),ndir);
[n n2] = subplotn(ndir);
figure;
for idir = 1:ndir
    ind = setdiff(find(baseDir == dirs(idir)),ind_con);
    data_dfof_dir(:,:,idir) = nanmean(data_base_dfof(:,:,ind),3);
    subplot(n,n2,idir)
    imagesc(data_dfof_dir(:,:,idir))
    title(dirs(idir))
end

data_dfof_targ = zeros(sz(1),sz(2),nDelta);
[n n2] = subplotn(nDelta);
figure;
for idir = 1:nDelta
    ind = find(targetDelta == deltas(idir));
    data_dfof_targ(:,:,idir) = nanmean(data_targ_dfof(:,:,ind),3);
    subplot(n,n2,idir)
    imagesc(data_dfof_targ(:,:,idir))
    title(deltas(idir))
end
data_fs_dfof = cat(3,data_dfof_dir,data_dfof_targ);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_fs_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)

%% dir tuning responses

nOn = dir_input.nScansOn;
nOff = dir_input.nScansOff;
dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
nTrials = length(dir_mat);
dir_data = data_reg(:,:,offset+1:end);
if nTrials.*(nOn+nOff)<size(dir_data,3)
    dir_data = dir_data(:,:,1:nTrials.*(nOn+nOff));
end
down = 10;
data_reg_down  = stackGroupProject(dir_data,down);

%compute dir dF/F
sz = size(data_reg_down);
nOn = dir_input.nScansOn./down;
nOff = dir_input.nScansOff./down;
data_trial = reshape(data_reg_down,[sz(1),sz(2),nOn+nOff,nTrials]);

dirs = unique(dir_mat);
ndir = length(dirs);
data_dir = zeros(sz(1),sz(2),nOn+nOff,ndir);
for idir = 1:ndir
    ind = find(dir_mat == dirs(idir));
    data_dir(:,:,:,idir) = mean(data_trial(:,:,:,ind),4);
end
data_dir_f = squeeze(mean(data_dir(:,:,ceil(nOff/2):nOff,:),3));
data_dir_avg = squeeze(mean(data_dir(:,:,nOff+1:end,:),3));
data_dir_dfof = (data_dir_avg-data_dir_f)./data_dir_f;
[n n2] = subplotn(ndir);
figure;
for idir = 1:ndir
    subplot(n,n2,idir)
    imagesc(data_dir_dfof(:,:,idir))
end
figure;
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = max(data_dir_dfof,[],3);
rgb(:,:,2) = max(data_fs_dfof,[],3);
imagesc(rgb);

%% segment cells
data_dfof_all = cat(3,data_dir_dfof,data_fs_dfof);

mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof_all;

for iStim = 1:size(data_dfof_all,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)

%% neuropil subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_mask_cell.mat']), 'data_dfof_all', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

% neuropil subtraction
down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str], [expt(iexp).date '_i' expt(iexp).mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
fprintf('\n')

%% FS analysis
if iscell(fs_input.nFramesOn)
    nOn = fs_input.nFramesOn{1};
else
    nOn = fs_input.nFramesOn;
end
prewin_frames = 30;
postwin_frames = 90;
tCyc = cell2mat(fs_input.tCyclesOn);
cStart = celleqel2mat_padded(fs_input.cFirstStim);
cTarget = celleqel2mat_padded(fs_input.cTargetOn);
nTrials = length(tCyc);
nCells = size(npSub_tc,2);
maxCyc = max(tCyc,[],2);
data_trial = nan(prewin_frames+postwin_frames,nCells,maxCyc+1,nTrials);

tFramesOff = nan(nTrials,maxCyc);
SIx = strcmp(fs_input.trialOutcomeCell, 'success');
MIx = strcmp(fs_input.trialOutcomeCell, 'ignore');
FIx = strcmp(fs_input.trialOutcomeCell, 'failure');
nCyc = tCyc;
nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
for itrial = 1:nTrials
    if isfield(fs_input, 'tFramesOff')
        if length(fs_input.tFramesOff{itrial}>0)
            tempFramesOff = fs_input.tFramesOff{itrial};
        else
            tempFramesOff = fs_input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
            fs_input.tFramesOff{itrial} = tempFramesOff;
        end
    else
        if iscell(fs_input.nFramesOff)
            tempFramesOff = fs_input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
        else
            tempFramesOff = fs_input.nFramesOff.*(ones(1,tCyc(itrial)));
        end
    end

    tFramesOff(itrial,1:tCyc(itrial)) = tempFramesOff(1:tCyc(itrial));
    if ~isnan(cStart(itrial))
        for icyc = 1:nCyc(itrial)
            if icyc > 1
                cyc_add = ((icyc-1)*nOn)+sum(tempFramesOff(1:icyc-1));
            else
                cyc_add = 0;
            end
            if cStart(itrial)+postwin_frames-1+cyc_add <= size(npSub_tc,1)
                data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-prewin_frames+cyc_add:cStart(itrial)+postwin_frames+cyc_add-1,:);
            else
                data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
            end 
        end
    else
        data_trial(:,:,icyc,itrial) = NaN(prewin_frames+postwin_frames,nCells);
    end
end
data_f = nanmean(data_trial(1:prewin_frames,:,1,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);

targCon = celleqel2mat_padded(fs_input.tGratingContrast);
if isfield(fs_input,'doRandCon') & fs_input.doRandCon
	baseCon = nan(maxCyc,nTrials);
    for itrial = 1:nTrials
        baseCon(:,itrial) = fs_input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
    end
    ind_con = [];
else
    baseCon = celleqel2mat_padded(fs_input.tBaseGratingContrast);
    ind_con = intersect(find(targCon == 1),find(baseCon == 0));
end
baseDir = celleqel2mat_padded(fs_input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
tGratingDir = round(double(celleqel2mat_padded(fs_input.tGratingDirectionDeg)),0);
if sum(tGratingDir-baseDir) == 0
    targetDelta = tGratingDir-baseDir;
else
    targetDelta = tGratingDir;
end
deltas = unique(targetDelta);
nDelta = length(deltas);
offs = unique(tFramesOff(:,1));
noff = length(offs);
frameRateHz = fs_input.frameRateHz;

base_win =32:34;
resp_win =36:38; 

figure;
subplot(2,1,1)
plot(squeeze(nanmean(mean(data_dfof(:,:,1,:),2),4)));
vline(base_win,'k:')
vline(resp_win,'r:')
title('Baseline')
subplot(2,1,2)
sz = size(data_dfof);
data_targ = zeros(sz(1),sz(2),length([find(SIx)]));
for itrial = 1:sz(4);
    if find([find(SIx)] == itrial)
        data_targ(:,:,itrial) = data_dfof(:,:,nCyc(itrial),itrial);
    end
end
for idelta = 1:nDelta
    ind = find(targetDelta == deltas(idelta));
    plot(squeeze(nanmean(mean(data_targ(:,:,ind),2),3)));
    hold on
end
title('Target')
vline(base_win,'k:')
vline(resp_win,'r:')
legend(num2str(deltas'),'Location','northeast')

%% save FS data
fs_str = 'runs';
for irun = 1:nrun
    fs_str = [fs_str '-' expt(iexp).runs(irun,:)];
end
if ~exist(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str]))
    mkdir(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str]))
end
data_tc = npSub_tc(offset+1:end,:);
npSub_tc = npSub_tc(1:offset,:);
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str '_TCs.mat']), 'npSub_tc')
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str '_dfofData.mat']), 'data_dfof', 'prewin_frames', 'postwin_frames', 'base_win','resp_win')
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str], [expt(iexp).date '_i' expt(iexp).mouse '_' fs_str '_stimData.mat']),'prewin_frames','baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'baseCon', 'ind_con', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart', 'base_win','resp_win')

%% Direction analysis
nOn = dir_input.nScansOn;
nOff = dir_input.nScansOff;
dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
nTrials = length(dir_mat);
dir_input.trialSinceReset = nTrials;
if nTrials.*(nOn+nOff)<dir_input.counter{end}(end)
    data_tc = data_tc(1:nTrials.*(nOn+nOff),:);
end
down = 10;
nframes = size(data_tc,1)./down;
data_tc_down = squeeze(mean(reshape(data_tc, [down,nframes,nCells]),1));

tuningDownSampFactor = down;
[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuning(data_tc_down,dir_input,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));

dir_str = ['runs-' expt(iexp).dirtuning];
if ~exist(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str]))
    mkdir(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str]))
end
save(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')

        
%%
dir_mat = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(find(dir_mat>=180))-180;
oris = unique(ori_mat);
figure; 
if nCells<49
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
start = 1;
x = 0;
for ic = 1:nCells
    if start > 49
        suptitle([expt(iexp).mouse ' ' expt(iexp).date ' n = ' num2str(length(find(fitReliability<22.5))) ' well-fit'])
        print(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    title([num2str(fitReliability(1,ic)) '- r2= ' num2str(chop(R_square(1,ic),2))])
    start = start+1;
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' n = ' num2str(length(find(fitReliability<22.5))) ' well-fit'])
print(fullfile(LG_base, 'Analysis\2P', [expt(iexp).date '_i' expt(iexp).mouse], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str], [expt(iexp).date '_i' expt(iexp).mouse '_' dir_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
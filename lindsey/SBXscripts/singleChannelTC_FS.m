%% get path names
date = '180416';
ImgFolder = strvcat('002');
time = strvcat('1227');
mouse = 'i843';
doFromRef = 0;
ref = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    %CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,1000);
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
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc
%% For behavior experiments
% Plot outcome by trial number
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');

figure;
plot(smooth(SIx,10));
hold on
plot(smooth(MIx,10));

% Crop data and input struct
input = trialChopper(input,[1 200]);
data = data(:,:,input.counterValues{1}(1):input.counterValues{end}(end));

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

data_avg = mean(data(:,:,40001:40500),3);
%% Register data

if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out_bx);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    %load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')

figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% find activated cells

tCyc = cell2mat(input.tCyclesOn);
cStart = cell2mat(input.cFirstStim);
cStim = cell2mat(input.cStimOn);
cTarget = celleqel2mat_padded(input.cTargetOn);
nTrials = length(tCyc);
sz = size(data_reg);
data_f = zeros(sz(1),sz(2),nTrials);
data_base = zeros(sz(1),sz(2),nTrials);
data_base2 = zeros(sz(1),sz(2),nTrials);
data_targ = zeros(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        data_f(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
        data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+10:cStart(itrial)+20),3);
        if cStim(itrial) > cStart(itrial) 
            data_base2(:,:,itrial) = mean(data_reg(:,:,cStim(itrial)+9:cStim(itrial)+19),3);
        else
            data_base2(:,:,itrial) = nan(sz(1),sz(2));
        end
    else
        data_f(:,:,itrial) = nan(sz(1),sz(2));
        data_base(:,:,itrial) = nan(sz(1),sz(2));
        data_base2(:,:,itrial) = nan(sz(1),sz(2));
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+19 < sz(3)
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+5:cTarget(itrial)+10),3);
        else
            data_targ(:,:,itrial) = nan(sz(1),sz(2));
        end
    else
        data_targ(:,:,itrial) = nan(sz(1),sz(2));
    end
end
data_base_dfof = (data_base-data_f)./data_f;
data_base2_dfof = (data_base2-data_f)./data_f;
data_targ_dfof = (data_targ-data_f)./data_f;
targCon = celleqel2mat_padded(input.tGratingContrast);
if input.doRandCon
        baseCon = ones(size(targCon));
else
    baseCon = celleqel2mat_padded(input.tBaseGratingContrast);
end
ind_con = intersect(find(targCon == 1),find(baseCon == 0));
baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
targetDelta = round(celleqel2mat_padded(input.tGratingDirectionDeg),0);
deltas = unique(targetDelta);
nDelta = length(deltas);
data_dfof_dir = zeros(sz(1),sz(2),ndir);
data_dfof2_dir = zeros(sz(1),sz(2),ndir);
[n n2] = subplotn(ndir);
figure;
for idir = 1:ndir
    ind = setdiff(find(baseDir == dirs(idir)),ind_con);
    data_dfof_dir(:,:,idir) = nanmean(data_base_dfof(:,:,ind),3);
    data_dfof2_dir(:,:,idir) = nanmean(data_base2_dfof(:,:,ind),3);
    subplot(n,n2,idir)
    imagesc(data_dfof_dir(:,:,idir))
    title(dirs(idir))
end
if sum(~isnan(data_dfof2_dir))>1
    data_dfof_dir_all = cat(3, data_dfof_dir, data_dfof2_dir);
else
    data_dfof_dir_all = data_dfof_dir;
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
data_dfof = cat(3,data_dfof_dir_all,data_dfof_targ);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
figure;
imagesc(data_dfof_max)

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)
% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);

%% neuropil mask and subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell
%% FS cycle analysis

if iscell(input.nFramesOn)
    nOn = input.nFramesOn{1};
else
    nOn = input.nFramesOn;
end

tCyc = cell2mat(input.tCyclesOn);
cStart = celleqel2mat_padded(input.cFirstStim);
cTarget = celleqel2mat_padded(input.cTargetOn);
nTrials = length(tCyc);
nCells = size(npSub_tc,2);
maxCyc = max(tCyc,[],2);
data_trial = nan(120,nCells,maxCyc+1,nTrials);
if iscell(input.nFramesOff)
    data_notarg = nan((nOn+max(celleqel2mat_padded(input.nFramesOff),[],2))*maxCyc,nCells,nTrials);
else
    data_notarg = nan((nOn+input.nFramesOff)*maxCyc,nCells,nTrials);
end
tFramesOff = nan(nTrials,maxCyc);
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');
FIx = strcmp(input.trialOutcomeCell, 'failure');
nCyc = tCyc;
nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
for itrial = 1:nTrials
    if isfield(input, 'tFramesOff')
        if length(input.tFramesOff{itrial}>0)
            tempFramesOff = input.tFramesOff{itrial};
        else
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
            input.tFramesOff{itrial} = tempFramesOff;
        end
    else
        if iscell(input.nFramesOff)
            tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
        else
            tempFramesOff = input.nFramesOff.*(ones(1,tCyc(itrial)));
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
            if cStart(itrial)+99+cyc_add <= size(npSub_tc,1)
                data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-20+cyc_add:cStart(itrial)+99+cyc_add,:);
            else
                data_trial(:,:,icyc,itrial) = NaN(120,nCells);
            end 
            if icyc == nCyc(itrial)
                if cStart(itrial)+99+cyc_add <= size(npSub_tc,1)
                    temp  = npSub_tc(cStart(itrial)-20:cStart(itrial)+cyc_add-1,:);
                    data_notarg(1:size(temp,1),:,itrial) = temp;
                end
            end
        end
    else
        data_trial(:,:,icyc,itrial) = NaN(120,nCells);
    end
end
data_f = nanmean(data_trial(1:20,:,1,:),1);
data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
data_notarg_dfof = bsxfun(@rdivide,bsxfun(@minus,data_notarg,mean(data_notarg(1:20,:,:),1)),mean(data_notarg(1:20,:,:),1));
targCon = celleqel2mat_padded(input.tGratingContrast);
if isfield(input,'doRandCon') & input.doRandCon
	baseCon = nan(maxCyc,nTrials);
    for itrial = 1:nTrials
        baseCon(:,itrial) = input.tBaseGratingContrast{itrial}(1:tCyc(itrial));
    end
    ind_con = [];
else
    baseCon = celleqel2mat_padded(input.tBaseGratingContrast);
    ind_con = intersect(find(targCon == 1),find(baseCon == 0));
end
baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
dirs = unique(baseDir);
ndir = length(dirs);
tGratingDir = round(double(celleqel2mat_padded(input.tGratingDirectionDeg)),0);
if sum(tGratingDir-baseDir) == 0
    targetDelta = tGratingDir-baseDir;
else
    targetDelta = tGratingDir;
end
deltas = unique(targetDelta);
nDelta = length(deltas);
offs = unique(tFramesOff(:,1));
noff = length(offs);
frameRateHz = input.frameRateHz;

base_win =21:23;
resp_win =27:29; 

figure;
if nCells<25
    ii = nCells;
else
    ii = 25;
end
for i = 1:ii
    subplot(5,5,i)
if length(ind_con)>10
    plot(squeeze(nanmean(mean(data_dfof(20:50,i,2,ind_con),2),4)))
elseif noff>1
    ind = find(tFramesOff(:,1) == offs(noff));
    plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
else
    plot(squeeze(nanmean(mean(data_dfof(20:50,i,1,:),2),4)))
end
vline(base_win-19)
vline(resp_win-19)
end

figure;
subplot(2,1,1)
plot(squeeze(nanmean(mean(data_dfof(:,:,1,:),2),4)));
vline(base_win)
vline(resp_win)
title('Baseline')
subplot(2,1,2)
sz = size(data_dfof);
data_targ = zeros(sz(1),sz(2),length([find(SIx)]));
for itrial = 1:sz(4);
    if find([find(SIx)] == itrial)
        data_targ(:,:,itrial) = data_dfof(:,:,nCyc(itrial),itrial);
    end
end
plot(squeeze(nanmean(mean(data_targ,2),3)));
title('Target')
vline(base_win)
vline(resp_win)

%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof', 'data_notarg_dfof')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'baseCon', 'ind_con', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart', 'base_win','resp_win')
%% FS analysis
%find good cells- responsive to all base or at least one base direction
[x,y] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
h = zeros(ndir,nCells);
p = zeros(ndir,nCells);


if ndir > 1
    data_dfof_dir = nan(40, nCells, maxCyc, ndir);
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
        for icyc = 1:maxCyc
            ind2 = intersect(ind,find(nCyc >= icyc));
            data_dfof_dir(:,:,icyc,idir)= nanmean(nanmean(data_dfof(:,:,icyc,ind2),3),4);
        end
        for iCell = 1:nCells
            [h(idir,iCell), p(idir,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,1,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,1,ind),1)),'tail','left','alpha',0.05./(ndir-1));
        end
    end
    good_ind = unique([find(x) find(sum(h,1)>0)]);
else
    data_dfof_dir = nan(40, nCells, maxCyc);
    ind_cyc = zeros(1,maxCyc);
    for icyc = 1:maxCyc
        if icyc == 1
            %ind = intersect(find(SIx),find(nCyc>icyc));
            ind = find(tCyc>icyc);
            ind_cyc(1,icyc)= length(ind);
        else
            %ind = intersect(find(SIx), find(nCyc >= icyc));
            ind = find(tCyc >= icyc);
            ind_cyc(1,icyc)= length(ind);
        end
        data_dfof_dir(:,:,icyc)= squeeze(nanmean(data_dfof(:,:,icyc,ind),4));
    end
    good_ind = find(x);
end

%plot ori tuning of all good cells
if ndir > 1
    figure; 
    [n n2] = subplotn(length(good_ind));
    data_dir_avg = zeros(nCells,ndir);
    max_dir = NaN(nCells,1);
    for i = 1:length(good_ind)
        subplot(n,n2,i); 
        iC = good_ind(i);
        for idir = 1:ndir
            ind = find(baseDir == dirs(idir));
            data_dir_avg(iC,idir) = squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1));
            if h(idir,iC)
                errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'or');
            else
                errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'ok');
            end
            hold on
            h_ind = find(h(:,iC));
            if length(h_ind)>0
                [max_val max_ind] = max(data_dir_avg(iC,h_ind),[],2);
                max_dir(iC,:) = h_ind(max_ind);
            else
                [max_val max_ind] = max(data_dir_avg(iC,:),[],2);
                max_dir(iC,:) = max_ind;
            end
        end
        title(['Cell ' num2str(iC) '- Pref ' num2str(dirs(max_dir(iC,:)))])
        %ylim([0 max(mean(data_dfof_dir(31:35,iC,:),1),[],3)*1.5])
        %xlim([-dirs(2) dirs(end)+dirs(2)])
    end
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
elseif ~doFromRef
    max_dir = ones(nCells,1);
end
%response by cycle- avg all cells
tt = [-19:20].*frameRateHz;
data_dfof_base = zeros(40,nCells,maxCyc);
if ~doFromRef
    for i = 1:length(good_ind)
        iC = good_ind(i);
        data_dfof_base(:,iC,:) = (data_dfof_dir(:,iC,:,max_dir(iC,:)));
    end
else
    %data_dfof_base = data_dfof_dir(:,:,1:maxCyc+1);
end

if maxCyc>1
    figure;
    subplot(2,1,1)
    for icyc = 1:maxCyc
        plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,icyc),nanmean(data_dfof_base(base_win,good_ind,icyc),1)),2))
        hold on
    end
    ylim([-0.02 0.15])
    xlabel('Time (ms)')
    legend([repmat('cycle ',[maxCyc 1]) num2str([1:maxCyc]')])
    legend('Location','Northwest')
    subplot(2,1,2)
    data_dfof_base_diff = zeros(maxCyc,length(good_ind));
    for icyc = 1:maxCyc
        data_dfof_base_diff(icyc,:) = nanmean(data_dfof_base(resp_win,good_ind,icyc),1)- nanmean(data_dfof_base(base_win,good_ind,icyc),1);
        errorbar(icyc,nanmean(data_dfof_base_diff(icyc,:),2),nanstd(data_dfof_base_diff(icyc,:),[],2)./sqrt(length(good_ind)),'o')
        hold on
    end
    ylim([0 0.15])
    ylabel('dF/F')
    xlabel('cycle #')
    set(gca,'XTick',1:maxCyc)
    suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells- ' num2str(ind_cyc)]) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc.pdf']),'-dpdf','-bestfit')
end

figure
[n n2] = subplotn(length(good_ind));
if ndir>1
    for iCell = 1:length(good_ind)
        iC = good_ind(iCell);
        subplot(n,n2,iCell)
        for idir = 1:ndir
             plot(tt,squeeze(data_dfof_dir(:,iC,1,idir)))
             hold on
        end
    end
    suptitle([mouse ' ' date])
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byDir_byCell.pdf']),'-dpdf','-bestfit')
else
    for iCell = 1:length(good_ind)
        iC = good_ind(iCell);
        subplot(n,n2,iCell)
        for icyc = 1:2
            plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,iC,icyc),nanmean(data_dfof_base(base_win,iC,icyc),1)),2))
            hold on
        end
    end
    suptitle([mouse ' ' date])
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_2Cyc_byCell.pdf']),'-dpdf','-bestfit')
end

if maxCyc>1
    off_str{1} = 'base';
    for ioff = 1:noff
        off_str{ioff+1} = num2str(offs(ioff).*frameRateHz);
    end
    %response by off frames, avg all cycles all cells
    data_dfof_off = zeros(40,nCells,noff);
    for i = 1:length(good_ind)
        iC = good_ind(i);
        for ioff = 1:noff
            temp = [];
            for icyc = 1:maxCyc-1
                if ~doFromRef
                    ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
                else
                    ind = find(tFramesOff(:,icyc) == offs(ioff));
                end
                temp = [temp squeeze(data_dfof(:,iC,icyc+1,ind))];
            end
            data_dfof_off(:,iC,ioff) = nanmean(temp,2);
        end
    end
    figure;
    subplot(2,1,1)
    plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
    hold on
    for ioff = 1:noff
        plot(tt,nanmean(bsxfun(@minus,data_dfof_off(:,good_ind,ioff), nanmean(data_dfof_off(base_win,good_ind,ioff),1)),2));
    end
    legend(off_str)
    legend('Location','Northwest')
    ylim([-0.02 0.15])
    xlabel('Time (ms)')
    subplot(2,1,2)
    errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
    hold on
    data_dfof_base_off = zeros(noff,length(good_ind));
    for ioff = 1:noff
        temp = nanmean(data_dfof_off(28:31,good_ind,ioff),1)- nanmean(data_dfof_off(20:23,good_ind,ioff),1);
        errorbar((offs(ioff)*frameRateHz)./1000,mean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
        hold on
    end
    set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
    xlim([0.05 10])
    ylim([0 0.15])
    ylabel('dF/F')
    xlabel('Off Interval (ms)')
    suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt.pdf']),'-dpdf', '-bestfit')

    %response by off frames by cell, avg all cycles
    figure;
    [n n2] = subplotn(length(good_ind));
    for iCell = 1:length(good_ind)
        iC = good_ind(iCell);
        subplot(n, n2, iCell)
        plot(tt,bsxfun(@minus,data_dfof_base(:,iC,1),nanmean(data_dfof_base(base_win,iC,1),1)))
        hold on
        for ioff = 1:noff
            plot(tt,bsxfun(@minus,data_dfof_off(:,iC,ioff), nanmean(data_dfof_off(base_win,iC,ioff),1)));
        end
    end
    suptitle([mouse ' ' date]) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_byCell.pdf']),'-dpdf', '-bestfit')


    %response by off frames by cycle, avg all cells
    data_dfof_off_cyc = zeros(40,nCells,noff,maxCyc);
    ind_n = zeros(maxCyc,noff);
    for i = 1:length(good_ind)
        iC = good_ind(i);
        for ioff = 1:noff
            for icyc = 1:maxCyc
                if ~doFromRef
                    ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
                    ind_n(icyc,ioff) = length(ind);
                else
                    ind = find(tFramesOff(:,icyc) == offs(ioff));
                end
                data_dfof_off_cyc(:,iC,ioff,icyc) = nanmean(data_dfof(:,iC,icyc+1,ind),4);
            end
        end
    end
    figure;
    [n n2] = subplotn(double(maxCyc-1));
    for icyc = 1:maxCyc-1
        subplot(n,n2,double(icyc))
        plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
        hold on
        for ioff = 1:noff
            plot(tt,mean(bsxfun(@minus,data_dfof_off_cyc(:,good_ind,ioff,icyc),mean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1)),2));
        end
        if icyc == 1
            legend(off_str,'Location','Northwest')
        end
        title(['Cycle ' num2str(icyc+1) '- ' num2str(ind_n(icyc,:)) ' trials'])
        ylim([-0.02 0.15])
        xlabel('Time (ms)')
    end
    suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt.pdf']),'-dpdf', '-bestfit')

    figure;
    for icyc = 1:(maxCyc-1)
        subplot(n,n2,double(icyc))
        errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
        hold on
        for ioff = 1:noff
            temp = bsxfun(@minus, nanmean(data_dfof_off_cyc(resp_win,good_ind,ioff,icyc),1), nanmean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1));
            errorbar((offs(ioff)*frameRateHz)./1000,nanmean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
            hold on
        end
        set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
        xlim([0.05 10])
        ylim([-0.02 0.15])
        ylabel('dF/F')
        xlabel('Off Interval (s)')
        title(['Cycle ' num2str(icyc+1)])
    end
    suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt_quant.pdf']),'-dpdf', '-bestfit')
end

%find cells responsive to target or at least one target direction
data_dfof_targ = nan(40,nCells,nTrials);
tFramesOff_targ = nan(1,nTrials);
for itrial = 1:nTrials
    if ~isnan(cTarget(itrial))
        icyc = nCyc(itrial);
        data_dfof_targ(:,:,itrial) = data_dfof(:,:,icyc,itrial);
        tFramesOff_targ(1,itrial) = tFramesOff(itrial,icyc-1);
    end
end
if ~input.doPairedPulse
    [x_targ y_targ] = ttest(squeeze(nanmean(data_dfof_targ(base_win,:,:),1))',squeeze(nanmean(data_dfof_targ(resp_win,:,:),1))','tail','left','alpha',0.05);
    if ndir > 1
        h_targ = zeros(ndir,nCells);
        p_targ = zeros(ndir,nCells);
        for idir = 1:ndir
            ind1 = intersect(find(baseDir == dirs(idir)),find(targetDelta== deltas(1)));
            shift = idir-2;
            if shift<1
                shift = shift+input.baseGratingDirectionStepN;
            end
            ind2 = intersect(find(baseDir == dirs(shift)),find(targetDelta== deltas(2)));
            ind = [ind1 ind2];
            for iCell = 1:nCells
                [h_targ(idir,iCell), p_targ(idir,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof_targ(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
            end
        end
    else
        h_targ = zeros(nDelta,nCells);
        p_targ = zeros(nDelta,nCells);
        for itarg = 1:nDelta
            ind = find(targetDelta== deltas(itarg));
            for iCell = 1:nCells
                [h_targ(itarg,iCell), p_targ(itarg,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
            end
        end
    end
else
    h_targ = h;
    x_targ = x;
end
good_targ_ind = intersect(good_ind, unique([find(sum(h_targ,1)>0) find(x_targ)]));
good_resp_ind = unique([find(sum(h_targ,1)>0) find(x_targ) good_ind]);


save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']), 'data_dfof', 'max_dir','good_ind', 'good_targ_ind', 'good_resp_ind')

if ndir > 1
    %figure for each cell with base, target and target by change
    for iCell = 1:nCells
        figure;
        for idir = 1:ndir
            ind = find(baseDir == dirs(idir));
            subplot(2,2,1)
            plot(mean(data_dfof(:,iCell,1,ind),4))
            hold on
            if length(deltas)>1
                subplot(2,2,2)
                ind1 = intersect(ind,find(targetDelta== deltas(1)));
                shift = idir-2;
                if shift<1
                    shift = shift+input.baseGratingDirectionStepN;
                end
                ind2 = intersect(find(baseDir == dirs(shift)),find(targetDelta== deltas(2)));
                ind_t = [ind1 ind2];
                plot(mean(data_dfof(:,iCell,6,ind_t),4)-mean(mean(data_dfof(base_win,iCell,6,ind_t),1),4))
                hold on
            end
            for itarg = 1:length(deltas)
                ind3 = intersect(ind, find(targetDelta == deltas(itarg)));
                subplot(2,2,2+itarg)
                plot(mean(data_dfof(:,iCell,nCyc(1),ind3),4)-mean(mean(data_dfof(base_win,iCell,nCyc(1),ind3),1),4))
                hold on
            end
            subplot(2,2,1)
            title(['Base- ' num2str(dirs(find(h(:,iCell))))])
            subplot(2,2,2)
            title(['Target- ' num2str(dirs(find(h_targ(:,iCell))))])
            subplot(2,2,3)
            title(['Target- ' num2str(deltas(1)) ' deg change'])
            if length(deltas)>1
            subplot(2,2,4)
            title(['Target- ' num2str(deltas(2)) ' deg change'])
            end
        end
        base = num2str(length(find(good_ind == iCell)));
        if base == '1'
            pref = num2str(dirs(max_dir(iCell,:)));
            targ =  num2str(length(find(good_targ_ind == iCell)));
        else
            pref = 'NaN';
            targ = 'NaN';
        end
        suptitle(['Cell ' num2str(iCell) '- Pref is ' pref '- Base is ' base '- Target is ' targ])
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'AllCells',[date '_' mouse '_' run_str '_avgRespCell' num2str(iCell) '.pdf']),'-dpdf')
        close all
    end

    data_dfof_targ = zeros(40,nCells,ndir,nDelta);
    for idir = 1:ndir
        ind_dir = find(baseDir == dirs(idir));
        for itarg = 1:nDelta
            ind_targ = intersect(ind_dir, find(targetDelta == deltas(itarg)));
            if length(ind_targ)>4
                data_dfof_targ(:,:,idir,itarg) = squeeze(mean(data_dfof(:,:,end,ind_targ),4));
            else
                data_dfof_targ(:,:,idir,itarg) = NaN(40,1);
            end
        end
    end

    %average target response for absolute preferred stim
    data_dfof_delta = zeros(40,sz(2),nDelta);
    for i = 1:length(good_targ_ind)
        iC = good_targ_ind(i);
        idir = max_dir(iC,:);
        shift1 = idir-1;
        if shift1<1
            shift1 = shift1+ndir;
        end
        ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));    
        data_dfof_delta(:,iC,1) = mean(data_dfof(:,iC,6,ind1),4);
        shift2 = idir-3;
        if shift2<1
            shift2 = shift2+ndir;
        end
        ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));    
        data_dfof_delta(:,iC,2) = mean(data_dfof(:,iC,6,ind2),4);
    end
    figure;
    subplot(2,1,1)
    plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_targ_ind,1),mean(data_dfof_base(base_win,good_targ_ind,1),1)),2))
    hold on
    for itarg = 1:nDelta
        plot(tt, mean(bsxfun(@minus,data_dfof_delta(:,good_targ_ind,itarg),mean(data_dfof_delta(base_win,good_targ_ind,itarg),1)),2))
        hold on
    end
    legend('base', [num2str(deltas(1)) ' deg target'], [num2str(deltas(2)) ' deg target'])
    legend('Location','Northwest')
    xlabel('Time (ms)')
    ylabel('dF/F')
    subplot(2,1,2)
    temp = mean(data_dfof_base(resp_win,good_targ_ind,1),1)-mean(data_dfof_base(base_win,good_targ_ind,1),1);
    errorbar(0,mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
    hold on
    for itarg = 1:nDelta
        temp  = mean(data_dfof_delta(resp_win,good_targ_ind,itarg),1)-mean(data_dfof_delta(base_win,good_targ_ind,itarg),1);
        errorbar(deltas(itarg), mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
        hold on
    end
    xlabel('Degree change')
    ylabel('dF/F')
    ylim([0 0.1])
    xlim([-10 100])
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta.pdf']),'-dpdf')

    %tuning curves for baseline
    data_dfof_base_tuning = zeros(40,sz(2),ndir,nCyc);
    for i = 1:length(good_ind)
        iC = good_ind(i);
        for idir = 1:ndir
            dir_temp = max_dir(iC,:)+(idir-1);
            if dir_temp>ndir
                dir_temp = rem(dir_temp,ndir);
            end
            ind = find(baseDir == dirs(dir_temp));
            for icyc = 1:nCyc
                data_dfof_base_tuning(:,iC,idir,icyc) = mean(data_dfof(:,iC,icyc,ind),4);
            end
        end
    end
    %for base responsive cells
    figure;
    for icyc = 1:nCyc
        for idir = 1:ndir
            subplot(2,double(nCyc),double(icyc))
            plot(tt,mean(bsxfun(@minus,data_dfof_base_tuning(:,good_ind,idir,icyc),mean(data_dfof_base_tuning(base_win,good_ind,idir,icyc),1)),2))
            hold on
            subplot(2,double(nCyc),double(icyc+nCyc))
            temp = bsxfun(@minus,mean(data_dfof_base_tuning(resp_win,good_ind,idir,icyc),1),mean(data_dfof_base_tuning(base_win,good_ind,idir,icyc),1));
            errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_ind)),'o');
            hold on
        end
        ylim([0 0.1])
        xlim([-30 180])
        xlabel('Ori (deg)')
        subplot(2,double(nCyc),double(icyc))
        ylim([-0.02 0.1])
        title(['Cycle # ' num2str(icyc)])
        xlabel('Time (ms)')
    end
    suptitle('Baseline responsive cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_BaseOriTuning_goodCells.pdf']),'-dpdf')

    %for target responsive cells
    figure;
    for icyc = 1:nCyc
        for idir = 1:ndir
            subplot(2,double(nCyc),double(icyc))
            plot(tt,mean(bsxfun(@minus,data_dfof_base_tuning(:,good_targ_ind,idir,icyc),mean(data_dfof_base_tuning(base_win,good_targ_ind,idir,icyc),1)),2))
            hold on
            subplot(2,double(nCyc),double(icyc+nCyc))
            temp = bsxfun(@minus,mean(data_dfof_base_tuning(resp_win,good_targ_ind,idir,icyc),1),mean(data_dfof_base_tuning(base_win,good_targ_ind,idir,icyc),1));
            errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
            hold on
        end
        ylim([0 0.1])
        xlim([-30 180])
        xlabel('Ori (deg)')
        subplot(2,double(nCyc),double(icyc))
        ylim([-0.02 0.1])
        title(['Cycle # ' num2str(icyc)])
        xlabel('Time (ms)')
    end
    suptitle('Target responsive cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_BaseOriTuning_goodTargCells.pdf']),'-dpdf')

    %tuning curves for targets
    data_dfof_targ_tuning = zeros(40,sz(2),ndir,nDelta+4);
    for i = 1:length(good_targ_ind)
        iC = good_targ_ind(i);
        for idir = 1:ndir
            dir_temp = max_dir(iC,:)+(idir-1);
            if dir_temp>ndir
                dir_temp = rem(dir_temp,ndir);
            end
            ind = find(baseDir == dirs(dir_temp));
            data_dfof_targ_tuning(:,iC,idir,1) = mean(data_dfof(:,iC,1,ind),4);
            data_dfof_targ_tuning(:,iC,idir,2) = mean(mean(data_dfof(:,iC,3:5,ind),3),4);
            shift1 = max_dir(iC,:)-1+(idir-1);
            if shift1<1
                shift1 = shift1+ndir;
            elseif shift1>ndir
                shift1 = rem(shift1,ndir);
            end
            ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));    
            data_dfof_targ_tuning(:,iC,idir,3) = mean(data_dfof(:,iC,6,ind1),4);
            shift2 = max_dir(iC,:)-1+(idir-3);
            if shift2<1
                shift2 = shift2+ndir;
            elseif shift2>ndir
                shift2 = rem(shift2,ndir);
            end
            ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));    
            data_dfof_targ_tuning(:,iC,idir,4) = mean(data_dfof(:,iC,6,ind2),4);
        end
    end
    figure;
    for i = 1:4
        for idir = 1:ndir
            subplot(2,4,i)
            plot(tt,mean(bsxfun(@minus,data_dfof_targ_tuning(:,good_targ_ind,idir,i),mean(data_dfof_targ_tuning(base_win,good_targ_ind,idir,i),1)),2))
            hold on
            subplot(2,4,i+4)
            temp = bsxfun(@minus,mean(data_dfof_targ_tuning(resp_win,good_targ_ind,idir,i),1),mean(data_dfof_targ_tuning(base_win,good_targ_ind,idir,i),1));
            errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
            hold on
        end
        ylim([0 0.1])
        ylabel('dF/F')
        xlim([-30 180])
        xlabel('Ori (deg)')
        subplot(2,4,i)
        ylim([-0.02 0.1])
        ylabel('dF/F')
        xlabel('Time (ms)')
    end
    subplot(2,4,1)
    title('Base 1')
    subplot(2,4,2)
    title('Base 3-5')
    subplot(2,4,3)
    title(['Target- ' num2str(deltas(1)) ' deg change'])
    subplot(2,4,4)
    title(['Target- ' num2str(deltas(2)) ' deg change'])
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_OriTuning.pdf']),'-dpdf')

    %tuning curves for targets by interval
    data_dfof_targ_tuning_off = zeros(40,sz(2),ndir,nDelta,noff);
    for i = 1:length(good_targ_ind)
        iC = good_targ_ind(i);
        for idir = 1:ndir
            shift1 = max_dir(iC,:)-1+(idir-1);
            if shift1<1
                shift1 = shift1+ndir;
            elseif shift1>ndir
                shift1 = rem(shift1,ndir);
            end
            shift2 = max_dir(iC,:)-1+(idir-3);
            if shift2<1
                shift2 = shift2+ndir;
            elseif shift2>ndir
                shift2 = rem(shift2,ndir);
            end
            for ioff = 1:noff
                ind1 = intersect(find(tFramesOff(:,5) == offs(ioff)),intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1))));    
                data_dfof_targ_tuning_off(:,iC,idir,1,ioff) = mean(data_dfof(:,iC,6,ind1),4);
                ind2 = intersect(find(tFramesOff(:,5) == offs(ioff)),intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2))));    
                data_dfof_targ_tuning_off(:,iC,idir,2,ioff) = mean(data_dfof(:,iC,6,ind2),4);
            end
        end
    end

    figure; 
    for itarg = 1:nDelta
        for ioff = 1:noff
            subplot(2,3,ioff +((itarg-1)*noff))
            for idir = 1:ndir
                temp = bsxfun(@minus,mean(data_dfof_targ_tuning_off(resp_win,good_targ_ind,idir,itarg,ioff),1),mean(data_dfof_targ_tuning_off(base_win,good_targ_ind,idir,itarg,ioff),1));
                errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
                hold on
            end
            xlim([-30 180])
            xlabel('Ori (deg)')
            ylim([-0.04 0.12])
            ylabel('dF/F')
            hline(0, '--k')
            title([num2str(offs(ioff)*frameRateHz) ' ms; ' num2str(deltas(itarg)) ' deg change'])
        end
    end
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_OriTuning_byInt.pdf']),'-dpdf')

    %target response (for absolute pref stim) by interval
    data_dfof_delta_off = zeros(40,sz(2),nDelta,noff);
    for i = 1:length(good_targ_ind)
        iC = good_targ_ind(i);
        idir = max_dir(iC,:);
        shift1 = idir-1;
        if shift1<1
            shift1 = shift1+ndir;
        end
        ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));
        for ioff = 1:noff
            ind1_off = intersect(ind1, find(tFramesOff(:,5) == offs(ioff)));
            if length(ind1_off)>4
                data_dfof_delta_off(:,iC,1,ioff) = mean(data_dfof(:,iC,6,ind1_off),4);
            else
                data_dfof_delta_off(:,iC,1,ioff) = NaN(40,1);
            end
        end
        shift2 = idir-3;
        if shift2<1
            shift2 = shift2+ndir;
        end
        ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));
        for ioff = 1:noff
            ind2_off = intersect(ind2, find(tFramesOff(:,5) == offs(ioff)));
            if length(ind2_off)>4
                data_dfof_delta_off(:,iC,2,ioff) = mean(data_dfof(:,iC,6,ind2_off),4);
            else
                data_dfof_delta_off(:,iC,2,ioff) = NaN(40,1);
            end
        end
    end
    figure;
    for itarg = 1:nDelta
        subplot(2,2,itarg)
        plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_targ_ind,1),mean(data_dfof_base(base_win,good_targ_ind,1),1)),2))
        hold on
        for ioff = 1:noff
            plot(tt, nanmean(bsxfun(@minus, data_dfof_delta_off(:,good_targ_ind,itarg,ioff),nanmean(data_dfof_delta_off(base_win,good_targ_ind,itarg,ioff),1)),2))
            hold on
        end
        xlabel('Time (ms)')
        title([num2str(deltas(itarg)) ' deg target'])
        subplot(2,2,itarg+2)
        temp = mean(data_dfof_base(resp_win,good_targ_ind,1),1)-mean(data_dfof_base(base_win,good_targ_ind,1),1);
        errorbar(8,mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
        hold on
        for ioff = 1:noff
            temp  = bsxfun(@minus, nanmean(data_dfof_delta_off(resp_win,good_targ_ind,itarg,ioff),1),nanmean(data_dfof_delta_off(base_win,good_targ_ind,itarg,ioff),1));
            errorbar((offs(ioff)*frameRateHz)./1000, nanmean(temp,2),nanstd(temp,[],2)./sqrt(length(~isnan(temp))),'o')
            hold on
        end
        set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
        xlim([0.05 10])
        ylim([-0.02 0.1])
        ylabel('dF/F')
        xlabel('Off Interval (s)')
    end
    subplot(2,2,2)
    legend('base', num2str(offs(1).*frameRateHz),num2str(offs(2).*frameRateHz),num2str(offs(3).*frameRateHz))
    legend('Location','Northwest')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byInt.pdf']),'-dpdf')
else
    figure;
    [n n2] = subplotn(length(good_resp_ind));
    for iCell = 1:length(good_resp_ind)
        iC = good_resp_ind(iCell);
        subplot(n,n2,iCell)
        plot(tt,data_dfof_base(:,iC,1))
        hold on
        for itarg = 1:nDelta
            ind = find(targetDelta == deltas(itarg));
            plot(tt,nanmean(data_dfof_targ(:,iC,ind),3)-nanmean(nanmean(data_dfof_targ(base_win,iC,ind),1),3));
            hold on;
        end
    end
    suptitle([mouse ' ' date]) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byCell.pdf']),'-dpdf','-bestfit')

    targ_only_ind = unique([find(sum(h_targ,1)>0) find(x_targ)]);
    indn = zeros(nDelta, noff);
    figure;
    for itarg = 1:nDelta
        for ioff = 1:noff
            ind = intersect(find(targetDelta == deltas(itarg)), find(tFramesOff_targ== offs(ioff)));
            indn(itarg,ioff) = length(ind);
            subplot(2,nDelta,itarg)
            plot(tt, nanmean(nanmean(bsxfun(@minus,data_dfof_targ(:,targ_only_ind,ind), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3),2))
            hold on
            subplot(2,nDelta,itarg+nDelta)
            temp = nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,targ_only_ind,ind),1), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3);
            errorbar(offs(ioff)*frameRateHz, nanmean(temp,2), nanstd(temp,[],2)./sqrt(length(targ_only_ind)),'o')
            hold on
        end
        ylim([0 0.2])
        ylabel('dF/F')
        xlabel('Off interval (ms)')
        subplot(2,nDelta,itarg)
        ylim([-0.05 0.2])
        ylabel('dF/F')
        xlabel('Time (ms)')
        title([num2str(deltas(itarg)) ' deg change- ' num2str(indn(itarg,:))])
    end
    suptitle([mouse ' ' date '- ' num2str(length(targ_only_ind)) ' cells']) 
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byInt_targInd.pdf']),'-dpdf','-bestfit')

    figure; 
    [n n2] = subplotn(length(good_resp_ind));
    dirs = [0:30:150];
    for iCell = 1:length(good_resp_ind)
        iC = good_resp_ind(iCell);
        subplot(n, n2, iCell)
        errorbar(dirs(1),mean(bsxfun(@minus,nanmean(data_dfof(resp_win,iC,1,:),1),nanmean(data_dfof(base_win,iC,1,:),1)),4), nanstd(bsxfun(@minus,nanmean(data_dfof(resp_win,iC,1,:),1),nanmean(data_dfof(base_win,iC,1,:),1)),[],4)./sqrt(nTrials), 'o')
        hold on
        %errorbar(dirs(1),nanmean(bsxfun(@minus,mean(data_dfof_base(resp_win,iC,5,:),1),mean(data_dfof(base_win,iC,5,:),1)),4), std(bsxfun(@minus,mean(data_dfof(resp_win,iC,5,:),1),mean(data_dfof(base_win,iC,5,:),1)),[],4)./sqrt(nTrials), 'o')
        for itarg = 1:nDelta
            ind = find(targetDelta == deltas(itarg));
            errorbar(deltas(itarg),nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),3), nanstd(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),[],3)./sqrt(length(ind)), 'o')
        end
        xlim([-30 100])
        if length(unique(max_dir)) > 1
            title(num2str(dirs(max_dir(iC,:))))
        end
    end

    baseV90 = zeros(2,nCells);
    ind = find(targetDelta == deltas(itarg));
    for iCell = 1:length(good_targ_ind)
        iC = good_targ_ind(iCell);
        baseV90(:,iC) = [mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,1,:),1),mean(data_dfof(base_win,iC,1,:),1)),4) mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,6,ind),1),mean(data_dfof(base_win,iC,6,ind),1)),4)];
    end


    for i = 1:length(good_targ_ind)
        figure;
        iC = good_targ_ind(i);
        for ioff = 1:noff
            subplot(2, noff, ioff)
            ind = find(tFramesOff(:,5) == offs(ioff));
            plot(tt,mean(mean(bsxfun(@minus,data_dfof(:,iC,5,ind),mean(data_dfof(base_win,iC,5,ind),1)),4),2))
            hold on
            subplot(2, noff, ioff+noff)
            errorbar(dirs(1),mean(mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,5,ind),1),mean(data_dfof(base_win,iC,5,ind),1)),4),2), std(mean(bsxfun(@minus,mean(data_dfof(resp_win,:,5,:),1),mean(data_dfof(base_win,:,5,:),1)),4),[],2)./sqrt(nCells), 'o')
            hold on
            for itarg = 1:nDelta
                subplot(2, noff, ioff)
                ind = intersect(find(tFramesOff(:,5) == offs(ioff)), find(targetDelta == deltas(itarg)));
                plot(tt,mean(mean(bsxfun(@minus,data_dfof(:,iC,6,ind),mean(data_dfof(base_win,iC,6,ind),1)),4),2))
                hold on
                subplot(2, noff, ioff+noff)
                errorbar(deltas(itarg),mean(mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,6,ind),1),mean(data_dfof(base_win,iC,6,ind),1)),4),2), std(mean(bsxfun(@minus,mean(data_dfof(resp_win,:,6,ind),1),mean(data_dfof(base_win,:,6,ind),1)),4),[],2)./sqrt(nCells), 'o')
                hold on
            end
            ylim([-0.02 0.3])
            xlim([-30 100])
            subplot(2, noff, ioff)
            title([num2str(offs(ioff)*frameRateHz) ' ms interval'])
            ylim([-0.02 0.3])
        end
        suptitle(['Cell #' num2str(iC) '; Pref- ' num2str(dirs(max_dir(iC,:))) ' deg'])
    end

    figure;
    delt_col = strvcat('k', 'b');
    [n n2] = subplotn(length(good_resp_ind));
    for iCell = 1:length(good_resp_ind)
        iC = good_resp_ind(iCell);
        subplot(n,n2,iCell)
        for ioff = 1:noff
            for idelta = 1:nDelta
                ind = intersect(find(tFramesOff_targ == offs(ioff)), find(targetDelta == deltas(idelta)));
                errorbar(offs(ioff).*frameRateHz,nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),3),nanstd(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),[],3), ['o' delt_col(idelta)])
                hold on
                ylim([-0.05 0.2])
            end
        end
    end
end

ind_n = zeros(ndir,nDelta,noff);
for idir = 1:ndir
    ind_dir = find(baseDir == dirs(idir));
    for itarg = 1:nDelta
        ind_targ = intersect(ind_dir, find(targetDelta == deltas(itarg)));
        for ioff = 1:noff
            ind_off = intersect(ind_targ,find(tFramesOff_targ == offs(ioff)));
            ind_n(idir,itarg,ioff) = length(ind_off);
        end
    end
end
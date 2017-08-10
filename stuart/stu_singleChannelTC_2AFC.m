%% get path names
date = '170217';%change to date you want
ImgFolder = strvcat('001');%change to which runs to pull from
time = strvcat('1320');
mouse = 'i553';%change to which mouse
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

[s,tUsername] = dos('ECHO %USERNAME%');
if tUsername(1:4) == 'ryan'
    tDir = 'ryan';
elseif tUsername(1:4) == 'lind'
    tDir = 'lindsey';
elseif strcmp(tUsername(1:4),'stua')
    tDir = 'stuart';  %used for file path
end

%% load and register
tic %~11min
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = ['\\CRASH.dhe.duke.edu\data\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)]; %2P data/image stack path
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW14\two-photon imaging\' date '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes); %loads frames
    
    if size(time,1) >= irun
        fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat']; %behavior data path
        load(fName);
        temp(irun) = input;
        if irun>1
            if isfield(input, 'tLeftTrial')
                for itrial = 1:ntrials
                    temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset; %counter of frame for trial start
                    temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset; %counter of frame for stim on
                    temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset; %counter of frame for decision
                end
            end
        end
    end
    offset = offset+nframes;
    
    data_temp = squeeze(data_temp); %get rid of dimensions with one level
    data = cat(3,data,data_temp);
end
input = concatenateDataBlocks(temp); %combine mat files
clear data_temp
clear temp %clear important for memory
toc/60
%% Choose register interval
nep = floor(size(data,3)./10000); %how many 10,000 frame chunks
[n,n2] = subplotn(nep);
figure
for i = 1:nep;
    %subplot(n,n2,i); 
    figure
    imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); %average of data along x and y through certain time chunks
    title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]);
end

figure
imagesc(mean(data(:,:,60001:60100),3))


%% Register data
tic %~2hr for first time, ~45 second time round
data_avg = mean(data(:,:,60001:60100),3); %put time chunk chosen from last section
%if already registered take from file
if exist(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [~, data_reg]=stackRegister_MA(data,[],[],out);
    clear out
    1
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
    %%%save 'data_avg' or 'data_reg'?
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    2
else %register the data, then save it
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    3
end
clear data %important for memory
toc/60
%% test stability

figure
for i = 1:nep;
    subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3));
    title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]);
end
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
%print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']))
%save as pdf, but gives error in adobe
%% find activated cells

%try different strategies

tLeftTrial = celleqel2mat_padded(input.tLeftTrial); %binary of if left trial in vector form
%tGratingContrast = celleqel2mat_padded(input.tGratingContrast);%vector of contrast
cStimOn = celleqel2mat_padded(input.cStimOn);%vector of frame start times
cDecision = celleqel2mat_padded(input.cDecision);%vector of frame decision times
SIx = strcmp(input.trialOutcomeCell, 'success');%trials that were correct (success index)
nTrials = length(tLeftTrial);
sz = size(data_reg); %y, x, z(time)

data_f = cell(9,1);

data_f{1} = nan(sz(1),sz(2),nTrials); %baseline for before all trials, stim, dim( y, x, trials)
data_f{2} = nan(sz(1),sz(2),1); %baseline for before any trials, beware shifts throughout recording
%data_f{3} = nan(sz(1),sz(2),nTrials); %baseline for before all trials, choice
%data_f{4} = data_f{3}; %baseline for prestim and choice

%data_f{5} = nan(1,1,nTrials); %baseline of whole frame for each pre choice
%data_f{6} = data_f{5}; %baseline of whole frame for each prestim
%data_f{7} = data_f{5}; %for before choice and stims
%data_f{8} = nan(1,1,sz(3)); %for each frame
%data_f{9} = nan(1); %for before anything

data_f{3} = nan(sz(1),sz(2),1);

data_targ = nan(sz(1),sz(2),nTrials); %
data_targ_late = nan(sz(1),sz(2),nTrials); %
data_resp = nan(sz(1),sz(2),nTrials); %
data_acells = nan(sz(1),sz(2),1);

data_f{2}(:,:,1) = mean(data_reg(:,:,1:cStimOn(1)-1),3); %average before first stim, shows some activity, can pick out cells
data_f{3}(:,:,1) = mean(data_reg(:,:,:),3);
%data_f{9}(1) = mean(mean(mean(data_reg(:,:,1:cStimOn(1)-1))));
for itrial = 1:nTrials %assign average fluourescent values
    data_f{1}(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3); %average of 20 frames before stim on
    %data_f{3}(:,:,itrial) = mean(data_reg(:,:,cDecision(itrial)-20:cDecision(itrial)-1),3);
    %data_f{4}(:,:,itrial) = mean(cat(3,data_f{1}(:,:,itrial),data_f{3}(:,:,itrial)),3); %average of f1 and f3
    
    %data_f{6}(itrial) = mean(mean(mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1))));
    %data_f{5}(itrial) = mean(mean(mean(data_reg(:,:,cDecision(itrial)-20:cDecision(itrial)-1))));
    %data_f{7}(itrial) = mean([data_f{6}(itrial),data_f{5}(itrial)]); %average of f6 and f5
    
    data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3); %average of 5 to 25 frames after stim on
    data_targ_late(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5+ceil(input.stimOnTimeMs./frame_rate):cStimOn(itrial)+25+ceil(input.stimOnTimeMs./frame_rate)),3); %5 to 25 frames after stim off
    if cDecision(itrial)+25<nframes
        data_resp(:,:,itrial) = mean(data_reg(:,:,cDecision(itrial)+5:cDecision(itrial)+25),3);%average after decision frame
    end
end

for ii=1:sz(3)/100
    data_acells=max(cat(3, mean(data_reg(:,:,ii*100-99:ii*100),3), data_acells),[],3); 
    
end
data_acells_dfof = (data_acells - data_f{3})./data_f{3};
figure;imagesq(data_acells_dfof)

%for ii=1:sz(3)
%    data_f{8}(ii)=mean(mean(mean(data_reg(:,:,ii))));
%end

%calculate variance for baseline too, not just mean

%change in fluourescence over baseline fluourescence (deltaF/F_b)(dfof)
data_targ_dfof=cell(9,1);
data_resp_dfof=cell(9,1);
data_targ_late_dfof=cell(9,1);
for ii=1:3%9
    if length(size(data_f{ii}))<3
        temp=nan(sz(1),sz(2),nTrials);
        for jj=1:nTrials
            temp(:,:,jj)=data_f{ii}(:,:,1);
        end
    elseif size(data_f{ii},1)==1
        if ii==9
            temp=data_f{ii}*ones(sz(1),sz(2),nTrials);
        elseif ii==8
            %continue
            for jj=1:nTrials
                temp(:,:,jj,1)=mean(data_f{ii}(cStimOn(itrial)+5:cStimOn(itrial)+25))*ones(sz(1),sz(2));
                temp(:,:,jj,2)=mean(data_f{ii}(cStimOn(itrial)+5+ceil(input.stimOnTimeMs./frame_rate):cStimOn(itrial)+25+ceil(input.stimOnTimeMs./frame_rate)))*ones(sz(1),sz(2));
                temp(:,:,jj,3)=mean(data_f{ii}(cDecision(itrial)+5:cDecision(itrial)+25))*ones(sz(1),sz(2));
            end
        else
            for jj=1:nTrials
                temp(:,:,jj)=data_f{ii}(jj)*ones(sz(1),sz(2),1);
            end
        end
    else
        temp=data_f{ii};
    end
    if ii==8
        data_targ_dfof{ii} = (data_targ-temp(:,:,:,1))./temp(:,:,:,1);
        data_targ_late_dfof{ii} = (data_targ_late-temp(:,:,:,2))./temp(:,:,:,2);
        data_resp_dfof{ii} = (data_resp-temp(:,:,:,3))./temp(:,:,:,3);
    else
        data_targ_dfof{ii} = (data_targ-temp)./temp;%need way to compute across dimensions
        data_resp_dfof{ii} = (data_resp-temp)./temp;
        data_targ_late_dfof{ii} = (data_targ_late-temp)./temp;
    end
end

data_dfof_L=cell(9,1);
data_dfof_late_L=cell(9,1);
data_dfof_R=cell(9,1);
data_dfof_late_R=cell(9,1);
data_dfof_resp=cell(9,1);
data_dfof=cell(9,1);
data_dfof_max=cell(9,1);

indL = find(tLeftTrial); %index of left trials
indR = find(~tLeftTrial);%right
for ii=1:3%9
    %dfof for only left trials
    data_dfof_L{ii} = mean(data_targ_dfof{ii}(:,:,indL),3);
    data_dfof_late_L{ii} = mean(data_targ_late_dfof{ii}(:,:,indL),3);
    %dfof for only right trials
    %indR = intersect(find(tGratingContrast==1),find(~tLeftTrial)); %and also only high contrast
    data_dfof_R{ii} = mean(data_targ_dfof{ii}(:,:,indR),3);
    data_dfof_late_R{ii} = mean(data_targ_late_dfof{ii}(:,:,indR),3);
    data_dfof_resp{ii} = nanmean(data_resp_dfof{ii}(:,:,SIx),3);
    figure;
    subplot(2,3,1)
    imagesc(data_dfof_L{ii})
    title('Left Stim')
    subplot(2,3,2)
    imagesc(data_dfof_R{ii})
    title('Right Stim')
    subplot(2,3,3)
    imagesc(data_dfof_late_L{ii})
    title('Left Stim Late')
    subplot(2,3,4)
    imagesc(data_dfof_late_R{ii})
    title('Right Stim Late')
    subplot(2,3,5)
    imagesc(data_dfof_resp{ii})
    title('Decision')
    %print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '__' num2str(ii) '_FOV_dFoF.pdf']), '-dpdf')
    data_dfof{ii} = cat(3, data_dfof_resp{ii}, data_dfof_late_L{ii}, data_dfof_late_R{ii}, data_dfof_L{ii}, data_dfof_R{ii}); %all the responses
    data_dfof_max{ii} = max(data_dfof{ii},[],3);
    figure; imagesc(data_dfof_max{ii})
end

%only f1 and f2 look promising, second f3 looks best
data_dfof{4}=data_acells_dfof; %all cells not seperated by stim type
data_dfof_max{4}=data_acells_dfof;
%% cell segmentation, semimanual
mask_cell=cell(4,1);
mask_np=cell(4,1);
for ii=1:1%4
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof{ii};
%data_acells_dfof
for iStim = 1:size(data_dfof{ii},3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(mask_all >= 1) = 0;
    bwout = imCellEditInteractive(mask_data_temp);%click on cells, (image toolbox to automate?)
    mask_2 = bwlabel(bwout);
    mask_all = mask_all+mask_2;
    %close all
end
mask_cell{ii} = bwlabel(mask_all);

% bwout = imCellEditInteractive(data_dfof_max);
% mask_cell = bwlabel(bwout);

mask_np{ii} = imCellNeuropil(mask_cell{ii}, 3, 5);
end
%save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
%number 4 definitely best
clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir
%% look at mask options
figure
for ii=1:length(mask_cell)
    subplot(2,2,ii)
    imagesq(mask_cell{ii})
end
%fourth looks strongest
%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell{1}); %Calculate pixel intensities inside multiple ROIs
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell{1});%Downsample stack array by ratio then calculate pixel intensities, reduces noise
nCells = size(data_tc,2);
%np_tc = stackGetTimeCourses(data_reg,mask_np);
clear np_tc np_tc_down
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);

for i = 1:nCells
    np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np{1}(:,:,i));
    np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np{1}(:,:,i));
    fprintf(['Cell #' num2str(i) '%s/n'])
end
%get weights by maximizing skew, different for each cell, you want signal
%to be sparse and positive
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew, ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);%Binary Singleton Expansion Function

%%
%save time courses
save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_tc np_tc
clear data_reg data_reg_down

%% %% 2AFC analysis

tLeftTrial = celleqel2mat_padded(input.tLeftTrial); %vector of if left trial or not
nside = 2;
tGratingContrast = celleqel2mat_padded(input.tGratingContrast);%vector of contrasts
cons = unique(tGratingContrast);%possible contrast values
ncon = length(cons);
cStimOn = celleqel2mat_padded(input.cStimOn);%stim vector
cDecision = celleqel2mat_padded(input.cDecision);%decision vector
nTrials = length(tLeftTrial);
nCells = size(npSub_tc,2);
nframes = size(npSub_tc,1);
data_stim = nan(50,nCells,nTrials);
data_choice = nan(50,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes %if within range of recording:
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:); %collect 20 frames before and 29 after
    end
    if cDecision(itrial) %if there was a decision:
        if cDecision(itrial)+29 < nframes
            data_choice(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:); %collect 20 frames before and 29 after
        end
    end
end
%think about how you want to get baseline for this
dataf{1} = mean(data_stim(1:5,:,:),1);
dataf{2} = mean(data_stim(1:10,:,:),1); %average of first 20 (or 10) rows across x and cells
dataf{3} = mean(data_stim(1:20,:,:),1);
g=3;

%for ii=1:3
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf{g}), dataf{g}); % Binary Singleton Expansion Function, subtraction then division, dF/F
data_choice_dfof = bsxfun(@rdivide, bsxfun(@minus, data_choice, dataf{g}), dataf{g});
%end

ind = intersect(find(~tLeftTrial),find(tGratingContrast ==1));%index of right high contrast trials
figure; plot(mean(mean(data_stim_dfof(:,:,ind),2),3)); %data_stim_dfof(frames,cells,trial), graph of F collapsed over cells and trials
hold on; scatter(21,0)%frame of when stim or choice happens
figure; plot(mean(mean(data_choice_dfof(:,:,ind),2),3)); %data_stim_dfof(frames,cells,trial), graph of F collapsed over cells and trials
hold on; scatter(21,0)%frame of when stim or choice happens

SIx = strcmp(input.trialOutcomeCell, 'success');%success index
MIx = strcmp(input.trialOutcomeCell, 'ignore');%miss index
FIx = strcmp(input.trialOutcomeCell, 'incorrect');%fail index

h_stim = zeros(nside, nCells);
p_stim = zeros(nside, nCells);
h_choice = zeros(nside, nCells);
p_choice = zeros(nside, nCells);

for iS = 1:nside
    indS = intersect(find(SIx), find(tLeftTrial == iS-1));%success and right, then left
    indC = intersect(indS, find(tGratingContrast == 1));%plus high contrast
    [h_stim(iS,:), p_stim(iS,:)] = ttest(squeeze(mean(data_stim_dfof(36:50,:,indC),1))', squeeze(mean(data_stim_dfof(1:15,:,indC),1))', 'tail', 'right');% ttest for before and after stim for  right(row1) and left(row2)
    [h_choice(iS,:), p_choice(iS,:)] = ttest(squeeze(mean(data_choice_dfof(36:50,:,indC),1))', squeeze(mean(data_choice_dfof(1:15,:,indC),1))', 'tail', 'right'); %or choice
end
good_ind_stim = find(sum(h_stim,1)>0);%which cells, left or right trials, that reject null (5%)
good_ind_choice = find(sum(h_choice,1)>0);

% %
tt = (-20:29)*1000/frame_rate;%frames into msec for the 50 frames
figure;
[n, n2] = subplotn(nCells);%gets minimum dimensions for plotting
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        indS = intersect(find(tGratingContrast == 1), intersect(find(SIx), find(tLeftTrial == iS-1)));%high contrast and succesful and right, then left
        plot(tt,nanmean(data_stim_dfof(:,iCell,indS),3));%each cell's dF/F across time for all selected trials
        xlabel('msec')
        ylabel('dF/F')
        hold on
        if find(good_ind_stim == iCell)
            good_str = ' resp';%responsive based off ttest h value
        else
            good_str = ' not resp';
        end
        title(['Cell # ' num2str(iCell) ' is' good_str])
    end
end
suptitle([mouse ' ' date '- stimAlign- blue is right; red is left'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

cmap = [blues(ncon); reds(ncon)];%color map
%tt = (-19:30)*1000/frame_rate;%redundant
figure;
[n, n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        for icon = 1:ncon
            indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));%as above but with variable contrast
            plot(tt,nanmean(data_stim_dfof(:,iCell,indS),3),'Color',cmap(icon+((iS-1)*ncon),:));
            xlabel('msec')
            ylabel('dF/F')
            hold on;
            if find(good_ind_stim == iCell)
                good_str = ' resp';
            else
                good_str = ' not resp';
            end
            title(['Cell # ' num2str(iCell) ' is' good_str])
        end
    end
end
suptitle([mouse ' ' date '- stimAlign- by side and contrast- left is red; right is blue'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide&Contrast_byCell.pdf']), '-dpdf', '-bestfit')

%plot all trials for cells' response to different sides and contrasts
for iCell = 1:nCells
    figure
    for iS = 1:nside
        for icon = 1:ncon
            subplot(2,3,icon+((iS-1)*ncon))
            indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));
            plot(tt, squeeze(data_stim_dfof(:,iCell,indS)));
            xlabel('msec')
            ylabel('dF/F')
            title(['side ',num2str(iS),' contrast ',num2str(cons(icon))])
            ylim([-1,3])
        end
    end
    suptitle(['Cell #' num2str(iCell)])
end


figure
[n, n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        indS = intersect(find(SIx), find(tLeftTrial == iS-1));%successful and right then left
        plot(tt,nanmean(data_choice_dfof(:,iCell,indS),3));%plots response to choice, includes all contrasts
        xlabel('msec')
        ylabel('dF/F')
        hold on;
        if find(good_ind_choice == iCell)
            good_str = ' resp';
        else
            good_str = ' not resp';
        end
        title(['Cell # ' num2str(iCell) ' is' good_str])
    end
end
suptitle([mouse ' ' date '- choiceAlign- blue is right; red is left'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

figure
[n, n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        for icon = 1:ncon
            indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));
            plot(tt,nanmean(data_choice_dfof(:,iCell,indS),3),'Color',cmap(icon+((iS-1)*ncon),:));%with all contrasts seperate
            xlabel('msec')
            ylabel('dF/F')
            hold on;
            if find(good_ind_choice == iCell)
                good_str = ' resp';
            else
                good_str = ' not resp';
            end
            title(['Cell # ' num2str(iCell) ' is' good_str])
        end
    end
end
suptitle([mouse ' ' date '- choiceAlign- by side and contrast- left is purple; right is blue'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide&Contrast_byCell.pdf']), '-dpdf', '-bestfit')

figure %stim response
for icon = 1:ncon
    subplot(2,2,1)
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));%successful, right, contrast loop
    indRH(1,icon) = length(ind1); %how many fit the above parameters
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3))
    hold on
    subplot(2,2,2)
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));%fail, right, contrast loop
    indRM(1,icon) = length(ind2); %how many
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3))
    hold on
    subplot(2,2,3)
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));%successful, left, contrast loop
    indLH(1,icon) = length(ind3); %how many
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind3),2),3))
    hold on
    subplot(2,2,4)
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));%fail, left, contrast loop
    indLM(1,icon) = length(ind4);%how many
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind4),2),3))
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
xlabel('msec')
ylabel('dF/F')
ylim([-.05 .2])
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
xlabel('msec')
ylabel('dF/F')
ylim([-.05 .2])
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
xlabel('msec')
ylabel('dF/F')
ylim([-.05 .2])
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
xlabel('msec')
ylabel('dF/F')
ylim([-.05 .2])
suptitle([mouse ' ' date '- stimAlign'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byOutcome.pdf']), '-dpdf', '-bestfit')

figure %choice response
for icon = 1:ncon
    subplot(2,2,1)
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));%successful,right,contrast loop
    indRH(1,icon) = length(ind1); %how many
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3))
    hold on
    subplot(2,2,2)
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));%fail,right,contrast loop
    indRM(1,icon) = length(ind2); %how many
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3))
    hold on
    subplot(2,2,3)
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));%successful,left,contrast
    indLH(1,icon) = length(ind3); %how many
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind3),2),3))
    hold on
    subplot(2,2,4)
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));%fail,left,contrast
    indLM(1,icon) = length(ind4); %how many
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind4),2),3))
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
ylim([-.05 .2])
xlabel('msec')
ylabel('dF/F')
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
ylim([-.05 .2])
xlabel('msec')
ylabel('dF/F')
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
ylim([-.05 .2])
xlabel('msec')
ylabel('dF/F')
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
ylim([-.05 .2])
xlabel('msec')
ylabel('dF/F')
suptitle([mouse ' ' date '- choiceAlign'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byOutcome.pdf']), '-dpdf')


%tt = (-19:30)*1000/frame_rate;%redundant
%--------unclear A:probability of left, attention task, get used to always
%looking one way then switch
if input.doRandProb
    figure  %?
    tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);  %?
    nprob = length(input.ProbList);
    probs = cell2mat(input.ProbList);
    indR = zeros(nprob,ncon);
    indL = zeros(nprob,ncon);
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));%prob and contrast loop, successful,right
            indR(iprob,icon) = length(ind1);
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));%prob and contrast loop, successful,left
            indL(iprob,icon) = length(ind2);
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')
    
    figure
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- choiceAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')
    
    figure
    tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
    nprob = length(input.ProbList);
    probs = cell2mat(input.ProbList);
    indR = zeros(nprob,ncon);
    indL = zeros(nprob,ncon);
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            indR(iprob,icon) = length(ind1);
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));
            indL(iprob,icon) = length(ind2);
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb_Miss.pdf']), '-dpdf', '-bestfit')
    
    figure
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- choiceAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byProb_Miss.pdf']), '-dpdf', '-bestfit')
end
%----------
%% Wheel analysis

Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));%ignore index
Tix = setdiff(1:length(input.trialOutcomeCell), Iix);%not ignore index
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);%max decision time in ms

qVals_final = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_act = nan(18001, uint16(length(input.trialOutcomeCell)));
qTimes_thresh = nan(1, uint16(length(input.trialOutcomeCell)));
cVals_thresh = nan(1, uint16(length(input.trialOutcomeCell)));

for trN = 1:length(input.trialOutcomeCell)-1
    if find(Tix == trN)
        qTimes = double([input.quadratureTimesUs{trN} input.quadratureTimesUs{trN+1}]./1000);
        qVals = double([input.quadratureValues{trN} input.quadratureValues{trN+1}]);
        cTimes = double(input.counterTimesUs{trN}./1000);
        cVals = double(input.counterValues{trN});
        stimTime = double(input.stimTimestampMs{trN});
        %stimVal = double(input.qStimOn{trN});
        qTimes_zero = qTimes-stimTime;
        qVals = qVals-qVals(1);
        time_ind = find(qTimes_zero>= -8000 & qTimes_zero<=10000);
        if length(time_ind)>2
            qTimes_sub = qTimes_zero(time_ind);
            qVals_sub = qVals(time_ind);
            qTimes_temp = qTimes(time_ind);
            rep_ind = find(diff(qTimes_sub)==0);
            qTimes_sub(rep_ind) = [];
            qVals_sub(rep_ind) = [];
            qTimes_temp(rep_ind) = [];
            qTimes_final = -8000:10000;
            qTimes_act(:,trN) = interp1(qTimes_temp, qTimes_temp, qTimes_final+stimTime)';
            qVals_final(:,trN) = interp1(qTimes_sub, qVals_sub, qTimes_final)';
            if input.tDecisionTimeMs{trN} < 10000
                if isnan(qVals_final(8000,trN))
                    qVals_final(8000,trN) = qVals_final(find(~isnan(qVals_final(:,trN)),1,'first'),trN);
                end
                qVal_off = qVals_final(:,trN)-qVals_final(8000,trN);
                qTimes_thresh(:,trN) = qTimes_act(8000+find(abs(qVal_off(8000:end,:))>5,1,'first'),trN);
                cVals_thresh(:,trN) = cVals(:,find(cTimes>qTimes_thresh(:,trN),1,'first'));
            end
        else
            return
        end
    end
end

minR = input.tooFastTimeMs;
figure;
for it = 1:30
    subplot(5,6,it)
    plot(qTimes_final, qVals_final(:,it))
    xlim([-100 input.tDecisionTimeMs{it}])
    vline(minR)
    if input.tLeftTrial{it}
        title(['Left ' input.trialOutcomeCell{it}])
    else
        title(['Right ' input.trialOutcomeCell{it}])
    end
end


SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'incorrect');
left = cell2mat(input.tLeftTrial);
maxR = input.reactionTimeMs;
lowR = intersect(find(cell2mat(input.tDecisionTimeMs)<maxR), find(cell2mat(input.tDecisionTimeMs)>minR));

qVals_offset = bsxfun(@minus, qVals_final, qVals_final(8001,:));
figure;
subplot(2,2,1)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(SIx),find(left)))))
xlim([-500 2000])
vline(minR)
title('Correct Left Trials')

subplot(2,2,2)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(SIx),find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Correct Right Trials')

subplot(2,2,3)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(MIx),find(left)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Left Trials')

subplot(2,2,4)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(find(MIx),find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Right Trials')

suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])

print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome.pdf']), '-dpdf')

figure
con_str = strvcat('b', 'r', 'y');
subplot(2,2,1)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(SIx), find(left)))),2), nanstd(qVals_offset(:,(intersect(find(SIx), find(left)))),[],2)./sqrt(length(intersect(find(SIx), find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(SIx), find(left==0)))),2), nanstd(qVals_offset(:,(intersect(find(SIx), find(left==0)))),[],2)./sqrt(length(intersect(find(SIx), find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all correct trials')

subplot(2,2,2)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(MIx), find(left)))),2), nanstd(qVals_offset(:,(intersect(find(MIx), find(left)))),[],2)./sqrt(length(intersect(find(MIx), find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(find(MIx), find(left==0)))),2), nanstd(qVals_offset(:,(intersect(find(MIx), find(left==0)))),[],2)./sqrt(length(intersect(find(MIx), find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all incorrect trials')

subplot(2,2,3)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(find(SIx), find(left))));
    shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
    hold on;
end
xlim([-500 2000])
vline(minR)
title('Avg left correct trials by contrast')

subplot(2,2,4)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(find(SIx), find(left==0))));
    shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
    hold on
end
xlim([-500 2000])
vline(minR)
title('Avg right correct trials by contrast')
suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome_avg.pdf']), '-dpdf')

%% eyetracking

calib = 1/26.6; %mm per pixel
pre_event_time = 1000;
post_event_time = 4000;
preevent_frames = ceil(pre_event_time*(frame_rate/1000));
postevent_frames = ceil(post_event_time*(frame_rate/1000));

% Load and combine eye tracking data
% Set current directory to crash folder
Area = {};
Centroid = {};
Eye_data = {};
for irun =  1:nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    fn = [ImgFolder(irun,:) '_000_000_eye.mat'];
    data = load(fn);          % should be a '*_eye.mat' file
    
    data = squeeze(data.data);      % the raw images...
    xc = size(data,2)/2;       % image center
    yc = size(data,1)/2;
    W=40;
    
    rad_range = [6 15];
    data = data(yc-W:yc+W,xc-W:xc+W,:);
    warning off;
    
    A = cell(size(data,3),1);
    B = cell(size(data,3),1);
    for n = 1:size(data,3)
        A{n} = [0,0];
        B{n} = [0];
    end
    eye = struct('Centroid',A,'Area',B);
    radii = [];
    for n = 1:size(data,3)
        [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
        if(isempty(center))
            eye(n).Centroid = [NaN NaN];    % could not find anything...
            eye(n).Area = NaN;
        else
            [~,idx] = max(metric);          % pick the circle with best score
            eye(n).Centroid = center(idx,:);
            eye(n).Area = pi*radii(idx)^2;
        end
        if mod(n,100)==0
            fprintf('Frame %d/%d\n',n,size(data,3));
        end
    end
    Centroid{irun} = cell2mat({eye.Centroid}');
    Area{irun} = cell2mat({eye.Area}');
    Eye_data{irun} = data;
end

%reset frame counter
run_trials = input.trialsSinceReset;

cStimOn = cell2mat(input.cStimOn);
cDecision = cell2mat(input.cDecision);

Area_temp = [];
Centroid_temp = [];
Eye_data_temp = [];
if nrun > 1
    for irun = 1:nrun
        if irun < nrun
            offset = size(Area{irun},1);
            startTrial = run_trials(irun)+1;
            endTrial = run_trials(irun)+run_trials(irun+1);
            cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
            cDecision(1,startTrial:endTrial) = cDecision(1,startTrial:endTrial)+offset;
        end
        Area_temp = [Area_temp; Area{irun}];
        Centroid_temp = [Centroid_temp; Centroid{irun}];
        Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
    end
else
    Area_temp = [Area_temp; Area{irun}];
    Centroid_temp = [Centroid_temp; Centroid{irun}];
    Eye_data_temp = cat(3, Eye_data_temp, Eye_data{irun});
end
clear Eye_data;
ntrials = length(input.trialOutcomeCell);

% no measurement frames
figure;
hist(sqrt(Area_temp./pi));
figure;
x = find(isnan(Area_temp));
if length(x)>25
    minx = 25;
else
    minx = length(x);
end
start = 1;
frames = sort(randsample(length(x),minx));
for i = 1:minx
    subplot(5,5,start);
    imagesq(Eye_data_temp(:,:,x(frames(i))));
    title(x(frames(i)))
    start = start+1;
end

%align eyetracking to
nanrun = ceil(500*(frame_rate/1000));
Rad_temp = sqrt(Area_temp./pi);
sz = size(Eye_data_temp);
rad_mat_start = zeros(preevent_frames+postevent_frames, ntrials);
centroid_mat_start = zeros(preevent_frames+postevent_frames,2, ntrials);
eye_mat_start = zeros(sz(1), sz(2), preevent_frames+postevent_frames, ntrials);
rad_mat_decide = zeros(preevent_frames+postevent_frames, ntrials);
centroid_mat_decide = zeros(preevent_frames+postevent_frames,2, ntrials);
eye_mat_decide = zeros(sz(1), sz(2), preevent_frames+postevent_frames, ntrials);
nframes = size(Rad_temp,1);
for itrial = 1:ntrials
    if itrial == ntrials
        crange = [double(cStimOn(itrial))-preevent_frames:nframes];
    else
        crange = [double(cStimOn(itrial))-preevent_frames: double(cStimOn(itrial+1)-preevent_frames-1)];
    end
    if sum(isnan(Rad_temp(crange,1)),2)>0
        if length(find(tsmovavg(isnan(Rad_temp(crange,1)), 's', nanrun, 1) == 1))> 0
            Rad_temp(crange,1) = NaN(length(crange),1);
        else
            nanind = find(isnan(Rad_temp(crange,1)));
            dataind = find(~isnan(Rad_temp(crange,1)));
            for inan = 1:length(nan_ind)
                gap = min(abs(nan_ind(inan)-data_ind),[],1);
                good_ind_stim = find(abs(nan_ind(inan)-data_ind) == gap);
                Rad_temp(nan_ind(inan),1) = mean(Rad_temp(data_ind(good_ind_stim),1),1);
                Centroid_temp(nan_ind(inan),:) = mean(Centroid_temp(data_ind(good_ind_stim),:),1);
            end
        end
    end
    if itrial < ntrials
        rad_mat_start(:,itrial) = Rad_temp(1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames,:);
        centroid_mat_start(:,:,itrial) = Centroid_temp(1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames,:);
        rad_mat_decide(:,itrial) = Rad_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
        centroid_mat_decide(:,:,itrial) = Centroid_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
        eye_mat_start(:,:,:,itrial) = Eye_data_temp(:,:,1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames);
        eye_mat_decide(:,:,:,itrial) = Eye_data_temp(:,:,1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames);
    else
        if (cStimOn(itrial)+postevent_frames)<nframes
            rad_mat_start(:,itrial) = Rad_temp(1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames,:);
            centroid_mat_start(:,:,itrial) = Centroid_temp(1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames,:);
            eye_mat_start(:,:,:,itrial) = Eye_data_temp(:,:,1+cStimOn(itrial)-preevent_frames:cStimOn(itrial)+postevent_frames);
        else
            rad_mat_start(:,itrial) = nan(preevent_frames+postevent_frames,1);
            centroid_mat_start(:,:,itrial) = nan(preevent_frames+postevent_frames,2,1);
            eye_mat_start(:,:,:,itrial) = nan(sz(1),sz(2),preevent_frames+postevent_frames,1);
        end
        if (cDecision(itrial)+postevent_frames)<nframes
            rad_mat_decide(:,itrial) = Rad_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
            centroid_mat_decide(:,:,itrial) = Centroid_temp(1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames,:);
            eye_mat_decide(:,:,:,itrial) = Eye_data_temp(:,:,1+cDecision(itrial)-preevent_frames:cDecision(itrial)+postevent_frames);
        else
            rad_mat_decide(:,itrial) = nan(preevent_frames+postevent_frames,1);
            centroid_mat_decide(:,:,itrial) = nan(preevent_frames+postevent_frames,2,1);
            eye_mat_decide(:,:,:,itrial) = nan(sz(1),sz(2),preevent_frames+postevent_frames,1);
        end
    end
    
end
rad_mat_start = bsxfun(@times, rad_mat_start, calib);
centroid_mat_start = bsxfun(@times,centroid_mat_start,calib);
rad_mat_decide = bsxfun(@times, rad_mat_decide, calib);
centroid_mat_decide = bsxfun(@times,centroid_mat_decide,calib);

save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str  '_pupil.mat']), 'Area', 'Centroid', 'frame_rate' , 'rad_mat_start','centroid_mat_start','rad_mat_decide','centroid_mat_decide', 'input', 'cDecision', 'cTrialStart' );

%% plot eyetracking data

tLeftTrial = cell2mat(input.tLeftTrial);
tLeftResponse  = cell2mat(input.tLeftResponse);
tRightResponse = cell2mat(input.tRightResponse);
centroid_mat_start = permute(centroid_mat_start,[1,3,2]);
centroid_mat_decide = permute(centroid_mat_decide,[1,3,2]);

rad_mat_norm = rad_mat_start./nanmean(nanmean(rad_mat_start(preevent_frames/2:preevent_frames,:),1),2);
centroid_mat_norm = bsxfun(@minus,centroid_mat_start,nanmean(nanmean(centroid_mat_start(preevent_frames/2:preevent_frames,:,:),1),2));
centroid_mat_norm = permute(centroid_mat_norm,[1,3,2]);

indLcorr = intersect(find(tLeftResponse),find(tLeftTrial));
indLincorr = intersect(find(tLeftResponse),find(~tLeftTrial));
indRcorr = intersect(find(tRightResponse),find(~tLeftTrial));
indRincorr = intersect(find(tRightResponse),find(tLeftTrial));

figure;
subplot(3,2,1)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLincorr),2)', nanstd(rad_mat_norm(:,indLincorr),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLcorr),2)', nanstd(rad_mat_norm(:,indLcorr),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,2)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRincorr),2)', nanstd(rad_mat_norm(:,indRincorr),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRcorr),2)', nanstd(rad_mat_norm(:,indRcorr),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,3)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLincorr)),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLcorr)),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,4)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRincorr)),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRcorr)),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,5)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLincorr)),[],2)./sqrt(length(indLincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLcorr)),[],2)./sqrt(length(indLcorr))','k');
title(['Left Responses- correct: ' num2str(length(indLcorr)) '; incorrect: ' num2str(length(indLincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,6)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRincorr)),[],2)./sqrt(length(indRincorr))','r');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRcorr)),[],2)./sqrt(length(indRcorr))','k');
title(['Right Responses- correct: ' num2str(length(indRcorr)) '; incorrect: ' num2str(length(indRincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])
suptitle([mouse ' ' date '- trial start align'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_byResp.pdf']),'-dpdf', '-bestfit')

indLcorr = intersect(find(tLeftResponse),find(tLeftTrial));
indLincorr = intersect(find(tLeftResponse),find(~tLeftTrial));
indRcorr = intersect(find(tRightResponse),find(~tLeftTrial));
indRincorr = intersect(find(tRightResponse),find(tLeftTrial));

figure;
subplot(3,2,1)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRcorr),2)', nanstd(rad_mat_norm(:,indRcorr),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLcorr),2)', nanstd(rad_mat_norm(:,indLcorr),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,2)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indLincorr),2)', nanstd(rad_mat_norm(:,indLincorr),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(rad_mat_norm(:,indRincorr),2)', nanstd(rad_mat_norm(:,indRincorr),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Radius')
xlabel('Time (ms)')
ylim([0.95 1.2])

subplot(3,2,3)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRcorr)),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLcorr)),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,4)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indLincorr)),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,1,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indRincorr)),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Horizontal Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,5)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRcorr)),[],2)./sqrt(length(indRcorr))','b');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLcorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLcorr)),[],2)./sqrt(length(indLcorr))','g');
title(['Left Trials- correct: ' num2str(length(indLcorr)) '; Right Trials- correct: ' num2str(length(indRcorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])

subplot(3,2,6)
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indLincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indLincorr)),[],2)./sqrt(length(indLincorr))','g');
hold on
shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate,nanmean(squeeze(centroid_mat_norm(:,2,indRincorr)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indRincorr)),[],2)./sqrt(length(indRincorr))','b');
title(['Left Trials- incorrect: ' num2str(length(indLincorr)) '; Right trials- incorrect: ' num2str(length(indRincorr))])
ylabel('Vertical Position')
xlabel('Time (ms)')
ylim([-.05 .15])
suptitle([mouse ' ' date '- trial start align'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_bySide.pdf']),'-dpdf', '-bestfit')

probList = cell2mat(input.ProbList);
tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
col_mat = strvcat('b', 'k', 'g');
indL_n = [];
indR_n = [];
figure;
for iprob = 1:length(probList)
    prob = probList(iprob);
    ind = find(tProbLeft == prob);
    ind = 1+((iprob-1)*80):80+((iprob-1)*80);
    indL = intersect(ind, intersect(find(tLeftResponse),find(tLeftTrial)));
    indR = intersect(ind, intersect(find(tRightResponse),find(~tLeftTrial)));
    indL_n = [indL_n length(indL)];
    indR_n = [indR_n length(indR)];
    subplot(3,2,1)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(rad_mat_norm(:,indL),2)', nanstd(rad_mat_norm(:,indL),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,2)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(rad_mat_norm(:,indR),2)', nanstd(rad_mat_norm(:,indR),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
    subplot(3,2,3)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,1,indL)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indL)),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,4)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,1,indR)),2)', nanstd(squeeze(centroid_mat_norm(:,1,indR)),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
    subplot(3,2,5)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,2,indL)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indL)),[],2)./sqrt(length(indL))',col_mat(iprob));
    hold on
    subplot(3,2,6)
    shadedErrorBar((1-preevent_frames:postevent_frames).*frame_rate, nanmean(squeeze(centroid_mat_norm(:,2,indR)),2)', nanstd(squeeze(centroid_mat_norm(:,2,indR)),[],2)./sqrt(length(indR))',col_mat(iprob));
    hold on
end
subplot(3,2,1)
title(['Left trials- ' num2str(indL_n)])
ylim([.95 1.25])
ylabel('Radius')
subplot(3,2,2)
title(['Right trials- ' num2str(indR_n)])
ylim([.95 1.25])
ylabel('Radius')
subplot(3,2,3)
ylim([-.1 .2])
ylabel('Horizontal Position')
subplot(3,2,4)
ylim([-.1 .2])
ylabel('Horizontal Position')
subplot(3,2,5)
ylim([-.1 .2])
ylabel('Vertical Position')
subplot(3,2,6)
ylim([-.1 .2])
ylabel('Vertical Position')
suptitle([mouse ' ' date '- trial start align; Block order: black, blue, cyan, green'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_byBlock.pdf']),'-dpdf', '-bestfit')

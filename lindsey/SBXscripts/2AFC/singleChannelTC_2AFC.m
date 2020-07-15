%% get path names
date = '181221';
ImgFolder = strvcat('003','004');
time = strvcat('1733');
mouse = 'i1103';
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
end
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\' ;
data_base = 'home\lindsey\Data\2P_images\';
%% load dataset
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [fn_base data_base date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ryan\2Pimages\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW14\two-photon imaging\' date '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    if size(time,1) >= irun
        fName = [fn_base 'Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
        load(fName);
        temp(irun) = input;
        if irun>1
            if isfield(input, 'tLeftTrial')
                ntrials = size(input.trialOutcomeCell,2);
                for itrial = 1:ntrials
                    temp(irun).counterValues{itrial} = bsxfun(@plus,temp(irun).counterValues{itrial},offset);
                    temp(irun).cTrialStart{itrial} = temp(irun).cTrialStart{itrial}+offset;
                    temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                    temp(irun).cDecision{itrial} = temp(irun).cDecision{itrial}+offset;
                end
            end
        end
    end
    
    nframes = min([info.config.frames input.counterValues{end}(end)]);
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    tic
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    toc
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Plot outcome by trial number
SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');

figure;
plot(smooth(SIx,10));
hold on
plot(smooth(MIx,10));

%% Crop data and input struct
%change trial range in trialChopper
input = trialChopper(input,[1 415]);
data = data(:,:,input.counterValues{1}(1):input.counterValues{end}(end));

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

%% Register data

data_avg = mean(data(:,:,20001:20500),3);

if exist(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    load(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability

figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
%% find activated cells

tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
cStimOn = celleqel2mat_padded(input.cStimOn);
cDecision = celleqel2mat_padded(input.cDecision);
SIx = strcmp(input.trialOutcomeCell, 'success');
nTrials = length(tLeftTrial);
sz = size(data_reg);
data_f = nan(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
data_targ_late = nan(sz(1),sz(2),nTrials);
data_resp = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials
    data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3);
    data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3);
    data_targ_late(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5+ceil(input.stimOnTimeMs./frame_rate):cStimOn(itrial)+25+ceil(input.stimOnTimeMs./frame_rate)),3);
    if cDecision(itrial)+25<nframes
        data_resp(:,:,itrial) = mean(data_reg(:,:,cDecision(itrial)+5:cDecision(itrial)+25),3);
    end
end
data_targ_dfof = (data_targ-data_f)./data_f;
data_resp_dfof = (data_resp-data_f)./data_f;
data_targ_late_dfof = (data_targ_late-data_f)./data_f;
indL = find(tLeftTrial);
data_dfof_L = mean(data_targ_dfof(:,:,indL),3);
data_dfof_late_L = mean(data_targ_late_dfof(:,:,indL),3);
indR = intersect(find(tGratingContrast==1),find(~tLeftTrial));
data_dfof_R = mean(data_targ_dfof(:,:,indR),3);
data_dfof_late_R = mean(data_targ_late_dfof(:,:,indR),3);
data_dfof_resp = nanmean(data_resp_dfof(:,:,find(SIx)),3);
figure;
subplot(2,3,1)
imagesc(data_dfof_L)
title('Left Stim')
subplot(2,3,2)
imagesc(data_dfof_R)
title('Right Stim')
subplot(2,3,3)
imagesc(data_dfof_late_L)
title('Left Stim Late')
subplot(2,3,4)
imagesc(data_dfof_late_R)
title('Right Stim Late')
subplot(2,3,5)
imagesc(data_dfof_resp)
title('Decision')
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_dFoF.pdf']), '-dpdf')
figure;
subplot(2,1,1)
imagesc(data_dfof_L)
title('Left Stim')
subplot(2,1,2)
imagesc(data_dfof_R)
title('Right Stim')
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_stimResp.pdf']), '-dpdf')
data_dfof = cat(3, data_dfof_resp,data_dfof_late_L, data_dfof_late_R,data_dfof_L,data_dfof_R);
data_dfof_max = max(data_dfof,[],3);
figure; imagesc(data_dfof_max)
    
%% cell segmentation 

mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell)

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
%np_tc = stackGetTimeCourses(data_reg,mask_np);
clear np_tc np_tc_down
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
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

save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_tc np_tc

%% 2AFC analysis

tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
nside = 2;
tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tGratingContrast);
ncon = length(cons);
cStimOn = celleqel2mat_padded(input.cStimOn);
cDecision = celleqel2mat_padded(input.cDecision);
nTrials = length(tLeftTrial);
nCells = size(npSub_tc,2);
nframes = size(npSub_tc,1);
data_stim = nan(50,nCells,nTrials);
data_choice = nan(50,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
    end
    if cDecision(itrial)
        if cDecision(itrial)+29 < nframes
            data_choice(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:);
        end
    end
end
dataf = mean(data_stim(1:20,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
data_choice_dfof = bsxfun(@rdivide, bsxfun(@minus, data_choice, dataf), dataf);

ind = intersect(find(~tLeftTrial),find(tGratingContrast ==1));
figure; plot(mean(mean(data_stim_dfof(:,:,ind),2),3));

SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');
FIx = strcmp(input.trialOutcomeCell, 'incorrect');

h_stim = zeros(nside, nCells);
p_stim = zeros(nside, nCells);
h_choice = zeros(nside, nCells);
p_choice = zeros(nside, nCells);
base_win = [19:23];
resp_win = [25:29];
for iS = 1:nside
    indS = intersect(find(SIx), find(tLeftTrial == iS-1));
    indC = intersect(indS, find(tGratingContrast == 1));
    [h_stim(iS,:) p_stim(iS,:)] = ttest(squeeze(mean(data_stim_dfof(resp_win,:,indC),1))', squeeze(mean(data_stim_dfof(base_win,:,indC),1))', 'tail', 'right');
    [h_choice(iS,:) p_choice(iS,:)] = ttest(squeeze(mean(data_choice_dfof(resp_win,:,indC),1))', squeeze(mean(data_choice_dfof(base_win,:,indC),1))', 'tail', 'right');
end
good_ind_stim = find(sum(h_stim,1)>0);
good_ind_stim_side{1} = find(h_stim(1,:)>0);
good_ind_stim_side{2} = find(h_stim(2,:)>0);
good_ind_choice = find(sum(h_choice,1)>0);

tt = (-20:29)*(1000/frame_rate);
figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
       indS = intersect(find(tGratingContrast == 1), intersect(find(SIx), find(tLeftTrial == iS-1)));
       plot(tt,nanmean(data_stim_dfof(:,iCell,indS),3));
       hold on;
    end
    vline(([resp_win(1) resp_win(end)]-20)*(1000/frame_rate))
    if find(find(sum(h_stim,1)== 2) == iCell)
        good_str = 'Resp Both';
    elseif find(good_ind_stim_side{1} == iCell)
        good_str = 'Resp Contra';
    elseif find(good_ind_stim_side{2} == iCell)
        good_str = 'Resp Ipsi';
    else
        good_str = 'Not resp';
    end
    title(good_str)
end
suptitle([mouse ' ' date '- stimAlign- blue is right; red is left'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byCell.jpg']), '-djpeg')

cmap = [blues(ncon); reds(ncon)];
figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        for icon = 1:ncon
           indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));
           plot(tt,nanmean(data_stim_dfof(:,iCell,indS),3),'Color',cmap(icon+((iS-1)*ncon),:));
           hold on;
           if find(find(sum(h_stim,1)== 2) == iCell)
                good_str = 'Resp Both';
            elseif find(good_ind_stim_side{1} == iCell)
                good_str = 'Resp Contra';
            elseif find(good_ind_stim_side{2} == iCell)
                good_str = 'Resp Ipsi';
            else
                good_str = 'Not resp';
            end
           %title(good_str)
           axis off
        end
    end
end
suptitle([mouse ' ' date '- stimAlign- by side and contrast- blue is right; red is left'])
h=gcf;
set(h,'PaperOrientation','landscape');
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide&Contrast_byCell.pdf']), '-dpdf', '-bestfit')

% for iCell = 1:nCells
% figure;
%     for iS = 1:nside
%         for icon = 1:ncon
%            subplot(nside,ncon,icon+((iS-1)*ncon))
%            indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));
%            plot(tt, squeeze(data_stim_dfof(:,iCell,indS)));
%         end
%     end
%     suptitle(['Cell #' num2str(iCell)])
% end

figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
       indS = intersect(find(SIx), find(tLeftTrial == iS-1));
       plot(tt,nanmean(data_choice_dfof(:,iCell,indS),3));
       hold on;
    end
    vline(([resp_win(1) resp_win(end)]-20)*(1000/frame_rate))
    if find(good_ind_choice == iCell)
        good_str = 'Resp';
    else
        good_str = 'Not resp';
    end
    title(good_str)
end
suptitle([mouse ' ' date '- choiceAlign- blue is right; red is left'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
        for icon = 1:ncon
           indS = intersect(find(tGratingContrast == cons(icon)), intersect(find(SIx), find(tLeftTrial == iS-1)));
           plot(tt,nanmean(data_choice_dfof(:,iCell,indS),3),'Color',cmap(icon+((iS-1)*ncon),:));
           hold on;
           if find(good_ind_choice == iCell)
                good_str = 'Resp';
           else
                good_str = 'Not resp';
           end
           title(good_str)
        end
    end
end
suptitle([mouse ' ' date '- choiceAlign- by side and contrast- left is red; right is blue'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide&Contrast_byCell.pdf']), '-dpdf', '-bestfit')

figure; 
for icon = 1:ncon
    ax(1) = subplot(2,2,1);
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRH(1,icon) = length(ind1); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:));
    hold on
    ax(2) = subplot(2,2,2);
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRM(1,icon) = length(ind2); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:));
    hold on
    ax(3) = subplot(2,2,3);
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLH(1,icon) = length(ind3); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind3),2),3),'Color',cmap(icon+ncon,:));
    hold on
    ax(4) = subplot(2,2,4);
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLM(1,icon) = length(ind4); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind4),2),3),'Color',cmap(icon+ncon,:));
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
linkaxes(ax, 'y')
suptitle([mouse ' ' date '- stimAlign'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byOutcome.pdf']), '-dpdf', '-bestfit')
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byOutcome.jpg']), '-djpeg')

figure; 
for icon = 1:ncon
    ax(1) = subplot(2,2,1);
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRH(1,icon) = length(ind1); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:))
    hold on
    ax(2) = subplot(2,2,2);
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRM(1,icon) = length(ind2); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:))
    hold on
    ax(3) = subplot(2,2,3);
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLH(1,icon) = length(ind3); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind3),2),3),'Color',cmap(icon+ncon,:))
    hold on
    ax(4) = subplot(2,2,4);
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLM(1,icon) = length(ind4); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind4),2),3),'Color',cmap(icon+ncon,:))
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
suptitle([mouse ' ' date '- choiceAlign'])
linkaxes(ax, 'y')
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byOutcome.pdf']), '-dpdf')


if input.doRandProb
    figure;
    tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
    nprob = length(input.ProbList);
    probs = cell2mat(input.ProbList);
    indR = zeros(nprob,ncon);
    indL = zeros(nprob,ncon);
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            indR(iprob,icon) = length(ind1); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            indL(iprob,icon) = length(ind2); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')

    figure;
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- choiceAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')

    figure;
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
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));           
            indL(iprob,icon) = length(ind2); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Left Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb_Miss.pdf']), '-dpdf', '-bestfit')

    figure;
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind1),2),3),'Color',cmap(icon,:))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind_stim,ind2),2),3),'Color',cmap(icon,:))
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

%% Wheel analysis

Iix = find(strcmp(input.trialOutcomeCell, 'ignore'));
Nix = find(celleqel2mat_padded(input.didNoGo));
Tix = setdiff(1:length(input.trialOutcomeCell), [Nix Iix]);
maxD = max(cell2mat(input.tDecisionTimeMs(Tix)),[],2);
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


SIx = setdiff(find(strcmp(input.trialOutcomeCell, 'success')),Nix);
MIx = setdiff(find(strcmp(input.trialOutcomeCell, 'incorrect')),Nix);
left = celleqel2mat_padded(input.tLeftTrial);
maxR = input.reactionTimeMs;
lowR = intersect(find(celleqel2mat_padded(input.tDecisionTimeMs)<maxR), find(celleqel2mat_padded(input.tDecisionTimeMs)>minR));

qVals_offset = bsxfun(@minus, qVals_final, qVals_final(8001,:));
figure;
subplot(2,2,1)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(SIx,find(left)))))
xlim([-500 2000])
vline(minR)
title('Correct Left Trials')

subplot(2,2,2)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(SIx,find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Correct Right Trials')

subplot(2,2,3)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(MIx,find(left)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Left Trials')

subplot(2,2,4)
plot(qTimes_final, qVals_offset(:,intersect(lowR, intersect(MIx,find(left==0)))))
xlim([-500 2000])
vline(minR)
title('Incorrect Right Trials')

suptitle(['Mouse ' num2str(input.subjectNum) '; React range: ' num2str(minR) '-' num2str(maxR) ' ms'])

print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Wheel_bySide_byOutcome.pdf']), '-dpdf')

figure
con_str = strvcat('g', 'c', 'b', 'k');
cons(find(cons<0.0001)) = [];
ncon = length(cons);
subplot(2,2,1)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(SIx, find(left)))),2), nanstd(qVals_offset(:,(intersect(SIx, find(left)))),[],2)./sqrt(length(intersect(SIx, find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(SIx, find(left==0)))),2), nanstd(qVals_offset(:,(intersect(SIx, find(left==0)))),[],2)./sqrt(length(intersect(SIx, find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all correct trials')

subplot(2,2,2)
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(MIx, find(left)))),2), nanstd(qVals_offset(:,(intersect(MIx, find(left)))),[],2)./sqrt(length(intersect(MIx, find(left)))),'b');
hold on;
shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,intersect(lowR,intersect(MIx, find(left==0)))),2), nanstd(qVals_offset(:,(intersect(MIx, find(left==0)))),[],2)./sqrt(length(intersect(MIx, find(left==0)))),'r');
xlim([-500 2000])
vline(minR)
title('Avg all incorrect trials')

subplot(2,2,3)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(SIx, find(left))));
    if length(ind)>0
        shadedErrorBar(qTimes_final, nanmean(qVals_offset(:,ind),2), nanstd(qVals_offset(:,ind),[],2)./sqrt(length(ind)),con_str(icon));
        hold on;
    end
end
xlim([-500 2000])
vline(minR)
title('Avg left correct trials by contrast')

subplot(2,2,4)
for icon = 1:ncon
    ind = intersect(find(tGratingContrast == cons(icon)), intersect(lowR,intersect(SIx, find(left==0))));
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
    
    save(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str  '_pupil.mat']), 'Area', 'Centroid', 'frame_rate' , 'rad_mat_start','centroid_mat_start','rad_mat_decide','centroid_mat_decide', 'input', 'cDecision', 'cStimOn' );

    %% plot eyetracking data
    
tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
tLeftResponse  = celleqel2mat_padded(input.tLeftResponse);
tRightResponse = celleqel2mat_padded(input.tRightResponse);
centroid_mat_start = permute(centroid_mat_start,[1,3,2]);
centroid_mat_decide = permute(centroid_mat_decide,[1,3,2]);

rad_mat_norm = rad_mat_start./nanmean(nanmean(rad_mat_start(preevent_frames/2:preevent_frames,:),1),2);
centroid_mat_norm = bsxfun(@minus,centroid_mat_start,nanmean(nanmean(centroid_mat_start(preevent_frames/2:preevent_frames,:,:),1),2));
centroid_mat_norm = permute(centroid_mat_norm,[1,3,2]);

indLcorr = setdiff(intersect(find(tLeftResponse),find(tLeftTrial)),Nix);
indLincorr = setdiff(intersect(find(tLeftResponse),find(~tLeftTrial)),Nix);
indRcorr = setdiff(intersect(find(tRightResponse),find(~tLeftTrial)),Nix);
indRincorr = setdiff(intersect(find(tRightResponse),find(tLeftTrial)),Nix);

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

indLcorr = setdiff(intersect(find(tLeftResponse),find(tLeftTrial)),Nix);
indLincorr = setdiff(intersect(find(tLeftResponse),find(~tLeftTrial)),Nix);
indRcorr = setdiff(intersect(find(tRightResponse),find(~tLeftTrial)),Nix);
indRincorr = setdiff(intersect(find(tRightResponse),find(tLeftTrial)),Nix);

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
suptitle([mouse ' ' date '- trial start align; left: green; right: blue'])
print(fullfile(['\\CRASH.dhe.duke.edu\data\home\' tDir '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_startAlign_EyeTC_bySide.pdf']),'-dpdf', '-bestfit')

if doRandProb
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
end
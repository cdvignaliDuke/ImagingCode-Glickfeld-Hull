%% get path names
close all;clear all;clc;

ds = 'MovingDotDetect_ExptList';
iexp = 7; 
rc = behavConstsAV;
eval(ds)

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).dotFolder;
time = expt(iexp).dotTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging retinotopy analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
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
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp
toc

%% For behavior experiments
% % Plot outcome by trial number
% SIx = strcmp(input.trialOutcomeCell, 'success');
% MIx = strcmp(input.trialOutcomeCell, 'ignore');
% 
% figure;
% plot(smooth(SIx,10));
% hold on
% plot(smooth(MIx,10));
% 
% % Crop data and input struct
% input = trialChopper(input,[1 200]);
% data = data(:,:,input.counterValues{1}(1):input.counterValues{end}(end));

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end

data_avg = mean(data(:,:,20001:20500),3);
%% Register data

if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(double(data),[],[],out_bx);
    clear out outs
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

cTarget = celleqel2mat_padded(input.cTargetOn);
cStart = celleqel2mat_padded(input.cLeverDown);
nTrials = length(cTarget);
sz = size(data_reg);
data_f_base = zeros(sz(1),sz(2),nTrials);
data_base = nan(sz(1),sz(2),nTrials);
data_f_targ = zeros(sz(1),sz(2),nTrials);
data_targ = nan(sz(1),sz(2),nTrials);
for itrial = 1:nTrials    
    if ~isnan(cStart(itrial))
        if cStart(itrial)+19 < sz(3)
            data_f_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
            data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+5:cStart(itrial)+10),3);
        end
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+19 < sz(3)
            data_f_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)-20:cTarget(itrial)-1),3);
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+5:cTarget(itrial)+10),3);
        end
    end
    
end
data_base_dfof = (data_base-data_f_base)./data_f_base;
data_targ_dfof = (data_targ-data_f_targ)./data_f_targ;
targSpeed = celleqel2mat_padded(input.tDotSpeedDPS);
spds = unique(targSpeed);
nSpd = length(spds);
baseCon = celleqel2mat_padded(input.tBaseDotContrast);
cons = unique(baseCon);
nCon = length(cons);
data_dfof_spd = zeros(sz(1),sz(2),nSpd,nCon);
[n n2] = subplotn(nSpd);
for iCon = 1:nCon
    figure;
    ind_con = find(baseCon == cons(iCon));
    for ispd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds(ispd)));
        data_dfof_spd(:,:,ispd,iCon) = nanmean(data_targ_dfof(:,:,ind_spd),3);
        subplot(n,n2,ispd)
        imagesc(data_dfof_spd(:,:,ispd,iCon))
        title(spds(ispd))
    end
    suptitle([mouse ' ' date '- Base con = ' num2str(cons(iCon))])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVbySpd_Con' num2str(cons(iCon)) '.pdf']),'-dpdf','-bestfit')
end

ind_con = find(baseCon == 1);
data_dfof_base = nanmean(data_base_dfof(:,:,ind_con),3);

data_dfof = cat(3, reshape(data_dfof_spd,[sz(1), sz(2), nCon*nSpd]), data_dfof_base);
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
%% target analysis
data_targ = nan(120,nCells,nTrials);
data_base = nan(120,nCells,nTrials);
for itrial = 1:nTrials
    if ~isnan(cStart(itrial))
        if cStart(itrial)+99 <= size(npSub_tc,1)
            data_base(:,:,itrial) = npSub_tc(cStart(itrial)-20:cStart(itrial)+99,:);
        end
    end
    if ~isnan(cTarget(itrial))
        if cTarget(itrial)+99 <= size(npSub_tc,1)
            data_targ(:,:,itrial) = npSub_tc(cTarget(itrial)-20:cTarget(itrial)+99,:);
        end
    end
end
data_f_targ = nanmean(data_targ(1:20,:,:),1);
data_dfof_targ = bsxfun(@rdivide,bsxfun(@minus,data_targ,data_f_targ),data_f_targ);
data_f_base = nanmean(data_base(1:20,:,:),1);
data_dfof_base = bsxfun(@rdivide,bsxfun(@minus,data_base,data_f_base),data_f_base);

frameRateHz = double(input.frameRateHz);

base_win =19:21;
resp_win =25:29; 

figure;
if nCells<25
    ii = nCells;
else
    ii = 25;
end
for i = 1:ii
    subplot(5,5,i)
    plot(squeeze(nanmean(mean(data_dfof_targ(20:50,i,:),2),3)))
    vline(base_win-19)
    vline(resp_win-19)
end

tt = [-20:99].*(1000/frameRateHz);
figure;
subplot(2,1,1)
plot(tt,squeeze(nanmean(mean(data_dfof_base,2),3)));
vline((base_win-20).*(1000/frameRateHz))
vline((resp_win-20).*(1000/frameRateHz))
title('Baseline')
subplot(2,1,2)
plot(tt,squeeze(nanmean(mean(data_dfof_targ,2),3)));
vline((base_win-20).*(1000/frameRateHz))
vline((resp_win-20).*(1000/frameRateHz))
title('Target')


%%

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'data_dfof_base', 'data_dfof_targ')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'baseCon', 'cons', 'nCon', 'targSpeed','spds','nSpd', 'nCells', 'frameRateHz', 'nTrials', 'cTarget', 'cStart', 'base_win','resp_win')

%% Speed analysis
late_win = resp_win+(2.*frameRateHz);
full_win = resp_win(1): resp_win(1)+(2.*frameRateHz);
figure;
spd_mat = zeros(nSpd,nCon,nCells,2);
spd_mat_late = zeros(nSpd,nCon,nCells,2);
h_early = zeros(nCon,nSpd,nCells);
p_early = zeros(nCon,nSpd,nCells);
h_late = zeros(nCon,nSpd,nCells);
p_late = zeros(nCon,nSpd,nCells);
h_supp = zeros(nCon,nSpd,nCells);
p_supp = zeros(nCon,nSpd,nCells);
data_dfof_spd = zeros(size(data_dfof_targ,1), nCells, nCon,nSpd);
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    Ax(iCon+(iCon-1)) = subplot(2,2,iCon+(iCon-1));
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,:,ind_spd),3),2))
        data_dfof_spd(:,:,iCon,iSpd) = nanmean(data_dfof_targ(:,:,ind_spd),3);
        hold on
        for iCell = 1:nCells
            [h_early(iCon,iSpd,iCell), p_early(iCon,iSpd,iCell)] = ttest(mean(permute(data_dfof_targ(resp_win,iCell,ind_spd),[1 3 2]),1),mean(permute(data_dfof_targ(base_win,iCell,ind_spd),[1 3 2]),1),'tail','right','alpha',0.05./(nSpd*nCon));
            [h_late(iCon,iSpd,iCell), p_late(iCon,iSpd,iCell)] = ttest(mean(permute(data_dfof_targ(late_win,iCell,ind_spd),[1 3 2]),1),mean(permute(data_dfof_targ(base_win,iCell,ind_spd),[1 3 2]),1),'tail','right','alpha',0.05./(nSpd*nCon));
            [h_supp(iCon,iSpd,iCell), p_supp(iCon,iSpd,iCell)] = ttest(mean(data_dfof_targ(late_win,iCell,ind_spd),1),mean(data_dfof_targ(base_win,iCell,ind_spd),1),'tail','left','alpha',0.05./(nSpd*nCon));

            spd_mat(iSpd,iCon,iCell,1) = mean(nanmean(data_dfof_targ(resp_win,iCell,ind_spd),3),1);
            spd_mat(iSpd,iCon,iCell,2) = nanstd(nanmean(data_dfof_targ(resp_win,iCell,ind_spd),1),[],3)./sqrt(length(ind_spd));
            spd_mat_late(iSpd,iCon,iCell,1) = mean(nanmean(data_dfof_targ(late_win,iCell,ind_spd),3),1);
            spd_mat_late(iSpd,iCon,iCell,2) = nanstd(nanmean(data_dfof_targ(late_win,iCell,ind_spd),1),[],3)./sqrt(length(ind_spd));
        end
    end
    y(iCon+(iCon-1),:) = Ax(iCon+(iCon-1)).YLim;
    title(['BaseCon = ' num2str(cons(iCon))]) 
    ylabel('dF/F')
    xlabel('Time from target (ms)')
    if iCon == 1
        legend(num2str(chop(spds,2)'))
    end
    Ax(iCon+(iCon-1)+1) = subplot(2,2,iCon+(iCon-1)+1);
    errorbar(spds, mean(spd_mat(:,iCon,:,1),3), std(spd_mat(:,iCon,:,1),[],3)./sqrt(nCells))
    hold on
    errorbar(spds, mean(spd_mat_late(:,iCon,:,1),3), std(spd_mat_late(:,iCon,:,1),[],3)./sqrt(nCells))
    y(iCon+(iCon-1)+1,:) = Ax(iCon+(iCon-1)+1).YLim;
    ylabel('dF/F')
    xlabel('Speed (DPS)')
    if iCon == 1
        legend({'early', 'late'})
    end
    title(['BaseCon = ' num2str(cons(iCon))]) 
end
min_y = min(y(:,1),[],1);
max_y = max(y(:,2),[],1);
set(Ax(1:4),'YLim',[min_y max_y])
suptitle([date ' ' mouse '- Avg speed tuning- n = ' num2str(nCells)])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgSpeedResp.pdf']),'-dpdf', '-bestfit')

good_ind = unique([find(sum(sum(h_early,1),2)); find(sum(sum(h_late,1),2))]);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cellResp.mat']),'data_dfof_spd', 'good_ind', 'tt')

for iCell = 1:length(good_ind)
    figure;
    iC = good_ind(iCell);
    y = zeros(4,2);
    clear Ax
    for iCon = 1:nCon
        ind_con = find(baseCon == cons(iCon));
        Ax(iCon+(iCon-1)) = subplot(2,2,iCon+(iCon-1));
        for iSpd = 1:nSpd
            ind_spd = intersect(ind_con, find(targSpeed==spds(iSpd)));
            plot(tt, nanmean(data_dfof_targ(:,iC,ind_spd),3))
            hold on
        end
        y(iCon+(iCon-1),:) = Ax(iCon+(iCon-1)).YLim;
        title(['BaseCon = ' num2str(cons(iCon))]) 
        ylabel('dF/F')
        xlabel('Time from target (ms)')
        if iCon == 1
            legend(num2str(chop(spds,2)'))
        end
        Ax(iCon+(iCon-1)+1) = subplot(2,2,iCon+(iCon-1)+1);
        errorbar(spds, spd_mat(:,iCon,iC,1), spd_mat(:,iCon,iC,2))
        hold on
        errorbar(spds, spd_mat_late(:,iCon,iC,1), spd_mat_late(:,iCon,iC,2))
        y(iCon+(iCon-1)+1,:) = Ax(iCon+(iCon-1)+1).YLim;
        ylabel('dF/F')
        xlabel('Speed (DPS)')
        if iCon == 1
            legend({'early', 'late'})
        end
        title(['BaseCon = ' num2str(cons(iCon))]) 
    end
    min_y = min(y(:,1),[],1);
    max_y = max(y(:,2),[],1);
    set(Ax(1:4),'YLim',[min_y max_y])
    suptitle([date ' ' mouse '- Cell #' num2str(iC)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgSpeedResp_Cell#' num2str(iC) '.pdf']),'-dpdf', '-bestfit')
end

supp_ind = find(h_supp(end,end,:));
high_resp_ind = find(h_late(end,end,:));
low_resp_ind = find(h_late(end,1,:));

figure;
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,low_resp_ind,ind_spd),3),2))
        hold on
    end
    title('Resp to low speed increment')
end
title(['n = ' num2str(size(low_resp_ind,1))])
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon+2)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,high_resp_ind,ind_spd),3),2))
        hold on
    end
    title('Resp to high speed increment')
end
title(['n = ' num2str(size(high_resp_ind,1))])
for iCon = 1:nCon
    ind_con = find(baseCon == cons(iCon));
    subplot(3,2,iCon+4)
    for iSpd = 1:nSpd
        ind_spd = intersect(ind_con, find(targSpeed==spds(iSpd)));
        plot(tt, mean(nanmean(data_dfof_targ(:,supp_ind,ind_spd),3),2))
        hold on
    end
    title('Supp to high speed increment')
end
title(['n = ' num2str(size(supp_ind,1))])
        

%% Parameters
expt = 'V1_5SFx5TF';
date = '140411';
mouse = 'G023';
run_time = ['21_17_44'; '21_41_29'; '21_17_44'; '21_41_29'];
input_time = ['1618'; '1641'; '1618'; '1641'];

orig_rate = 30;
final_rate = 3;
nrun = [3 4];

%% load input structures
base_path = 'Z:\home\ashley\2P imaging\Data';
mat_path = 'Z:\home\lindsey\Data\2P_images\VisStimFiles';
save_path = 'Z:\home\lindsey\Analysis\2P';

trial_SF = [];
trial_TF = [];
trial_Dir = [];
trial_loc = [];

for iRun = nrun
    input_path = fullfile(mat_path,['data-i' mouse(2:end) '-' date '-' input_time(iRun,:) '.mat']);
    load(input_path);
    trial_SF = [trial_SF cell2mat(input.tGratingSpatialFreqCPD)];
    trial_TF = [trial_TF cell2mat(input.tGratingTemporalFreqCPS)];
    trial_Dir = [trial_Dir cell2mat(input.tGratingDirectionDeg)];
    trial_loc = cat(2,trial_loc, [cell2mat(input.tStimWheelCounter); cell2mat(input.tITIWheelCounter)]);
end

SF_mat = unique(trial_SF);
nSF = length(SF_mat);
TF_mat = unique(trial_TF);
nTF = length(TF_mat);
Dir_mat = unique(trial_Dir);
nDir = length(Dir_mat);

%% Calculated parameters
down = orig_rate./final_rate;
nON = input.nScansOn./down;
nOFF = input.nScansOn./down;

nCond = length(SF_mat)*length(TF_mat)*length(Dir_mat);
%% load data
data_cat = [];
for iRun = nrun
    run_path = fullfile(base_path,[date '_' mouse],['20' date(1:2) '_' date(3:4) '_' date(5:6)],['(20' date '_' run_time(iRun,:) ')'],['(20' date '_' run_time(iRun,:) ')_' expt '_tr' num2str(iRun) '_XYT.raw']);
    data = readrawfile(run_path);
    data_down = stackGroupProject(data,down);
    clear data
    data_cat = cat(3,data_cat,data_down);
    clear data_down
end

%remove negative data (by addition)
data_sub = data_cat-min(min(min(data_cat,[],1),[],2),[],3);
clear data_cat

data_avg = mean(data_sub(:,:,1:50),3);
figure; imagesq(data_avg); colormap(gray)

%register
[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub
siz = size(data_reg);

data_avg = squeeze(mean(mean(data_reg,2),1));
nRep = siz(3)./(nON+nOFF);
data_rep = zeros(nON+nOFF,nRep);
start = 1;
for iRep = 1:nRep
    data_rep(:,iRep) = data_avg(start:start+nON+nOFF-1);
    start= start+nON+nOFF;
end
data_rep_avg = mean(data_rep,2);
figure; plot(data_rep_avg);


%% create dF/F stack


nRep_mat = zeros(1,nCond);
data_sort = zeros(size(data_reg));
data_avg = zeros(siz(1),siz(2), nON+nOFF, nCond);
data_dFoverF = zeros(siz(1),siz(2),nCond);
trial_loc_sort = zeros(size(trial_loc));

iCond = 1;
start = 1;
begin = floor((2/3)*nOFF);
rep = 1;

pre_win = [1 nOFF-begin];
post_win = [nOFF-begin+1 (nOFF+nON-begin)];

for iTF = 1:nTF
    TF_ind = find(trial_TF == TF_mat(iTF));
    for iSF = 1:nSF
        SF_ind = find(trial_SF == SF_mat(iSF));
        SFTF_ind = intersect(TF_ind, SF_ind);
        nRep_mat(1,iCond) = length(SFTF_ind);
        stim_rep = zeros(siz(1),siz(2),nON+nOFF,length(SFTF_ind));
        for iRep = 1:length(SFTF_ind)
            trial = SFTF_ind(iRep);
            if trial == nRep
                data_sort(:,:,start:start+nON+nOFF-1-begin) = data_reg(:,:,begin+1+((trial-1)*(nON+nOFF)):(trial*(nON+nOFF)));
                data_sort(:,:,start+nON+nOFF-begin:start+nON+nOFF-1) = data_reg(:,:,1:begin);
                stim_rep(:,:,1:nON+nOFF-begin,iRep) = data_reg(:,:,begin+1+((trial-1)*(nON+nOFF)):(trial*(nON+nOFF)));
                stim_rep(:,:,nON+nOFF+1-begin:nON+nOFF,iRep) = data_reg(:,:,1:begin);
                start = start+nON+nOFF;
                trial_loc_sort(:,rep) = trial_loc(:,trial);
                rep = rep+1;
            else
                data_sort(:,:,start:start+nON+nOFF-1) = data_reg(:,:,begin+1+((trial-1)*(nON+nOFF)):begin+(trial*(nON+nOFF)));
                stim_rep(:,:,:,iRep) = data_reg(:,:,begin+1+((trial-1)*(nON+nOFF)):begin+(trial*(nON+nOFF)));
                start = start+nON+nOFF;
                trial_loc_sort(:,rep) = trial_loc(:,trial);
                rep = rep+1;
            end
        end
        data_avg(:,:,:,iCond) = mean(stim_rep,4);
        data_dFoverF(:,:,iCond) = squeeze((mean(data_avg(:,:,post_win(1):post_win(2),iCond),3)-mean(data_avg(:,:,pre_win(1):pre_win(2),iCond),3))./(mean(data_avg(:,:,pre_win(1):pre_win(2),iCond),3)));
        iCond = iCond+1;
    end
end
clear data_reg
max_dF = max(data_dFoverF,[],3);
figure; imagesq(max_dF); colormap(gray)
data_avg_long = reshape(data_avg, [siz(1) siz(2) (nOFF+nON)*nCond]);

writetiff(data_dFoverF, fullfile(save_path, [date '_' mouse], [date '_' mouse '_dFoverF.tif']));
%% use max dF/F to find ROIS

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);
figure; imagesq(mask_cell);

%timecourses
data_TC = stackGetTimeCourses(data_sort,mask_cell);
figure; tcOffsetPlot(data_TC)

data_TCavg = stackGetTimeCourses(data_avg_long,mask_cell);
figure; tcOffsetPlot(data_TCavg)

nCell = size(data_TC,2);

% data_dF = bsxfun(@minus,data_sort,avgF);
% data_dFoverF = bsxfun(@rdivide,data_dF,avgF);
% clear data_dF
% data_TCdF = stackGetTimeCourses(data_dFoverF,mask_cell);
% figure; tcOffsetPlot(data_TCdF)
% 
% data_TCdFavg = stackGetTimeCourses(avgdFoverF,mask_cell);
% figure; tcOffsetPlot(data_TCdFavg)
%% analyze by stimulus type
resp = struct;
start = 1;
for iCond = 1:nCond
    reps = nRep_mat(1,iCond);
    resp(iCond).off = zeros(reps, nCell);
    resp(iCond).on = zeros(reps, nCell);
    resp(iCond).dF = zeros(reps, nCell);
    resp(iCond).dFoverF = zeros(reps,nCell);
    resp(iCond).dFoverF_avg = zeros(1,nCell);
    for iRep = 1:reps
        resp(iCond).off(iRep,:) = mean(data_TC(start-1+pre_win(1):start-1+pre_win(2),:),1);
        resp(iCond).on(iRep,:) = mean(data_TC(start-1+post_win(1):start+post_win(2),:),1);
        resp(iCond).dF(iRep,:) = resp(iCond).on(iRep,:)-resp(iCond).off(iRep,:);
        resp(iCond).dFoverF(iRep,:) = resp(iCond).dF(iRep,:)./resp(iCond).off(iRep,:); 
        start= start+nON+nOFF;
    end
    resp(iCond).dFoverF_avg = mean(resp(iCond).dFoverF,1);
end

max_resp = zeros(1,nCell);
for iCell = 1:nCell
    resp_mat_temp = zeros(1,nCond);
    for iCond = 1:nCond
        resp_mat_temp(1,iCond) = resp(iCond).dFoverF_avg(:,iCell);
    end
    max_resp(:,iCell) = max(resp_mat_temp,[],2);
end

alpha = 0.05;
Info_ttest_mat = zeros(nCell, nCond);

for iCell = 1:nCell
    for iCond = 1:nCond
        [H P] = ttest(resp(iCond).on(:,iCell), resp(iCond).off(:,iCell), 'tail', 'right');
        Info_ttest_mat(iCell, iCond) = P;
    end 
end

figure;
start= 1;
for iCell = 1:nCell
    if start>25
        figure;
        start=1;
    end
    subplot(5,5,start)
    for iCond = 1:nCond
        resp_norm = resp(iCond).dFoverF_avg(:,iCell)./max_resp(:,iCell);
        x = ceil(iCond/nTF);
        y = rem(iCond,nSF);
        if y == 0
            y =5;
        end
        if resp_norm>0
            scatter(x,y, ((resp_norm*10).^2), 'ok', 'filled');
            hold on;
        elseif isnan(resp_norm)
            scatter(x,y, 10, 'or', 'filled');
            hold on;
        elseif resp_norm < 0
            scatter(x,y, 10, 'ob', 'filled');
            hold on;
        end
    end
    if min(Info_ttest_mat(iCell,:),[],2)< (alpha/nCond);
    title('**')
    end
    xlim([0 7])
    ylim([0 7])
    axis square
    axis off
    start= start+1;
end
    
%plot data
for iCell = 3;
    figure;
    for iStim = 1:nStim
        subplot(2,3,iStim)
        rep_mat = zeros(nON+nOFF,nRep);
        for iRep = 1:nRep
            plot(1-nOFF:nON,data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell);
        end
        plot(1-nOFF:nON, mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_Az = rem(iStim,size(Az,2));
        if stim_Az == 0
            stim_Az = size(Az,2);
        end
        stim_El= ceil(iStim./size(Az,2));
        title(['Az = ' num2str(Az(1,stim_Az)) ' El = ' num2str(El(1,stim_El))]);
    end
end

        




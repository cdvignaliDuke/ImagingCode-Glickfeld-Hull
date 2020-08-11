%% get path names
close all;clear all;clc;

ds = 'MovingDotDetect_ExptList';
eval(ds)
nexp = length(expt);
iexp = 19;
rc = behavConstsAV;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).gratFolder;
dotFolder = expt(iexp).dotFolder;
time = expt(iexp).gratTime;
nrun = length(ImgFolder);
frameRateHz = params.frameRate;

run_str = catRunName(cell2mat(ImgFolder), nrun);
dot_run_str = catRunName(cell2mat(dotFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging TF analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
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
    %CD = [LG_base '\Data\2P_images\' date '_' mouse '\' ImgFolder{irun}];
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = ntrials*(nOn+nOff);
    
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,nframes);
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    fprintf('Complete')
end
input = concatenateDataBlocks(temp);
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp
clear temp

toc

%% register to dots experiment

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_reg_shifts.mat']))
[out, data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data
%% test stability
figure; 
subplot(2,2,1);
imagesc(data_reg_avg);
title('Grating run avg')
subplot(2,2,2);
imagesc(data_avg)
title('Dot run avg')
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
rgb(:,:,2) = data_avg./max(data_avg(:));
subplot(2,2,3);
image(rgb)
title('Overlay')

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% use dots mask to get TCs

fprintf(['Loading masks from retinotopy runs: ' cell2mat(dotFolder) '\n'])

% loads 'mask_cell', 'mask_np'
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

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

fprintf('Neuropil subtraction complete\n')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%% TF analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(input.tGratingTemporalFreqCPS);
nCells = size(npSub_tc,2);

data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff./2:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

TF_mat = celleqel2mat_padded(input.tGratingTemporalFreqCPS);
TFs = unique(TF_mat);
nTF = length(TFs);
grating_spds = TFs./input.gratingSpatialFreqCPD;

base_win = [nOff-5:nOff];
resp_win = [nOff+4:nOff+nOn];
tt_tf = (1-nOff:nOn).*(1000/frameRateHz);

data_dfof_tf = zeros(nOn+nOff, nCells, nTF);
h_tf = zeros(nCells,nTF);
tf_resp_avg = zeros(nCells,nTF);
tf_resp_mat = [];
tf_list = [];
for iTF = 1:nTF
    ind = find(TF_mat == TFs(iTF));
    data_dfof_tf(:,:,iTF) = mean(data_dfof(:,:,ind),3);
    tf_resp_avg(:,iTF) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
    tf_resp_mat = [tf_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
    tf_list = [tf_list; iTF.*ones(length(ind),1)];
    h_tf(:,iTF) = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nTF-1));
end
TFresp_ind = find(sum(h_tf,2));
h_tfanova = zeros(1,nCells); 
for iCell = 1:nCells
    h_tfanova(iCell) = anova1(tf_resp_mat(:,iCell),tf_list,'off');
end
TFtuned_ind = find(h_tfanova<0.05);
[max_val, max_ind_tf] = max(tf_resp_avg,[],2);

figure; 
cdfplot(grating_spds(max_ind_tf))
hold on
cdfplot(grating_spds(max_ind_tf(TFresp_ind)))
cdfplot(grating_spds(max_ind_tf(TFtuned_ind)))
xlabel('Speed (deg/s)')
ylabel('F(x)')
legend({['all- ' num2str(nCells)], ['responsive- ' num2str(length(TFresp_ind))], ['tuned- ' num2str(length(TFtuned_ind))]},'location', 'southeast')
title([mouse ' ' date ' ' area])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PrefTFdist.pdf']),'-dpdf', '-bestfit')

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'max_ind_tf', 'data_dfof', 'resp_win', 'base_win', 'tt_tf', 'data_dfof_tf', 'h_tf', 'TFresp_ind', 'TFtuned_ind')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'TF_mat', 'TFs', 'nTF', 'grating_spds', 'nOn', 'nOff','frameRateHz')


%% Compare TF and dots

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_trResp.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_singleCellDecode.mat']))

tf_resp_mat = squeeze(mean(data_dfof_tf(resp_win,:,:),1));
nCells = size(tf_resp_mat,1);

r_dots = zeros(1,nCells);
r_tf = zeros(1,nCells);
r_tf2 = zeros(1,nCells);
for iCell = 1:nCells
    nobase_resp = squeeze(spd_resp_mat(iCell,1,:));
    base_resp = squeeze(spd_resp_mat(iCell,2,:));
    if sum(isnan(base_resp))
        base_resp(find(isnan(base_resp))) = [];
    end
    tf_resp = tf_resp_mat(iCell,1:5);
    r_dots(1,iCell) = triu2vec(corrcoef(nobase_resp,base_resp));
    r_tf(1,iCell) = triu2vec(corrcoef(base_resp,tf_resp));
    r_tf2(1,iCell) = triu2vec(corrcoef(nobase_resp,tf_resp));
end


range = (-20:0.001:70);
r_dots_nan = nan(1,length(range));
r_tf_nan = nan(1,length(range));
r_tf2_nan = nan(1,length(range));
cellW_chop = round(cellW,3);
for i = 1:length(range)
    ind = find(cellW_chop==range(i));
    if length(ind)>0 
        if length(ind)>1
            ind = ind(1);
        end
        r_dots_nan(1,i) =r_dots(1,ind);
        r_tf_nan(1,i) = r_tf(1,ind);
        r_tf2_nan(1,i) = r_tf2(1,ind);
    end
end
% fprintf('Smoothing correlation')
% r_dots_smooth = smooth(r_dots_nan,20);
% r_tf_smooth = smooth(r_tf_nan,20);
% r_tf2_smooth = smooth(r_tf2_nan,20);


figure
subplot(4,3,1)
hist(r_dots)
ylabel('Number of cells')
xlabel('Correlation')
title('Dots vs speed increment')
subplot(4,3,4)
hist(r_dots(TFresp_ind))
ylabel('Number of cells')
xlabel('Correlation')
subplot(4,3,7)
scatter(cellW, r_dots)
hold on 
scatter(cellW(TFresp_ind), r_dots(TFresp_ind))
%plot(range,r_dots_smooth)
xlabel('Cell weights')
ylabel('Correlation')
subplot(4,3,2)
hist(r_tf)
ylabel('Number of cells')
xlabel('Correlation')
title('TF vs speed increment')
subplot(4,3,5)
hist(r_tf(TFresp_ind))
ylabel('Number of cells')
xlabel('Correlation')
subplot(4,3,8)
scatter(cellW, r_tf)
hold on
scatter(cellW(TFresp_ind), r_tf(TFresp_ind))
%plot(range,r_tf_smooth)
xlabel('Cell weights')
ylabel('Correlation')
subplot(4,3,3)
hist(r_tf2)
ylabel('Number of cells')
xlabel('Correlation')
title('TF vs Dots')
subplot(4,3,6)
hist(r_tf2(TFresp_ind))
ylabel('Number of cells')
xlabel('Correlation')
subplot(4,3,9)
scatter(cellW, r_tf2)
hold on
scatter(cellW(TFresp_ind), r_tf2(TFresp_ind))
%plot(range,r_tf2_smooth)
xlabel('Cell weights')
ylabel('Correlation')
subplot(4,3,10)
scatter(spds(max_ind), r_dots)
hold on
scatter(spds(max_ind(TFresp_ind)), r_dots(TFresp_ind))
xlabel('Speed Pref (deg/s)')
ylabel('Correlation')
subplot(4,3,11)
scatter(spds(max_ind), r_tf)
hold on
scatter(spds(max_ind(TFresp_ind)), r_tf(TFresp_ind))
xlabel('Speed Pref (deg/s)')
ylabel('Correlation')
subplot(4,3,12)
scatter(spds(max_ind), r_tf2)
hold on
scatter(spds(max_ind(TFresp_ind)), r_tf2(TFresp_ind))
xlabel('Speed Pref (deg/s)')
ylabel('Correlation')


suptitle([mouse ' ' date ' ' area])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_weightVsCorr.pdf']),'-dpdf', '-bestfit')
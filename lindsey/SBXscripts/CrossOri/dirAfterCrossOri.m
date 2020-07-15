%% get path names
close all;clear all;clc;

ds = 'CrossOriSingleStimAdapt_ExptList';
eval(ds)
nexp = length(expt);
iexp = 5;
rc = behavConstsAV;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).dirFolder;
coFolder = expt(iexp).coFolder;
time = expt(iexp).dirTime;
nrun = length(ImgFolder);
frameRateHz = params.frameRate;

run_str = catRunName(cell2mat(ImgFolder), nrun);
dot_run_str = catRunName(cell2mat(coFolder), nrun);

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

%% register to cross-ori experiment

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
title('Direction run avg')
subplot(2,2,2);
imagesc(data_avg)
title('Cross-ori run avg')
sz = size(data_avg);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_reg_avg./max(data_reg_avg(:));
rgb(:,:,2) = data_avg./max(data_avg(:));
subplot(2,2,3);
image(rgb)
title('Overlay')

print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')

%% use cross-ori mask to get TCs

fprintf(['Loading masks from cross-ori runs: ' cell2mat(coFolder) '\n'])

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

fprintf('\nNeuropil subtraction complete\n')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%% Direction analysis
% load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
% load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(input.tGratingDirectionDeg);
nCells = size(npSub_tc,2);

data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff./2:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
dirs = unique(dir_mat);
nDir = length(dirs);

base_win = [nOff-5:nOff];
resp_win = [nOff+4:nOff+nOn];
tt_dir = (1-nOff:nOn).*(1000/frameRateHz);

data_dfof_dir = zeros(nOn+nOff, nCells, nDir);
h_dir = zeros(nCells,nDir);
p_dir = zeros(nCells,nDir);
dir_resp_avg = zeros(nCells,nDir,2);
dir_resp_mat = [];
dir_list = [];
for iDir = 1:nDir
    ind = find(dir_mat == dirs(iDir));
    data_dfof_dir(:,:,iDir) = mean(data_dfof(:,:,ind),3);
    dir_resp_avg(:,iDir,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind),1),3));
    dir_resp_avg(:,iDir,2) = squeeze(std(mean(data_dfof(resp_win,:,ind),1),[],3))./sqrt(length(ind));
    dir_resp_mat = [dir_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind),1))'];
    dir_list = [dir_list; iDir.*ones(length(ind),1)];
    [h_dir(:,iDir), p_dir(:,iDir)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind),1)), squeeze(mean(data_dfof(base_win,:,ind),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./(nDir-1));
end
h_dir_all = zeros(1,nCells);
p_diranova = zeros(1,nCells); 
k1_hat = nan(1,nCells);
R1_hat = nan(1,nCells);
R2_hat = nan(1,nCells);
u1_hat = nan(1,nCells);
u2_hat = nan(1,nCells);
sse = nan(1,nCells);
R_square = nan(1,nCells);
range = 1:1:360;
figure;
start = 1;
n = 0;
for iCell = 1:nCells
    if find(p_dir(iCell,:)<0.05/16)
        h_dir_all(iCell) = 1;
    elseif length(find(p_dir(iCell,:)<0.05/8))>=2
        h_dir_all(iCell) = 1;
    elseif length(find(p_dir(iCell,:)<0.05/4))>=4
        h_dir_all(iCell) = 1;
    end
    p_diranova(iCell) = anova1(dir_resp_mat(:,iCell),dir_list,'off');
    if start>25
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuningAndFitsCells' num2str(start+((n-1).*25)) '-' num2str(25+((n-1).*25)) '.pdf']))
        figure;
        n = n+1;
        start = 1;
    end
    subplot(5,5,start)
    errorbar(dirs, dir_resp_avg(iCell,:,1), dir_resp_avg(iCell,:,2), 'ok')
    hold on
    if h_dir_all(iCell) || p_diranova(iCell)<0.05
        data = [dir_resp_avg(iCell,:,1) dir_resp_avg(iCell,1,1)];
        theta = [deg2rad(dirs) 2.*pi];
        [b_hat(:,iCell),k1_hat(:,iCell),R1_hat(:,iCell),R2_hat(:,iCell),u1_hat(:,iCell),u2_hat(:,iCell),sse(:,iCell),R_square(:,iCell)] ...
            = miaovonmisesfit_dir(theta,data);
        y_fit = b_hat(:,iCell)+R1_hat(:,iCell).*exp(k1_hat(:,iCell).*(cos(deg2rad(range)-u1_hat(:,iCell))-1))+R2_hat(:,iCell).*exp(k1_hat(:,iCell).*(cos(deg2rad(range)-u1_hat(:,iCell)-pi)-1));
        plot(range,y_fit)
    end
    title(num2str(chop(p_diranova(iCell),2)))
    start = start+1;
end
dirresp_ind = find(h_dir_all);
dirtuned_ind = find(p_diranova<0.05);
[max_val, max_ind_dir] = max(dir_resp_avg(:,:,1),[],2);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'max_ind_dir', 'data_dfof', 'resp_win', 'base_win', 'tt_dir', 'data_dfof_dir', 'h_dir', 'dirresp_ind', 'dirtuned_ind', 'k1_hat', 'R1_hat', 'R2_hat', 'u1_hat', 'u2_hat', 'sse', 'R_square')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'dir_mat', 'dirs', 'nDir', 'nOn', 'nOff','frameRateHz')



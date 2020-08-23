%% get path names
close all;clear all;clc;

ds = 'CrossOriRandPhase_ExptList';
eval(ds)
nexp = length(expt);
for iexp = 5
rc = behavConstsAV;

%%
if ~isempty(expt(iexp).retFolder)
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).retFolder;
coFolder = expt(iexp).coFolder;
time = expt(iexp).retTime;
nrun = length(ImgFolder);
frameRateHz = params.frameRate;

run_str = catRunName(cell2mat(ImgFolder), nrun);
co_run_str = catRunName(cell2mat(coFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

fprintf(['2P imaging Ret analysis \nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
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

% register to cross-ori experiment

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_reg_shifts.mat']))
[out, data_reg] = stackRegister(data,data_avg);
data_reg_avg = mean(data_reg,3);
mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data
% test stability
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

% use cross-ori mask to get TCs

fprintf(['Loading masks from cross-ori runs: ' cell2mat(coFolder) '\n'])

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_mask_cell.mat']))
fprintf('Cell and neuropil masks loaded\n')

nCells = max(mask_cell(:)); % take max label of mask_cell, should circumvent bwlabel
fprintf([num2str(nCells) ' total cells selected\n'])
fprintf('Cell segmentation complete\n')

%neuropil subtraction
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

%% eyetracking
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_pupil.mat']),'rect')

irun =  1;
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
cd(CD);
fn = [ImgFolder{irun} '_000_000_eye.mat'];
data = load(fn);          % should be a '*_eye.mat' file

data = squeeze(data.data);      % the raw images...

data = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);

rad_range = [3 19];
warning off;
A = cell(size(data,3),1);
B = cell(size(data,3),1);
C = cell(size(data,3),1);
D = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
    C{n} = [0];
    D{n} = [0];
end
eye = struct('Centroid',A,'Area',B,'Val',C,'SNR',D);
radii = [];
for n = 1:size(data,3)
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
    [val,idx] = max(metric);          % pick the circle with best score
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
        eye(n).Val = NaN;
        eye(n).SNR = NaN;
    else
        t = double(data(:,:,n));
        vector_of_y_values = (1:size(data,1)) - center(idx,2);
        vector_of_x_values = (1:size(data,2)) - center(idx,1);
        [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
        idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
        idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
        snr = mean(t(idx1))./mean(t(idx2));
        eye(n).SNR = snr;
        eye(n).Val = val;
        eye(n).Centroid = center(idx,:);
        eye(n).Area = pi*radii(idx)^2;
    end
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');
Val = double(cell2mat({eye.Val}'));
SNR = double(cell2mat({eye.SNR}'));
Eye_data = data;
nanframes(1,iexp) = length(find(isnan(Area)));

% no measurement frames
figure; 
subplot(2,2,1)
hist(sqrt(Area./pi));
xlabel('radius')
subplot(2,2,2)
hist(SNR);
xlabel('SNR')
subplot(2,2,3)
hist(Val);
xlabel('Metric')

x1 = find(isnan(Area));
x2 = find(~isnan(Area));
x3 = find(Val<0.26 & SNR<1.9);

x = unique([x1; x3]);
if length(x)>25
    minx = 25;
else
    minx = length(x);
end

frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
suptitle('No pupil detected')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil2.pdf']),'-dpdf','-fillpage');

x = setdiff(x2,x3);
frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
suptitle('Pupil detected')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');
        
%align eyetracking to trials 
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(celleqel2mat_padded(input.tGratingDirectionDeg));
Rad_temp = sqrt(Area./pi);
Centroid_temp = Centroid;
Rad_temp(unique([x1; x3]),:) =nan(length(unique([x1; x3])),1);
Centroid_temp(unique([x1; x3]),:) = nan(length(unique([x1; x3])),2);
sz = size(Eye_data);
rad_mat_start = zeros(nOn+nOff, ntrials);
centroid_mat_start = zeros(nOn+nOff,2, ntrials);
eye_mat_start = zeros(sz(1), sz(2),nOn+nOff, ntrials);

nframes = size(Rad_temp,1);
for itrial = 1:ntrials
    crange_on = nOff+((nOn+nOff).*(itrial-1)):(nOn+nOff).*itrial;
    if sum(isnan(Rad_temp(crange_on,1)),1)>0
        if sum(isnan(Rad_temp(crange_on,1)),1)./length(crange_on)> 0.25
            Rad_temp(crange_on,1) = NaN(length(crange_on),1);
            Centroid_temp(crange_on,:) = NaN(length(crange_on),2);
        else
            nanind = intersect(crange_on,find(isnan(Rad_temp)));
            dataind = intersect(crange_on,find(~isnan(Rad_temp)));
            for inan = 1:length(nanind)
                gap = min(abs(nanind(inan)-dataind),[],1);
                good_ind_stim = find(abs(nanind(inan)-dataind) == gap);
                Rad_temp(nanind(inan),1) = mean(Rad_temp(dataind(good_ind_stim),1),1);
                Centroid_temp(nanind(inan),:) = mean(Centroid(dataind(good_ind_stim),:),1);
            end
        end
    end
    crange_all = 1+((nOn+nOff).*(itrial-1)):(nOn+nOff).*itrial;
    rad_mat_start(:,itrial) = Rad_temp(crange_all,:);
    centroid_mat_start(:,:,itrial) = Centroid_temp(crange_all,:);
    eye_mat_start(:,:,:,itrial) = Eye_data(:,:,crange_all);
end
calib = 1/26.6; %mm per pixel
rad_mat_calib = bsxfun(@times, rad_mat_start, calib);
centroid_mat_calib = bsxfun(@times,centroid_mat_start,calib);
t = mean(centroid_mat_calib(nOff+1:end,:,:),1);
rad_base = mean(rad_mat_calib(1:nOff,:),1);
rad_stim = mean(rad_mat_calib(nOff+1:end,:),1);
centroid_base = squeeze(mean(centroid_mat_calib(1:nOff,:,:),1))./0.025;
centroid_stim = squeeze(mean(centroid_mat_calib(nOff+1:end,:,:),1))./0.025;
ind = find(~isnan(centroid_stim(1,:)));
centroid_med = findMaxNeighbors(centroid_stim(:,ind),2);
figure; subplot(2,1,1)
scatter(centroid_stim(1,:),centroid_stim(2,:), [], rad_stim); colorbar
hold on;
scatter(centroid_med(1),centroid_med(2),'og')
centroid_dist = sqrt((centroid_stim(1,:)-centroid_med(1)).^2 + (centroid_stim(2,:)-centroid_med(2)).^2);
title('Color- radius')
xlabel('x-pos')
ylabel('y-pos')
subplot(2,1,2)
hist(centroid_dist,0:0.5:60)
suptitle([num2str(sum(~isnan(centroid_dist))) '/' num2str(ntrials) ' measurable trials; ' num2str(length(find(centroid_dist<2))) ' within 2deg'])
xlabel('Centroid distance from median')
movegui('center')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilPosDist.pdf']),'-dpdf','-fillpage');

[n edges bin] = histcounts(centroid_dist,[0:2:30]);
    
i = find(n);
[n1 n2] = subplotn(length(i)); 
figure;
for ii = 1:length(i)
    subplot(n1,n2,ii)
    ind = find(bin== i(ii),1);
    if ii == 1
        ind_i = ind;
    end
    imagesc(mean(eye_mat_start(:,:,nOff+1:end,ind),3))
    hold on
    plot(squeeze(nanmean(centroid_mat_start(nOff+1:end,1,ind),1)), squeeze(nanmean(centroid_mat_start(nOff+1:end,2,ind),1)),'or')
    plot(squeeze(nanmean(centroid_mat_start(nOff+1:end,1,ind_i),1)), squeeze(nanmean(centroid_mat_start(nOff+1:end,2,ind_i),1)),'ok')
    title(num2str(edges(ii)))
end
suptitle('Example eye image by distance from median')
movegui('center')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilImgByDist.pdf']),'-dpdf','-fillpage');
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'rad_mat_start','centroid_mat_start', 'rad_base','rad_stim','centroid_base', 'centroid_stim', 'centroid_dist', 'centroid_med' );

%% Retiotopy analysis
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = length(input.tGratingDirectionDeg);
nCells = size(npSub_tc,2);

data_trial = permute(reshape(npSub_tc,[nOn+nOff ntrials nCells]),[1 3 2]);
data_f = mean(data_trial(nOff./2:nOff,:,:),1);
data_dfof = (data_trial-data_f)./data_f;

el_mat = celleqel2mat_padded(input.tGratingElevationDeg);
els = unique(el_mat);
nEl = length(els);
az_mat = celleqel2mat_padded(input.tGratingAzimuthDeg);
azs = unique(az_mat);
nAz = length(azs);
pha_mat = celleqel2mat_padded(input.tGratingStartingPhaseDeg);
phas = unique(pha_mat);
nPha = length(phas);

base_win = [nOff-5:nOff];
resp_win = [nOff+4:nOff+nOn];
tt_dir = (1-nOff:nOn).*(1000/frameRateHz);

data_dfof_stim = zeros(nOn+nOff, nCells, nEl, nAz, nPha);
h_stim = zeros(nCells, nEl, nAz, nPha);
p_stim = zeros(nCells, nEl, nAz, nPha);
trial_n = zeros(nEl, nAz, nPha,2);
stim_resp_avg = zeros(nCells, nEl, nAz, nPha,2);
stim_resp_mat = [];
stim_list = [];
for iEl = 1:nEl
    ind_el = find(el_mat == els(iEl));
    for iAz = 1:nAz
        ind_az = find(az_mat == azs(iAz));
        for iP = 1:nPha
            ind_p = find(pha_mat == phas(iP));
            ind_stim = intersect(ind_p, intersect(ind_az,ind_el));
            ind_use = intersect(ind_stim, find(centroid_dist<2));
            trial_n(iEl,iAz,iP,1) = length(ind_stim);
            trial_n(iEl,iAz,iP,2) = length(ind_use);
            data_dfof_stim(:,:,iEl,iAz,iP) = mean(data_dfof(:,:,ind_use),3);
            stim_resp_avg(:,iEl,iAz,iP,1) = squeeze(mean(mean(data_dfof(resp_win,:,ind_use),1),3));
            stim_resp_avg(:,iEl,iAz,iP,2) = squeeze(std(mean(data_dfof(resp_win,:,ind_use),1),[],3))./sqrt(length(ind_use));
            stim_resp_mat = [stim_resp_mat; squeeze(mean(data_dfof(resp_win,:,ind_use),1))'];
            stim_list = [stim_list; iEl.*ones(length(ind_use),1) iAz.*ones(length(ind_use),1) iP.*ones(length(ind_use),1)];
            [h_stim(:,iEl,iAz,iP), p_stim(:,iEl,iAz,iP)] = ttest(squeeze(mean(data_dfof(resp_win,:,ind_use),1)), squeeze(mean(data_dfof(base_win,:,ind_use),1)),'dim', 2, 'tail', 'right', 'alpha', 0.05./((nAz.*nEl)-1));
        end
    end
end

h_stim_all = zeros(1,nCells);
p_stimanova = zeros(3,nCells); 
rf_all = zeros(nEl,nAz,3,nCells);
stim_list_cell = cell(1,3);
for i = 1:3
    stim_list_cell{1,i} = stim_list(:,i);
end
figure;
movegui('center')
start = 1;
n = 1;
for iCell = 1:nCells
    for iP = 1:nPha
        if find(p_stim(iCell,:,:,iP)<0.05/((nAz.*nEl)-1))
            h_stim_all(iCell) = 1;
        elseif length(find(p_stim(iCell,:)<0.05/((nAz.*nEl)/2)))>=2
            h_stim_all(iCell) = 1;
        elseif length(find(p_stim(iCell,:)<0.05/((nAz.*nEl)/4)))>=4
            h_stim_all(iCell) = 1;
        end
    end
    p_stimanova(:,iCell) = anovan(stim_resp_mat(:,iCell),stim_list_cell,'display','off');
    if start>25
        suptitle([mouse ' ' date ' Ret'])
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_retCells' num2str(start+((n-1).*25)) '-' num2str(25+((n-1).*25)) '.pdf']))
        figure;
        movegui('center')
        n = n+1;
        start = 1;
    end
    subplot(5,5,start)
    rgb_temp = zeros(nEl,nAz,3);
    rgb_temp(:,:,1) = squeeze(stim_resp_avg(iCell,:,:,1,1));
    rgb_temp(:,:,3) = squeeze(stim_resp_avg(iCell,:,:,2,1));
    rgb_temp_norm = rgb_temp./max(rgb_temp(:));
    imagesc(rgb_temp_norm)
    rf_all(:,:,:,iCell) = rgb_temp;
    start = start+1;
    title(['Cell ' num2str(iCell) '- ' num2str(h_stim_all(iCell))])
end
suptitle([mouse ' ' date ' Ret'])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_retCells' num2str(start+((n-1).*25)) '-' num2str(25+((n-1).*25)) '.pdf']))
stimresp_ind = find(h_stim_all);
h_stimanova = p_stimanova<0.05;
stimtuned_ind = find(sum(h_stimanova,1));

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']), 'rf_all', 'data_dfof', 'resp_win', 'base_win', 'data_dfof_stim', 'h_stim', 'stimresp_ind', 'stimtuned_ind', 'trial_n')
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']), 'el_mat', 'els', 'nEl', 'az_mat', 'azs', 'nAz', 'pha_mat','phas', 'nPha', 'nOn', 'nOff','frameRateHz')

end
end

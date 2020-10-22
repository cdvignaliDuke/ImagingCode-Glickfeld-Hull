mouse = 'i726';
date = '170712';
time_mat = strvcat('1318','1348');
pn = '\\crash.dhe.duke.edu\lindsey\Data\Widefield_images';
suf = '_MMStack_Pos0.ome.tif';
nrun = size(time_mat,1);

%%
data = [];
for irun = 1:nrun
    data_temp = readtiff(fullfile(pn, [date '_' mouse], [date '_' mouse '_' num2str(irun)], [date '_' mouse '_' num2str(irun) suf]),'single');
    data = cat(3,data,data_temp);
    load(['\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time_mat(irun,:) '.mat'])
    temp(irun) = input;
end
clear data_temp
input = concatenateDataBlocks(temp);

[out, reg] = stackRegister(data, data(:,:,100));

nOn = input.nScansOn;
nOff = input.nScansOff;
speed = celleqel2mat_padded(input.tGratingTemporalFreqCPS)./celleqel2mat_padded(input.tGratingSpatialFreqCPD);
ntrials = length(speed);
speeds = unique(speed);
nspeed = length(speeds);

sz = size(reg);

data_tr = reshape(reg,[sz(1) sz(2) nOn+nOff ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_dfof = bsxfun(@rdivide, bsxfun(@minus, data_tr, data_f), data_f);

figure
for is = 1:nspeed
    subplot(2,2,is)
    ind = find(speed == speeds(is));
    data_avg = mean(data_dfof(:,:,:,ind),4);
    imagesc(mean(data_avg(:,:,nOff:nOff+nOn/2),3));
    colormap gray
    %clim([-0.01 0.008])
    colorbar
    axis square
    axis off
    title(num2str(speeds(is)))
end
subplot(2,2,4)
data_avg = mean(data_dfof(:,:,:,:),4);
imagesc(mean(data_avg(:,:,nOff:nOff+nOn/2),3));
title('All')
colormap gray
%clim([-0.01 0.0008])
colorbar
axis square
axis off
suptitle([mouse ' ' date])
if ~exist(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Widefield_imaging',mouse))
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Widefield_imaging',mouse))
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\Widefield_imaging',mouse,[date '_' mouse '_speedFOVs.pdf']), '-dpdf', '-bestfit')

date = '160622';
ImgFolder = strvcat('001');
time = strvcat('1303');
mouse = 'i919';
nrun = size(ImgFolder,1);
frame_rate = 15;

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

data = [];
for irun = 1:nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end
clear data_temp

data_avg = mean(data(:,:,51:100),3);
[out, data_reg] = stackRegister(data,data_avg);
mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\liniposdsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
clear data

for irun = 1:nrun
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
end
input = concatenateDataBlocks(temp);

nOn = input.nScansOn;
nOff = input.nScansOff;
%ntrials = size(input.tGratingDirectionDeg,2);
ntrials = 45;

sz = size(data_reg);
data_mat = zeros(sz(1), sz(2), nOn+nOff, ntrials);
for itrial = 1:ntrials
    data_mat(:,:,:,itrial) = data_reg(:,:,1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)));
end

data_f = mean(data_mat(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, data_mat, data_f);
data_dfof = bsxfun(@rdivide, data_df, data_f);
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof = imfilter(data_dfof, myfilter);

clear data_mat data_f data_df

El = celleqel2mat_padded(input.tGratingElevationDeg(1:ntrials));
Az = celleqel2mat_padded(input.tGratingAzimuthDeg(1:ntrials));
Els = unique(El);
Azs = unique(Az);
nPos = length(Els).*length(Azs);
pos = zeros(nPos,2);
ipos = 1;
figure;
data_pos = zeros(sz(1), sz(2), nOn+nOff, nPos);
for iEl = 1:length(Els)
    ind_E = find(El == Els(iEl));
    for iAz = 1:length(Azs)
        ind_A = find(Az == Azs(iAz));
        ind_pos = intersect(ind_A, ind_E);
        data_pos(:,:,:,ipos) = mean(data_dfof(:,:,:,ind_pos),4);
        subplot(length(Els), length(Azs), ipos)
        imagesc(mean(data_pos(:,:,nOff+1:nOff+nOn),3))
        title(['Az: ' num2str(Azs(iAz)) '; El: ' num2str(Els(iEl))])
        pos(ipos,:) = [Azs(iAz) Els(iEl)];
        ipos = ipos+1;
    end
end

figure;
for ipos = 1:nPos
     subplot(length(Els), length(Azs), ipos)
     plot(squeeze(mean(mean(data_pos(:,:,:,ipos),1),2)))
     ylim([-0.05 0.1])
     title(['Az: ' num2str(pos(ipos,1)) '; El: ' num2str(pos(ipos,2))])
end

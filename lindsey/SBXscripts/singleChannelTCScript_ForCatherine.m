%% get path names
date = '200118';
ImgFolder = strvcat('002');
time = strvcat('1508');
mouse = 'i1312';
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
%change these as needed
lg_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
%% load and register
data = [];
clear temp
trial_n = [];
for irun = 1:nrun
    %change to image directory and load imaging info
    CD = fullfile(lg_fn, mouse, date, ImgFolder(irun,:));
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    
    %load mworks data
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);
    temp(irun) = input; %needed in case there are multiple runs
    
    %load image data
    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    %this deals with a mismatch in the number of frames/trials in the
    %imaging and mworks data
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
        
    data_temp = squeeze(data_temp); %removes singleton dimension if only 1 PMT
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
nframes = 2000;
nep = floor(size(data,3)./nframes);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end

%% Register data

data_avg = mean(data(:,:,6001:6500),3); %change this according to best image, closest to middle of run from above

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str])) %if rerunning analysis using old registration shifts
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg(:,:,1:10000),3);
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end
figure; imagesq(data_reg_avg); truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.mat']))

%% find activated cells
    
%max by trial type    

nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2);
if nOn>29
    sz = size(data_reg);
    data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(data_tr), data_f); 
    data_dfof = bsxfun(@rdivide,data_df, data_f); 
    clear data_f data_df data_tr
else
    sz = size(data_reg);
    data_tr = zeros(sz(1),sz(2), 100, ntrials-1);
    for itrial = 1:ntrials-1
        data_tr(:,:,:,itrial) = data_reg(:,:,((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)));
    end
    data_f = mean(data_tr(:,:,1:50,:),3);
    data_df = bsxfun(@minus, double(data_tr), data_f); 
    data_dfof = bsxfun(@rdivide,data_df, data_f); 
    clear data_f data_df data_tr
end

if input.doDirStim
    Dir = cell2mat(input.tGratingDirectionDeg);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
    [n n2] = subplotn(nDirs);
    figure;
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir))
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    
    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i));
        clim([0 max(data_dfof_avg_all(:))])
        title(num2str(Dirs(i)))
        colormap(gray)
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3);
        end
    end
    for i = 1:nStim/2
        subplot(2,2,i)
        imagesc(data_dfof_avg_ori(:,:,i));
        clim([0 max(data_dfof_avg_ori(:))])
    end
    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    
elseif input.doRetStim
    Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
    El = celleqel2mat_padded(input.tGratingElevationDeg);
    Azs = unique(Az);
    Els = unique(El);
    if min(Els,[],2)<0
        Els = fliplr(Els);
    end
    nStim = length(Azs).*length(Els);
    Stims = [];
    data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
    start = 1;
    for iEl = 1:length(Els)
        ind1 = find(El == Els(iEl));
        for iAz = 1:length(Azs)
            Stims = [Stims; Els(iEl) Azs(iAz)];
            ind2 = find(Az == Azs(iAz));
            ind = intersect(ind1,ind2);
            data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
            start = start +1;
        end
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

    figure; 
    for i = 1:nStim; 
        subplot(length(Els),length(Azs),i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        colormap(gray)
        title(num2str(Stims(i,:)))
        clim([0 max(data_dfof_avg_all(:))])
    end
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')

end

 save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')


%% cell segmentation 

mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_all >= 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_2 = bwlabel(bwout);
    mask_all = mask_all+mask_2;
    close all
end
mask_cell = bwlabel(mask_all);

mask_data_temp = data_dfof_max;
mask_data_temp(find(mask_cell >= 1)) = 0;
figure; imagesc(mask_data_temp)


bwout = imCellEditInteractive(data_dfof_max);
mask_cell1 = bwlabel(bwout);


figure; 
[n n2] = subplotn(nStim);
for i = 1:nStim; 
    subplot(n,n2,i); 
    shade_img = imShade(data_dfof_avg_all(:,:,i), mask_all);
    imagesc(shade_img)
    if input.doSizeStim
    title([num2str(szs(i)) ' deg'])
    elseif input.doRetStim
        title([num2str(Stims(i,:))])
    end
    clim([0 max(data_dfof_avg_all(:))])
    colormap(gray)
end
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')


mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 


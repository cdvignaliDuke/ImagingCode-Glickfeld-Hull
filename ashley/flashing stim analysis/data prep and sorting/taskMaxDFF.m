fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,dirFolder);
load(fullfile(fnout,'reg_img'))
xpix = size(data_avg,2);
ypix = size(data_avg,1);

for irun = 1:expt(iexp).nrun;
    % load data
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    
    if iexp == 2 & irun == 3
        continue
    else
tic
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName);    
toc
    fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,runFolder);
    end
    
    if ~exist(fnout,'dir')
        mkdir(fnout)
    end

    data_sub = data-min(min(min(data,[],1),[],2),[],3);
    data = data_sub;
    clear data_sub

    % need only trial-related frames: trial start, last trial, target
    if iscell(input.nFramesOn)
        on = unique(cell2mat(input.nFramesOn));
        off = unique(cell2mat(input.nFramesOff));
    else
        on = input.nFramesOn;
        off = input.nFramesOff;        
    end
    
    frame_rate = expt(iexp).frame_rate;
    nframes_1s = ceil(frame_rate./down)*down;
    nframes_2stim = ceil((on+off)*2./down)*down;

    tr_start_all = double(cell2mat(input.cLeverDown));
    tr_L = cell2mat(input.tCyclesOn);
    tr_tar_all = cell2mat_padded(input.cTargetOn);

    tr_bl_ind = linspaceNDim(tr_start_all-nframes_1s+1,tr_start_all,nframes_1s);

    tr_start = double(tr_start_all(tr_L >= 3));
    tr_start_ind = linspaceNDim(tr_start+1,tr_start+nframes_2stim,nframes_2stim);

    tr_long = double(tr_start_all(tr_L >= 6));
    tr_long_ind = linspaceNDim(tr_long+(2*nframes_2stim)+1,tr_long+(3*nframes_2stim),nframes_2stim);

    tr_tar_on = ismember(input.trialOutcomeCell,{'success';'ignore'});
    tr_tar = double(tr_tar_all(tr_tar_on));
    tr_tar_ind = linspaceNDim(tr_tar+1,tr_tar+nframes_1s,nframes_1s);

    tr_bl_frames = data(:,:,tr_bl_ind(:));
    tr_start_frames = data(:,:,tr_start_ind(:));
    tr_long_frames = data(:,:,tr_long_ind(:));
    tr_tar_frames = data(:,:,tr_tar_ind(:));
 
    
    tr_bl_down = stackGroupProject(tr_bl_frames,down);
    tr_start_down = stackGroupProject(tr_start_frames,down);
    tr_long_down = stackGroupProject(tr_long_frames,down);
    tr_tar_down = stackGroupProject(tr_tar_frames,down);


    % register needed frames using tuning outs
    [out_bl, tr_bl_reg] = stackRegister(tr_bl_down,data_avg);
    [out_start, tr_start_reg] = stackRegister(tr_start_down,data_avg);
    [out_long, tr_long_reg] = stackRegister(tr_long_down,data_avg);
    [out_tar, tr_tar_reg] = stackRegister(tr_tar_down,data_avg);

    % get dF/F for each
    tr_bl_reg = double(reshape(tr_bl_reg,ypix,xpix,nframes_1s./down,length(tr_start_all)));
    tr_start_reg = double(reshape(tr_start_reg,ypix,xpix,nframes_2stim./down,length(tr_start)));
    tr_long_reg = double(reshape(tr_long_reg,ypix,xpix,nframes_2stim./down,length(tr_long)));
    tr_tar_reg = double(reshape(tr_tar_reg,ypix,xpix,nframes_1s./down,length(tr_tar)));

    F = mean(tr_bl_reg,3);

    if irun == 1
        dF = bsxfun(@minus,tr_start_reg,F(:,:,:,tr_L >= 3));
        dFF_start = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 3));

        dF = bsxfun(@minus,tr_long_reg,F(:,:,:,tr_L >= 6));
        dFF_long = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 6));

        dF = bsxfun(@minus,tr_tar_reg,F(:,:,:,tr_tar_on));
        dFF_tar = bsxfun(@rdivide,dF,F(:,:,:,tr_tar_on));    
    else
        dF = bsxfun(@minus,tr_start_reg,F(:,:,:,tr_L >= 3));
        dFF_start_temp = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 3));
        dFF_start = cat(4,dFF_start,dFF_start_temp);

        dF = bsxfun(@minus,tr_long_reg,F(:,:,:,tr_L >= 6));
        dFF_long_temp = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 6));
        dFF_long = cat(4,dFF_long,dFF_long_temp);

        dF = bsxfun(@minus,tr_tar_reg,F(:,:,:,tr_tar_on));
        dFF_tar_temp = bsxfun(@rdivide,dF,F(:,:,:,tr_tar_on));
        dFF_tar = cat(4,dFF_tar,dFF_tar_temp);
    end

    clear data
end

dFF_start = reshape(dFF_start,ypix,xpix,[]);
dFF_long = reshape(dFF_long,ypix,xpix,[]);
dFF_tar = reshape(dFF_tar,ypix,xpix,[]);

% dFF_start(:,xcrop,:) = 0;
% dFF_start(ycrop,:,:) = 0;
% dFF_long(:,xcrop,:) = 0;
% dFF_long(ycrop,:,:) = 0;
% dFF_tar(:,xcrop,:) = 0;
% dFF_tar(ycrop,:,:) = 0;

start_max = max(dFF_start,[],3);
long_max = max(dFF_long,[],3);
tar_max = max(dFF_tar,[],3);

fnsave = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate);

figure;
setFigParams4Print('portrait')
suptitle([mouse '-' expDate ', max dF/F'])
colormap gray
subplot(3,1,1)
imagesc(start_max)
title('trial start')
subplot(3,1,2)
imagesc(long_max)
title('trial mid to end')
subplot(3,1,3)
imagesc(tar_max)
title('trial target')
print(fullfile(fnsave,'bx_max'),'-dpdf','-fillpage');

figure;
colormap gray
subplot(3,1,1)
imagesc(mean(dFF_start,3))
title('trial start')
subplot(3,1,2)
imagesc(mean(dFF_long,3))
title('trial mid to end')
subplot(3,1,3)
imagesc(mean(dFF_tar,3))
title('trial target')

%% save max images,reg outs
save(fullfile(fnsave,'bx_max_images.mat'),'start_max','long_max','tar_max');
% save(fullfile(fnout,'bx_maxImg_regOuts.mat'),'

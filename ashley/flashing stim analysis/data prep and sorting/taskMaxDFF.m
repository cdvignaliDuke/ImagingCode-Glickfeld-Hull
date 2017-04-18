
xpix = size(data_bx_reg,2);
ypix = size(data_bx_reg,1);


    % need only trial-related frames: trial start, last trial, target
    if iscell(input_bx.nFramesOn)
        on = unique(cell2mat(input_bx.nFramesOn));
        off = unique(cell2mat(input_bx.nFramesOff));
    else
        on = input_bx.nFramesOn;
        off = input_bx.nFramesOff;        
    end
    
    frame_rate = expt(iexp).frame_rate;
    nframes_1s = frame_rate;
    nframes_2stim = (on+off)*2;

    tr_start_all = double(cell2mat(input_bx.cLeverDown));
    tr_L = cell2mat(input_bx.tCyclesOn);
    tr_tar_all = cell2mat_padded(input_bx.cTargetOn);

    tr_bl_ind = linspaceNDim(tr_start_all-nframes_1s+1,tr_start_all,nframes_1s);

    tr_start = double(tr_start_all(tr_L >= 3));
    tr_start_ind = linspaceNDim(tr_start+1,tr_start+nframes_2stim,nframes_2stim);

    tr_long = double(tr_start_all(tr_L >= 6));
    tr_long_ind = linspaceNDim(tr_long+(2*nframes_2stim)+1,tr_long+(3*nframes_2stim),nframes_2stim);

    tr_tar_on = ismember(input_bx.trialOutcomeCell,{'success';'ignore'});
    tr_tar = double(tr_tar_all(tr_tar_on));
    tr_tar_ind = linspaceNDim(tr_tar+1,tr_tar+nframes_1s,nframes_1s);

    tr_bl_frames = data_bx_reg(:,:,tr_bl_ind(:));
    tr_start_frames = data_bx_reg(:,:,tr_start_ind(:));
    tr_long_frames = data_bx_reg(:,:,tr_long_ind(:));
    tr_tar_frames = data_bx_reg(:,:,tr_tar_ind(:));
    clear data_bx_reg
 

%     % register needed frames using tuning outs
%     [out_bl, tr_bl_reg] = stackRegister(tr_bl_down,data_corr);
%     [out_start, tr_start_reg] = stackRegister(tr_start_down,data_corr);
%     [out_long, tr_long_reg] = stackRegister(tr_long_down,data_corr);
%     [out_tar, tr_tar_reg] = stackRegister(tr_tar_down,data_corr);

    % get dF/F for each
    tr_bl_frames = double(reshape(tr_bl_frames,ypix,xpix,nframes_1s,length(tr_start_all)));
    tr_start_frames = double(reshape(tr_start_frames,ypix,xpix,nframes_2stim,length(tr_start)));
    tr_long_frames = double(reshape(tr_long_frames,ypix,xpix,nframes_2stim,length(tr_long)));
    tr_tar_frames = double(reshape(tr_tar_frames,ypix,xpix,nframes_1s,length(tr_tar)));

    F = mean(tr_bl_frames,3);

    dF = bsxfun(@minus,tr_start_frames,F(:,:,:,tr_L >= 3));
    dFF_start = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 3));
    clear tr_start_frames

    dF = bsxfun(@minus,tr_long_frames,F(:,:,:,tr_L >= 6));
    dFF_long = bsxfun(@rdivide,dF,F(:,:,:,tr_L >= 6));
    clear tr_long_frames

    dF = bsxfun(@minus,tr_tar_frames,F(:,:,:,tr_tar_on));
    dFF_tar = bsxfun(@rdivide,dF,F(:,:,:,tr_tar_on));  
    clear tr_tar_frames
    
    clear dF 

dFF_start = reshape(dFF_start,ypix,xpix,[]);
dFF_long = reshape(dFF_long,ypix,xpix,[]);
dFF_tar = reshape(dFF_tar,ypix,xpix,[]);

start_max = max(dFF_start,[],3);
long_max = max(dFF_long,[],3);
tar_max = max(dFF_tar,[],3);

figure;
setFigParams4Print('portrait')
suptitle([mouse '-' expDate ', max dF/F behavior'])
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
if ~exist(fullfile(fnbx,'max images'),'dir')
    mkdir(fullfile(fnbx,'max images'))
end
print(fullfile(fnbx,'max images',[mouse '_' expDate '_bx.pdf']),'-dpdf','-fillpage');
writetiff(max(cat(3,start_max,long_max,tar_max),[],3),fullfile(fnbx,'max images',[mouse '_' expDate '_bx']))
%% save max images,reg outs
save(fullfile(fnbx,'bx_max_images.mat'),'start_max','long_max','tar_max');

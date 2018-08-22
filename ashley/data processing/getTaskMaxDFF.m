xpix = size(data_bx_reg,2);
ypix = size(data_bx_reg,1);

nRuns = expt(iexp).nrun;

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
    
    nTrialsEaRun = input_bx.trialsSinceReset;
    cmlvTrialSum = cumsum(nTrialsEaRun);

    tr_start_all = double(cell2mat(input_bx.cFirstStim));
    tr_L = cell2mat(input_bx.tCyclesOn);
    tr_tar_all = cell2mat_padded(input_bx.cTargetOn);
    
    tr_start_offset = cell(1,nRuns);
    tr_L_offset = cell(1,nRuns);
    tr_tar_offset = cell(1,nRuns);
    for irun = 1:nRuns
        if irun == 1
            offset = 0;
            ind1 = 1:cmlvTrialSum(irun);
            tr_tar_thisrun = tr_tar_all(ind1);
            ind2 = tr_tar_thisrun < nFrPerRun(irun);
            tr_start_offset{irun} = tr_start_all(ind2);
            tr_L_offset{irun} = tr_L(ind2);
            tr_tar_offset{irun} = tr_tar_all(ind2);
        else
            offset = nFrPerRun(irun-1)+offset;
            ind1 = cmlvTrialSum(irun-1)+1:cmlvTrialSum(irun);
            tr_tar_thisrun = tr_tar_all(ind1);
            tr_start_thisrun = tr_start_all(ind1);
            tr_L_thisrun = tr_L(ind1);
            ind2 = tr_tar_thisrun < nFrPerRun(irun);
            tr_tar_offset{irun} = tr_tar_thisrun(ind2) + offset;
            tr_start_offset{irun} = tr_start_thisrun(ind2) + offset;
            tr_L_offset{irun} = tr_L_thisrun(ind2);
        end
    end
    tr_start_all = cell2mat(tr_start_offset);
    tr_tar_all = cell2mat(tr_tar_offset);
    tr_L = cell2mat(tr_L_offset);

%     if strcmp(expDate,'150508')
%         ind = find(tr_tar_all > size(data_bx_reg,3),1,'first');
%         tr_start_all = tr_start_all(1:(ind-2));
%         tr_L = tr_L(1:(ind-2));
%         tr_tar_all = tr_tar_all(1:(ind-2));
%         tr_tar_on = ismember(input_bx.trialOutcomeCell,{'success';'ignore'});
%         tr_tar_on = tr_tar_on(1:(ind-2));
%     else        
        tr_tar_on = ismember(input_bx.trialOutcomeCell,{'success';'ignore'});
%     end

tr_start = tr_L >= 3;
tr_long = tr_L >= 6;

nt = length(tr_tar_on);
tr_bl_frames = nan(ypix,xpix,nframes_1s,nt);
tr_start_frames = nan(ypix,xpix,nframes_2stim,sum(tr_start));
ind1 = 0;
tr_long_frames = nan(ypix,xpix,nframes_2stim,sum(tr_long));
ind2 = 0;
tr_tar_frames = nan(ypix,xpix,nframes_1s,sum(tr_tar_on));
ind3 = 0;
ind4 = 0;
for i = 1:nt
    tr_bl_frames(:,:,:,i) = data_bx_reg(:,:,...
        (tr_start_all(i)-nframes_1s+1):tr_start_all(i));
    if tr_start(i)
        ind1 = ind1+1;
        tr_start_frames(:,:,:,ind1) = data_bx_reg(:,:,...
        double((tr_start_all(i)+1):(tr_start_all(i)+nframes_2stim)));
    end
    if tr_long(i)
        ind2 = ind2+1;
        tr_long_frames(:,:,:,ind2) = data_bx_reg(:,:,...
        (tr_start_all(i)+(2*nframes_2stim)+1):(tr_start_all(i)+(3*nframes_2stim)));
    end
    if tr_tar_on(i) && (tr_tar_all(i)+nframes_1s) < nfr_bx
        ind3 = ind3+1;
        tr_tar_frames(:,:,:,ind3) = data_bx_reg(:,:,...
        (tr_tar_all(i)+1):(tr_tar_all(i)+nframes_1s));
    elseif (tr_tar_all(i)+nframes_1s) < nfr_bx
        tr_tar_frames = tr_tar_frames(:,:,:,1:ind3);
        if ind4 == 0
            tr_tar_on = tr_tar_on(1:nt);
            ind4 = 1;
        end
    end
end


%     tr_bl_ind = linspaceNDim(tr_start_all-nframes_1s+1,tr_start_all,nframes_1s);
% 
%     tr_start = double(tr_start_all(tr_L >= 3));
%     tr_start_ind = linspaceNDim(tr_start+1,tr_start+nframes_2stim,nframes_2stim);
% 
%     tr_long = double(tr_start_all(tr_L >= 6));
%     tr_long_ind = linspaceNDim(tr_long+(2*nframes_2stim)+1,tr_long+(3*nframes_2stim),nframes_2stim);
% 
%     tr_tar = double(tr_tar_all(tr_tar_on));
%     tr_tar_ind = linspaceNDim(tr_tar+1,tr_tar+nframes_1s,nframes_1s);
%     
%     if tr_tar_ind(end) > size(data_bx_reg,3)
%         ind = find(tr_tar+nframes_1s > size(data_bx_reg,3),1,'first');
%         nt = ind-1;
%         tr_bl_ind = tr_bl_ind(1:nt,:);
%         tr_start_ind = tr_start_ind(1:nt,:);
%         tr_tar_ind = tr_tar_ind(1:nt,:);
%         tr_tar_on = tr_tar_on(1:nt);
%         tr_L = tr_L(1:nt);
%     end
%     
%     if tr_long_ind(end) > size(data_bx_reg,3)
%         ind = find(tr_long+nframes_1s > size(data_bx_reg,3),1,'first');
%         nt = ind-1;
%         tr_long_ind = tr_long_ind(1:nt,:);
%     end
%         
%     tr_bl_frames = data_bx_reg(:,:,tr_bl_ind(:));
%     tr_start_frames = data_bx_reg(:,:,tr_start_ind(:));
%     tr_long_frames = data_bx_reg(:,:,tr_long_ind(:));
%     tr_tar_frames = data_bx_reg(:,:,tr_tar_ind(:));

%     % get dF/F for each
%     tr_bl_frames = double(reshape(tr_bl_frames,ypix,xpix,nframes_1s,...
%         size(tr_bl_ind,1)));
%     tr_start_frames = double(reshape(tr_start_frames,ypix,xpix,nframes_2stim,...
%         size(tr_start_ind,1)));
%     tr_long_frames = double(reshape(tr_long_frames,ypix,xpix,nframes_2stim,...
%         size(tr_long_ind,1)));
%     tr_tar_frames = double(reshape(tr_tar_frames,ypix,xpix,nframes_1s,...
%         size(tr_tar_ind,1)));

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
    
    clear dF F


dFF_start_meanAllTrials = squeeze(mean(dFF_start,3));
dFF_long_meanAllTrials = squeeze(mean(dFF_long,3));
dFF_tar_meanAllTrials = squeeze(mean(dFF_tar,3));

start_max = max(dFF_start_meanAllTrials,[],3);
long_max = max(dFF_long_meanAllTrials,[],3);
tar_max = max(dFF_tar_meanAllTrials,[],3);

% remove mis-aligned frames
motionCutoff = 0.2;
maxDFF = max(cat(3,start_max,long_max,tar_max),[],3);
bwout = imCellEditInteractive(maxDFF);
mask = bwlabel(bwout);
nROI = length(unique(mask))-1;

start_tc = reshape(stackGetTimeCourses(reshape(dFF_start,ypix,xpix,[]),mask),...
    nframes_2stim,sum(tr_start),nROI);
start_motion_ind = max(diff(mean(start_tc,3),1)) < motionCutoff;
long_tc = reshape(stackGetTimeCourses(reshape(dFF_long,ypix,xpix,[]),mask),...
    nframes_2stim,sum(tr_long),nROI);
long_motion_ind = max(diff(mean(long_tc,3),1)) < motionCutoff;
tar_tc = reshape(stackGetTimeCourses(reshape(dFF_tar,ypix,xpix,[]),mask),...
    nframes_1s,sum(tr_tar_on),nROI);
tar_motion_ind = max(diff(mean(tar_tc,3),1)) < motionCutoff;

start_max = max(dFF_start_meanAllTrials(:,:,start_motion_ind),[],3);
long_max = max(dFF_long_meanAllTrials(:,:,long_motion_ind),[],3);
tar_max = max(dFF_tar_meanAllTrials(:,:,tar_motion_ind),[],3);

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

dFF_bxMax = cat(3,start_max,long_max,tar_max);

print(fullfile(fnout,[mouse '_' expDate '_maxImagesBx']),'-dpdf','-fillpage');
writetiff(max(dFF_bxMax,[],3),fullfile(fnout,[mouse '_' expDate '_maxImagesBx']))
%% save max images,reg outs
save(fullfile(fnout,'bx_max_images.mat'),'dFF_bxMax');
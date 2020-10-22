 %%variables and df/f
[nFrames_all nCells] = size(npSub_tc);
tGratingDirection = celleqel2mat_padded(input.tGratingDirectionDeg);
 nTrials = size(tGratingDirection,2);
 nOn = input.nScansOn;
 nOff = input.nScansOff;
 nFrames = nOn+nOff;
 sz = size(npSub_tc);
 
data_trial_beforereshape = reshape(npSub_tc,[nFrames nTrials nCells]);
 data_trial = permute(data_trial_beforereshape,[1 3 2]);
 
 data_f = mean(data_trial(ceil(nOff/2):nOff,:,:),1);
 data_dfof = (data_trial-data_f)./data_f;
 data_dfof_avg = squeeze(mean(data_dfof(nOff+1:(nOff+nOn),:,:),1));

 dirs = unique(tGratingDirection);
 nDir = length(dirs);

 
 %%divide into trial types
 iCell = 1:nCells;
[n n2] = subplotn(nDir);
 dfof_dir = zeros(nFrames,nCells,nDir);

figure;
for idir = 1:nDir
    subplot(2.*n, n2,idir)
    ind = find(tGratingDirection==dirs(idir));
    dfof_dir(:,:,idir) = squeeze(mean(data_dfof(:,:,ind),3));
    plot(dfof_dir(:,iCell,idir))
    title(num2str(dirs(idir)))
    ylabel('dF/F')
    xlabel('Frames')
end  
     
 

 
%% some variables from mworks
if iscell(input_tun.nScansOn)
    on = unique(cell2mat(input_tun.nScansOn))./down;
    off = unique(cell2mat(input_tun.nScansOff))./down;
else
    on = input_tun.nScansOn./down;
    off = input_tun.nScansOff./down; 
end
off_starts = 1:off+on:size(data_tun_reg,3);
on_starts = off+1:off+on:size(data_tun_reg,3);

tDir = cell2mat(input_tun.tGratingDirectionDeg);
[tDir_ind dir] = findgroups(tDir);
ntrials = length(tDir);
nstim = length(dir);

%% data by trial

if (off+on)*ntrials > nfr_tun 
    ntrials = floor(nfr_tun./(off+on));
    [tDir_ind dir] = findgroups(tDir(1:ntrials));
    data_tun_reg = data_tun_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_tun_reg,ypix,xpix,off+on,ntrials);
elseif (off+on)*ntrials < nfr_tun
    ntrials = length(tDir);
    [tDir_ind dir] = findgroups(tDir(1:ntrials));
    data_tun_reg = data_tun_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_tun_reg,ypix,xpix,off+on,ntrials);
else
    data_tr = reshape(data_tun_reg,ypix,xpix,off+on,ntrials);
end

%% dF/F by trial

F = nanmean(data_tr(:,:,floor(off/2):off,:),3);
dFF = bsxfun(@rdivide, bsxfun(@minus, data_tr, F), F);

motionCutoff = 0.2;
maxDFF = max(squeeze(mean(dFF(:,:,off+1:off+on,:),3)),[],3);
bwout = imCellEditInteractive(maxDFF);
mask = bwlabel(bwout);
nROI = length(unique(mask))-1;
tun_tc = reshape(stackGetTimeCourses(reshape(dFF,ypix,xpix,[]),mask),...
    off+on,ntrials,nROI);
motion_ind = max(diff(mean(tun_tc,3))) < motionCutoff;

dFF_dirmax = zeros(ypix,xpix,nstim);
for istim = 1:nstim
   dFF_dirmax(:,:,istim) = max(squeeze(mean(...
       dFF(:,:,off+1:off+on,tDir_ind == istim & motion_ind),3)),[],3); 
end
clear F dFF

%% figure

figure; setFigParams4Print('portrait')
colormap gray
imagesc(max(dFF_dirmax,[],3))
title([mouse '-' expDate '-tun'])
print(fullfile(fnout,[mouse '_' expDate '_maxImagesTun']),'-dpdf')
writetiff(max(dFF_dirmax,[],3), fullfile(fnout,[mouse '_' expDate '_maxImagesTun']))

%% save max images,reg outs
save(fullfile(fnout,'tun_max_images.mat'),'dFF_dirmax');

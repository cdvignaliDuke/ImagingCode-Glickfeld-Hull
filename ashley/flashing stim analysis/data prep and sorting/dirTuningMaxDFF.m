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

clear data_tun_reg
%% dF/F by trial

F = nanmean(data_tr(:,:,floor(off/2):off,:),3);
dFF = bsxfun(@rdivide, bsxfun(@minus, data_tr, F), F);

dFF_dirmax = zeros(ypix,xpix,nstim);
for istim = 1:nstim
   dFF_dirmax(:,:,istim) = max(squeeze(nanmean(dFF(:,:,off+1:off+on,tDir_ind == istim),3)),[],3); 
end
clear F dFF

%% figure

figure; setFigParams4Print('portrait')
colormap gray
imagesc(max(dFF_dirmax,[],3))
title([mouse '-' expDate '-tun'])
if ~exist(fullfile(fnout, 'max images'),'dir')
    mkdir(fullfile(fnout, 'max images'))
end
print(fullfile(fnout, 'max images',[mouse '_' expDate '_tun']),'-dpdf')
writetiff(max(dFF_dirmax,[],3), fullfile(fnout,'max images',[mouse '_' expDate '_tun']))

%% save max images,reg outs
if ~exist(fnout,'dir')
    mkdir(fnout)
end
save(fullfile(fnout,'tun_max_images.mat'),'dFF_dirmax');

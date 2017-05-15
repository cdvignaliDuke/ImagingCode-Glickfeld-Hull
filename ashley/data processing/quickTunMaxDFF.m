function maxDFF = quickTunMaxDFF(data_reg,input,down,tunVar,fnout,doSave)

%% image params
xpix = size(data_reg,2);
ypix = size(data_reg,1);
nfr = size(data_reg,3);

%% some variables from mworks
on = input.nScansOn./down;
off = input.nScansOff./down; 
off_starts = 1:off+on:nfr;
on_starts = off+1:off+on:nfr;

tStim = chop(cell2mat(eval(['input.' tunVar])),3);
[tStim_ind, stims] = findgroups(tStim);
ntrials = length(tStim);
nstim = length(stims);

%% data by trial

if (off+on)*ntrials > nfr 
    ntrials = floor(nfr./(off+on));
    [tStim_ind stims] = findgroups(tStim(1:ntrials));
    data_reg = data_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
elseif (off+on)*ntrials < nfr
    ntrials = length(tStim);
    [tStim_ind stims] = findgroups(tStim(1:ntrials));
    data_reg = data_reg(:,:,1:((off+on).*ntrials));
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
else
    data_tr = reshape(data_reg,ypix,xpix,off+on,ntrials);
end

clear data_reg
%% dF/F by trial

F = nanmean(data_tr(:,:,floor(off/2):off,:),3);
dFF = bsxfun(@rdivide, bsxfun(@minus, data_tr, F), F);

dFF_stimmax = zeros(ypix,xpix,nstim);
for istim = 1:nstim
   dFF_stimmax(:,:,istim) = max(squeeze(nanmean(dFF(:,:,off+1:off+on,tStim_ind == istim),3)),[],3); 
end
clear F dFF

%% figure
figure; setFigParams4Print('portrait')
colormap gray
imagesc(max(dFF_stimmax,[],3))

if doSave
if ~exist(fullfile(fnout, 'max images'),'dir')
    mkdir(fullfile(fnout, 'max images'))
end
print(fullfile(fnout, 'max images'),'-dpdf')
writetiff(max(dFF_stimmax,[],3), fullfile(fnout,'max images'))

%% save max images,reg outs
save(fullfile(fnout,'tun_max_images.mat'),'dFF_stimmax');

end
maxDFF = dFF_stimmax;
end

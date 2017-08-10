function maxDFF = quickExptMaxDFF(data_reg,input,down,exptVar,frRateHz,fnout,doSave)

%% image params
xpix = size(data_reg,2);
ypix = size(data_reg,1);
nfr = size(data_reg,3);

%% analysis params
bl_fr = frRateHz./down;
stim_fr = 2*bl_fr;
%% some trial variables
trStart_ind = double(cell2mat(eval(['input.' exptVar])));
trFrames = linspaceNDim(trStart_ind-bl_fr+1,trStart_ind+stim_fr,bl_fr+stim_fr);
ntrials = length(trStart_ind);
%% data by trial
if trFrames(end) > nfr
    trFrames = trFrames(:,ntrials-1);
end
data_tr = data_reg(:,:,trFrames(:));
clear data_reg
data_tr = reshape(data_tr,ypix,xpix,ntrials,bl_fr+stim_fr);

%% dF/F by trial

F = nanmean(data_tr(:,:,:,1:bl_fr),4);
dFF = bsxfun(@rdivide, bsxfun(@minus, data_tr, F), F);
dFFmean = mean(dFF,4);
clear F dFF
dffMax = max(dFFmean,[],3);

%% figure
figure; setFigParams4Print('portrait')
colormap gray
imagesc(dffMax)

if doSave
if ~exist(fullfile(fnout, 'max images'),'dir')
    mkdir(fullfile(fnout, 'max images'))
end
print(fullfile(fnout, 'max images'),'-dpdf')
writetiff(dffMax, fullfile(fnout,'max images'))

%% save max images,reg outs
save(fullfile(fnout,'tun_max_images.mat'),'dffMax');

end
maxDFF = dffMax;
end

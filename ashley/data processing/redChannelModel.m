ms = '1211';
dt = '190124';
t = '1109';
imgFolder = '003';
imgName = '003_000_000';
nfr = 5000;
doShortSession = true;

nBLFr = 30;
nPostStimFr = 30;
respwin = 35:37;
%%
fn = fullfile('Z:\home\ashley\data',ms,'two-photon imaging',dt,imgFolder);
fnout = fullfile('Z:\home\ashley\Analysis',ms,'two-photon imaging',dt,imgFolder);
cd(fn);
d = squeeze(sbxread(imgName,0,nfr));
tc = squeeze(mean(mean(d,1),2));

dg = squeeze(d(1,:,:,:));
dr = squeeze(d(2,:,:,:));

[ypix,xpix,~] = size(dg);

regimg = mean(dr(:,:,1500:1600),3);
figure;colormap gray;imagesc(regimg)
[outs,dr_reg] = stackRegister(dr,regimg);

%%
mw = loadMworksFile(ms,dt,t);

tTargetOn = cell2mat(mw.cTargetOn);
maxTrial = find(tTargetOn+nPostStimFr > nfr,1) - 1;

tTargetOn = tTargetOn(1:maxTrial);
tLeverDown = cell2mat(mw.cLeverDown(1:maxTrial));

firstAligned = nan(ypix,xpix,nBLFr+nPostStimFr,maxTrial);
targetAligned = nan(ypix,xpix,nBLFr+nPostStimFr,maxTrial);
for i = 1:maxTrial
    ind = (tLeverDown(i) - nBLFr):(tLeverDown(i)+nPostStimFr-1);
    firstAligned(:,:,:,i) = dg(:,:,ind);
    ind = (tTargetOn(i) - nBLFr):(tTargetOn(i)+nPostStimFr-1);
    targetAligned(:,:,:,i) = dg(:,:,ind);
end

F0 = mean(firstAligned(:,:,1:nBLFr,:),3);
dFF_first = (firstAligned - F0)./F0;
F0 = mean(targetAligned(:,:,1:nBLFr,:),3);
dFF_target = (targetAligned - F0)./F0;

maxDFF_first = max(squeeze(mean(dFF_first(:,:,respwin,:),4)),[],3);
figure;imagesc(maxDFF_first)
maxDFF_target = max(squeeze(mean(dFF_target(:,:,respwin,:),4)),[],3);
figure;imagesc(maxDFF_target)

mask_green = maskFromMultiMaxDFFStack(cat(3,maxDFF_first,maxDFF_target));
nCells = max(mask_green(:));
%% get average green and red values for each pixel in ROI
dg_vect = reshape(dg,[ypix*xpix,nfr]);
dr_vect = reshape(dr,[ypix*xpix,nfr]);
mask_vect = mask_green(:);

ROIs_red = cell(1,nCells);
ROIs_green = cell(1,nCells);
for i = 1:nCells
    ROIs_red{i} = mean(dg_vect(mask_vect == i,:),2);
    ROIs_green{i} = mean(dr_vect(mask_vect == i,:),2);
end
%%
mw_tun = loadMworksFile(ms,dt,'1211');
d_tun = squeeze(sbxread('004_000_000',0,28800));
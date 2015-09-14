mouse = 'AW33';
date = '150908';
expFolder = 'AW33_FS_Ret3pos_10hz_msLAT_1on3off_1';
fileName = 'AW33_FS_Ret3pos_10hz_msLAT_1on3off';


cd(['Z:\data\' mouse '\widefield imaging\' date '_' mouse '\' expFolder])
data = readtiff([fileName '.tif']);
% data = readtiff([expFolder '_MMStack.ome.tif']);
data = double(data);
%%
nON = 1;
nOFF = 3;
nStim = 3;
nRep = size(data,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);

% %% find dF/F accross images
% dataF = mean(data,3);
% dataDF = bsxfun(@minus,data,dataF);
% dataDFoverF = bsxfun(@rdivide,dataDF,dataF);

%%
%% initialize variables
offs = zeros(size(data,1),size(data,2), nOFF*nTrials);
ons = zeros(size(data,1),size(data,2),nTrials);

%% Pull out the off & on frames

for i = 1:nTrials
    offIND = 1 + ((nON+nOFF)*(i-1)); % indices for pulling off frames out
    offFrames = offIND:(offIND+(nOFF-1)); % off frame indicies from current trial
    offs(:,:,1+(nOFF*(i-1)):(1+(nOFF*(i-1)))+(nOFF-1)) = data(:,:,offFrames);    
end

offF = mean(offs,3);
dataDF = bsxfun(@minus,data,offF);
dataDFoverF = bsxfun(@rdivide,dataDF,offF);

for i = 1:nTrials
    onIND = nOFF+1 + ((nON+nOFF)*(i-1)); % indices for putting off frames in new structure
    onFrames = onIND:(onIND+nON-1); % on frame indicies from current trial
    ons(:,:,i) = mean(dataDFoverF(:,:,onFrames),3);
end


%%

onsF_sort = zeros(size(ons,1),size(ons,2),nRep,nStim);
for i = 1:nStim
    tempind = i:nStim:nTrials;
    onsF_sort(:,:,:,i) = ons(:,:,tempind);
end

%%
figure;
for i = 1:3
subplot(1,3,i)
imagesc(squeeze(mean(onsF_sort(:,:,:,i),3)))
colormap gray
img = squeeze(mean(onsF_sort(:,:,:,i),3));
imgname = ['img' num2str(i)];
writetiff(img,imgname);
end


%%
onINDstart = nOFF+1:nOFF+nON:size(data,3);

aroundOns = zeros(size(data,1),size(data,2),4,nTrials-1);
for i = 1:nTrials-1
    tempind = onINDstart(i)-1:onINDstart(i)+2;
    aroundOns(:,:,:,i) = dataDFoverF(:,:,tempind);
end

aroundOns_mean = mean(aroundOns,4);
aroundOnsTiffSize = reshape(aroundOns,[size(data,1),size(data,2),4*(nTrials-1)]);
writetiff(aroundOnsTiffSize,'framesaroundstimonsorted');
writetiff(aroundOns_mean,'framesaroundstimonmean');
%%
aroundOns_mean = squeeze(aroundOns_

%%
figure;
for i = 1:nTrials-1
    for iplot = 1:4
        subplot(4,29,iplot+((i-1)*4))
        imagesc(aroundOns(:,:,iplot,i));
        colormap gray
    end
end

    





















%% Script to Sort random trial according to coherence...






trialCohs = cell2mat(input.tDotCoherence); % lists coherences of trials in order
CoherenceSteps = unique(trialCohs); % coherence levels (ascending)
nCSteps = length(CoherenceSteps); % num of coherence levels

%% Figure out # of trials for each coherence level

for i=1:nCSteps
    nCohTrials(i) = length(find(trialCohs == CoherenceSteps(i)));
end

%% initialize variables

cohTrack = zeros(nCSteps,1); % # of times each level experienced
controls = zeros(size(data,1),size(data,2), size(data,3)/2);
cTrials = {}; % cell array holding matrix of frames for each coherence level
%zeros(size(data,1),size(data,2), size(data,3)/(2*nCSteps), nCSteps);

%% Pull out the control & trial frames

for i = 1:input.trialSinceReset
    offIND = 1 + (40*(i-1)); % indices for pulling control frames out
    onIND = 1 + (20*(i-1)); % indices for putting control frames in new structure
    offFrames = offIND:(offIND+19); % control frame indicies from current trial
    controls(:,:,onIND:(onIND+19)) = data(:,:,offFrames);
    cLevel = find(CoherenceSteps == trialCohs(i)); % check coherence level of trial
    
    l = 1 + (20*(cohTrack(cLevel)));
    % put 20 frames after control frames into cTrials
    cTrials{cLevel}(:,:,l:(l+19)) = data(:,:,(offFrames+20));
    
    cohTrack(cLevel) = cohTrack(cLevel) + 1;
end

%%
baselineF = mean(controls,3);
for i = 1:nCSteps
    dataDF{i} = bsxfun(@minus,double(cTrials{i}),baselineF);
    dataDFoverF{i} = bsxfun(@rdivide,double(dataDF{i}),baselineF);
end

figure;
for i = 1:nCSteps
    subplot(2,2,i)
    img = mean(dataDFoverF{i},3);
    imgname = ['img' num2str(i)];
    imagesc(img)
    colormap gray
    writetiff(img,imgname);
    title(['coherence = ' num2str(CoherenceSteps(i))])
end

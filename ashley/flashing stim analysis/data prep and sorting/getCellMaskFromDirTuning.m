iexp =4;
awFSAVdatasets_V1;
datasetStr = '_V1';
rc = behavConstsAV;
fnout = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'cells only');
%% load data get maxDFoverF image

    mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).dirtuning_time]; 
    load (fullfile(rc.behavData,mworks));
    nON = input.nScansOn/10;
    nOFF = input.nScansOff/10;
    nTrials = input.trialSinceReset;
    nStim = input.gratingDirectionStepN;
    data = readtiff(fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder,expt(iexp).date,expt(iexp).dirtuning,'DirectionTuning_V1.tif'));
    tStartInd = double(1:nON+nOFF:size(data,3));
    if ~isinteger(size(data,3)/nTrials)        
        nTrials = floor(size(data,3)/(nON+nOFF))
    end
    dataF = mean(reshape(data(:,:,linspaceNDim(tStartInd,tStartInd+(nOFF-1),nOFF)),size(data,1),size(data,2),nTrials,nOFF),4);
    dataDF = zeros(size(data));
    dFoverF = zeros(size(data));
    for itrial = 1:nTrials
        dataDF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)) = bsxfun(@minus,data(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)),dataF(:,:,itrial));
        dFoverF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)) = bsxfun(@rdivide,dataDF(:,:,tStartInd(itrial):tStartInd(itrial)+(nON+nOFF-1)),dataF(:,:,itrial));
    end


maxDFoverF = max(dFoverF,[],3);
% maxDF = max(dataDF,[],3);
figure; imagesc(maxDFoverF); colormap gray

writetiff(maxDFoverF,fullfile(fnout,'maxDFoverF'));
%% choose adjustment factor
adj = 1:3:18;
figure;
for i = 1:6
    subplot(3,2,i)
    
    image = maxDFoverF;
    x = adj(i)*median(image(:));
    image(image > x) = x;
    imagesc(image);colormap gray
    title(num2str(adj(i)));
end

%*****
ADJ = 4*median(image(:));
maxDFoverF(maxDFoverF > ADJ) = ADJ;

writetiff(maxDFoverF,fullfile(fnout,'maxDFoverF_ADJ'));
%% choose cells
bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data,mask_cell);
figure; tcOffsetPlot(data_TC)

save(fullfile(fnout,'mask&TCDir.mat'),'mask_cell','data_TC');
%% subtract neuropil
nCells = size(data_TC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);

npTC = zeros(size(data_TC));
for i = 1:size(data_TC,2)
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data,tempNPmask);
    end
end

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
npSubTC = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);

save(fullfile(fnout,'neuropil.mat'),'neuropil','npTC');

data_TC = npSubTC;
save(fullfile(fnout,'Timecourses.mat'),'data_TC');

down = 10;
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nStim = double(input.gratingAzimuthStepN)*double(input.gratingElevationStepN);

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,200:210),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% 

nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);


nOFF_ind = zeros(1,(nOFF*nStim*nRep));
nOFF_1 = zeros(1,nRep*nStim);
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    nOFF_1(1,iStim) = 1+((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end


F_data = mean(data_reg(:,:,nOFF_ind),3);

% dF average F
dF_data = bsxfun(@minus,data_reg, F_data);
dFoverF = bsxfun(@rdivide,dF_data,F_data);

maxDFoverF = max(dFoverF,[],3);
figure;imagesc(maxDFoverF);colormap(gray)

bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)

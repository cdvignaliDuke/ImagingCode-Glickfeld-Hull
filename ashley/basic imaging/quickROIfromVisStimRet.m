try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end

down = 10;
nON = double(input.nScansOn)./down;
nOFF = double(input.nScansOff)./down;
nTrials = length(input.tGratingDirectionDeg);

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,100:110),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

data_reg = data_reg(:,:,1:(nON+nOFF)*nTrials);
nOFF_ind = linspaceNDim(1:floor(nON+nOFF):size(data_reg,3),(1:floor(nON+nOFF):size(data_reg,3))+floor(nOFF-1),floor(nOFF));
nOFF_avg = mean(data_reg(:,:,nOFF_ind(:)),3);

% dF average F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF = bsxfun(@rdivide,dF_data,nOFF_avg);

maxDFoverF = max(dFoverF,[],3);
bwout = imCellEditInteractive(maxDFoverF);
mask_cell = bwlabel(bwout);
data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('mask&TCDir.mat','mask_cell','data_TC');
save('regImg.mat', 'data_avg');
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

data = data(:,:,1:4000);

% down = 10;
nTrials = sum(cell2mat(input.cLeverDown)<4000);

%average signals in time
data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data-down-min(min(min(data,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,301:400),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

leverDown = cell2mat(input.cLeverDown)
trialStart = double(leverDown(leverDown < 4000));
nOFF_ind = linspaceNDim(trialStart-29,trialStart,30);
nOFF_avg = double(mean(data_reg(:,:,nOFF_ind(:)),3));
data_reg = double(data_reg);

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
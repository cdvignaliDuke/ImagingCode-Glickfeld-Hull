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

writetiff(data,'rawF');
% data = readtiff('rawF.tif');


down = 10;

data_down = stackGroupProject(data,down);
% clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,100:110),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub
%%
F = mean(data_reg,3);
dFoverF = bsxfun(@rdivide,bsxfun(@minus,data_reg,F),F);

writetiff(data_reg,'dwnsmplF');
writetiff(dFoverF,'dFoverF');

max_dFoverF = max(dFoverF,[],3);

%%
bwout = imCellEditInteractive(max_dFoverF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data,mask_cell);
figure; tcOffsetPlot(data_TC)
save('rawCellTimecoursesAndMask.mat','mask_cell','data_TC')
writetiff(max_dFoverF,'max_dFoverF')


%% split into first and second half of experiment
L = size(data_reg,3)/2;

% first half dF/F
data_reg1 = data_reg(:,:,1:L);
F1 = mean(data_reg1,3);
dFoverF1 = bsxfun(@rdivide,bsxfun(@minus,data_reg1,F1),F1);

max_dFoverF1 = max(dFoverF1,[],3);

bwout1 = imCellEditInteractive(max_dFoverF1);
mask_cell1 = bwlabel(bwout1);

data_TC1 = stackGetTimeCourses(data_reg1,mask_cell1);
figure; tcOffsetPlot(data_TC1)

%second half
data_reg2 = data_reg(:,:,L+1:L*2);
F2 = mean(data_reg2,3);
dFoverF2 = bsxfun(@rdivide,bsxfun(@minus,data_reg2,F2),F2);

max_dFoverF2 = max(dFoverF2,[],3);

bwout2 = imCellEditInteractive(max_dFoverF2);
mask_cell2 = bwlabel(bwout2);

data_TC2 = stackGetTimeCourses(data_reg2,mask_cell);
figure; tcOffsetPlot(data_TC2)


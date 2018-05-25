clear all
close all
rc = behavConstsAV;

    subnum = '845';
    mouse = subnum;
    expDate = '180523';
    greenFolder = '001';
    nFrames = 1000;
    downSampleRate = 5;
    frameRateHz = 30;
    redFolder = '006';
    
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    %% load data & register
    fName = [greenFolder '_000_000'];
    retData = loadsbx_choosepmt(1,mouse,expDate,greenFolder,fName,nFrames);
    retData_down = stackGroupProject(retData,downSampleRate);
    clear retData
    
    regImg = mean(retData_down(:,:,1:10),3);
    figure;imagesc(regImg);colormap gray
    
    if ~exist(fullfile(fn,greenFolder),'dir')
        mkdir(fullfile(fn,greenFolder))
    end
    writetiff(regImg,fullfile(fn,'quickImg'))
    print(fullfile(fn,greenFolder,'quickImg'),'-dpdf','-fillpage')
    
    fName = [redFolder '_000_000'];
    redData = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName);
    redData_down = stackGroupProject(redData,downSampleRate);
    clear redData
    redImg = mean(redData_down,3);
    figure;imagesc(redImg);colormap gray
    
    writetiff(redImg,fullfile(fn,'quickImg_red'))
    print(fullfile(fn,greenFolder,'quickImg_red'),'-dpdf','-fillpage')
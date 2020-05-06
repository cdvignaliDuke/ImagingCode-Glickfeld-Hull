%% check behavior for running data


%% get path names
clear all;clc;

mouse = 'i838';
date = '180430';
ImgFolder = char('001');
time = char('2048');
doFromRef = 0;
ref = char('001');
nrun = size(ImgFolder,1);
frame_rate = 15;
run_str = catRunName(ImgFolder, nrun);

%% load input
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
load(fName);

%%
mouses = ['i838';'i840';'i842';'i843';'i844';'i880';'i881';'i883';'i884'];
track = [];
for iMouse = 1:length(mouses)
    mouse = mouses(iMouse,:)
    files = dir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-*-*.mat']);
    dates = string(ones(length(files),1));
    for iDate = 1:length(files)
        dates(iDate) = string(files(iDate).name(11:16));
    end
    dates = unique(dates)
    for iDate = 1:length(dates)
        date = char(dates(iDate));
        files = dir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-*.mat']);
        for iFile = 1:length(files)
            load(fullfile(files(iFile).folder,files(iFile).name));
            files(iFile).name
            gauss = input.doGaussianMask
            ellipse = input.doEllipseMask
            pause
        end
    end
end
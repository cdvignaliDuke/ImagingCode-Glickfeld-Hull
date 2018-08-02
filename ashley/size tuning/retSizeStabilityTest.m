%% get path names
clear all;clc;

ds = 'szTuning_dreadds_V1';
rc = behavConstsAV;
eval(ds)

for iexp = 1:size(expt,2)
mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
ret_runs = expt(iexp).retinotopyFolder;
ret_expTime = expt(iexp).retinotopyTime;
ret_nrun = length(ret_runs);
sz_runs = expt(iexp).sizeTuningFolder;
sz_expTime = expt(iexp).sizeTuningTime;
sz_nrun = length(sz_runs);
frame_rate = params.frameRate;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fnout)
    mkdir(fnout)
end

fprintf(['2P imaging retinotopy analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' expDate '\nExperiments:\n'])
for irun=1:ret_nrun
    fprintf([ret_runs{irun} ' - ' ret_expTime{irun} '\n'])
end

% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,ret_nrun);

fprintf(['\nBegin reading ' num2str(ret_nrun) ' runs...'])
for irun = 1:ret_nrun
    %load imaging data
    dataFolder = ret_runs{irun};
    fName = [dataFolder '_000_000'];
    switch expt(iexp).saveLoc
        case 'ashley'
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,500);
        case 'kevin'
            fdir = ['\\crash.dhe.duke.edu\data\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,500,fdir);
        otherwise
            error('identify data source')
    end
    data_avg = mean(data_temp,3);
    data = cat(3,data,data_avg);
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
clear data_temp

% size tuning

fprintf(['2P imaging size tuning analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' expDate '\nExperiments:\n'])
for irun=1:sz_nrun
    fprintf([sz_runs{irun} ' - ' sz_expTime{irun,:} '\n'])
end

% load data, read with sbxread, and concatenate selected runs

fprintf(['\nBegin reading ' num2str(sz_nrun) ' runs...'])
for irun = 1:sz_nrun
    %load imaging data
    dataFolder = sz_runs{irun};
    fprintf(['\nLoading run ' num2str(irun) '...'])
    fName = [dataFolder '_000_000'];
    switch expt(iexp).saveLoc
        case 'ashley'
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,500);
        case 'kevin'
            fdir = ['\\crash.dhe.duke.edu\data\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,500,fdir);
        otherwise
            error('identify data source')
    end
    data_avg = mean(data_temp,3);
    data = cat(3,data,data_avg);
    fprintf('Complete')
end
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])

writetiff(data,fullfile(fnout, [mouse '_' expDate '_stability.tif']));

end
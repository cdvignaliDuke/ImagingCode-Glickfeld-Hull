%% get path names
clear all;clc;

ds = 'szTuning_dreadds_V1';
rc = behavConstsAV;
eval(ds)

for iexp = 1:size(expt,1)
%%
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

%% load data, read with sbxread, and concatenate selected runs
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
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName);
        case 'kevin'
            fdir = ['\\crash.dhe.duke.edu\data\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,[],fdir);
        otherwise
            error('identify data source')
    end

    % load behavior data
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' mouse '-' expDate '-' time{irun} '.mat'];
    load(fName);

    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    ntrials = size(temp(irun).tGratingDirectionDeg,2);
    nframes = size(data_temp,3);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete\n')
end
fprintf('All runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
regIntv = 3000;
nep = floor(size(data,3)./regIntv);
fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

% plot 500 frame means at each register interval
[n, n2] = subplotn(nep);
figure(1);clf;
colormap(gray)
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

print(fullfile(fnout, [mouse '_' expDate '_ret_stability.pdf']))
%% size tuning

fprintf(['2P imaging size tuning analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' expDate '\nExperiments:\n'])
for irun=1:sz_nrun
    fprintf([sz_runs{irun} ' - ' sz_expTime{irun,:} '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,sz_nrun);

fprintf(['\nBegin reading ' num2str(sz_nrun) ' runs...'])
for irun = 1:nrun
    %load imaging data
    dataFolder = sz_runs{irun};
    fprintf(['\nLoading run ' num2str(irun) '...'])
    fName = [dataFolder '_000_000'];
    switch expt(iexp).saveLoc
        case 'ashley'
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName);
        case 'kevin'
            fdir = ['\\crash.dhe.duke.edu\data\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data_temp = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,[],fdir);
        otherwise
            error('identify data source')
    end

    % load behavior data
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' mouse '-' expDate '-' sz_expTime{irun} '.mat'];
    load(fName);
    ntrials = size(input.tGratingDirectionDeg,2);
    input.run = mat2cell(irun.*ones(1,ntrials),1,ntrials);
    temp(irun) = input;
    
    % store values on nOn + nOff, and measure number of trials
    nOn = temp(irun).nScansOn;
    nOff = temp(irun).nScansOff;
    nframes = size(data_temp,3);
    
    % squeeze because only 1 pmt channel
    data_temp = squeeze(data_temp);
    
    % if nframes =/= ntrials*(frames in trial), then resize
    if nframes>ntrials*(nOn+nOff)
        fprintf('Too many frames, truncating...\n')
        data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        nframes = (ntrials*(nOn+nOff))
    elseif nframes<ntrials*(nOn+nOff)
        fprintf('Too many trials, chop chop...\n')
        temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        ntrials = ceil(nframes./(nOn+nOff))
    end
    temp(irun).run = mat2cell(irun.*ones(1,ntrials),1,ntrials);
    data = cat(3,data,data_temp);
    trial_n(irun) = ntrials;
    fprintf('Complete')
end
fprintf('\nAll runs read\n')
fprintf([num2str(size(data,3)) ' total frames\n'])
input = concatenateDataBlocks(temp);
for i=1:length(input.tGratingContrast) % replace int64(con=1) with double
    if ~(class(input.tGratingContrast{i})=="double")
        input.tGratingContrast{i} = double(input.tGratingContrast{i});
    end
end
clear data_temp
clear temp

%% Check stability
regIntv = 3000;
nep = floor(size(data,3)./regIntv);
fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

% plot 500 frame means at each register interval
[n, n2] = subplotn(nep);
figure(1);clf;
colormap(gray)
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

print(fullfile(fnout, [mouse '_' expDate '_sz_stability.pdf']))

end
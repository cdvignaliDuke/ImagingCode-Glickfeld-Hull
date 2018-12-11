%% retOnly script by kevin 
% modified from singleChannelTCScript.m
% then modified from newScript.m
% then modified from retOnly.m 180801

% reads in Tuning Curves from single channel imaging data

%% get path names
close all;clear all;clc;

ds = 'szTuning_axons_AL';
iexp = 5; 
rc = behavConstsAV;
eval(ds)

%% conditionals for ashley analysis
doRegFrame = true;
doUsePreviousReg = false;
analyzer = 'ashley';

%%
mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
runs = expt(iexp).retinotopyFolder;
expTime = expt(iexp).retinotopyTime;
nrun = length(runs);
frame_rate = params.frameRate;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fnout)
    mkdir(fnout)
end

fprintf(['2P imaging retinotopy analysis - by KM, Glickfeld Lab\nSelected data:\nMouse: ' mouse '\nDate: ' expDate '\nExperiments:\n'])
for irun=1:nrun
    fprintf([runs{irun} ' - ' expTime{irun} '\n'])
end

%% load data, read with sbxread, and concatenate selected runs
data = [];
clear temp
trial_n = zeros(1,nrun);

fprintf(['\nBegin reading ' num2str(nrun) ' runs...'])
for irun = 1:nrun
    %load imaging data
    dataFolder = runs{irun};
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
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' mouse '-' expDate '-' expTime{irun} '.mat'];
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
if doRegFrame
    fprintf(['\nSplitting into ' num2str(nep) ' epochs of length ' num2str(regIntv) ' frames.\n'])

    % plot 500 frame means at each register interval
    [n, n2] = subplotn(nep);
    figure(1);clf;
    colormap(gray)
    for i = 1:nep
        subplot(n,n2,i);
        imagesc(mean(data(:,:,(1:500)+(i-1)*regIntv),3));
        title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
        clim([600 1500])
    end
end
%% Register data

chooseInt = 5; %nep/2 % interval chosen for data_avg =[epoch of choice]-1

fprintf('\nBegin registering...\n')
% if exist(fullfile(fnout,dataFolder), 'dir')
%     if exist(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']))
%     % checks if analysis already present
%     % load reg_shifts.mat (out, data_avg) and save the current input file
%     fprintf('Found previous analysis! Loading...\n')
%     
%     switch analyzer
%         case 'kevin'
%             load(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']))
%         case 'lindsey'
%             load(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']))
%         case 'ashley'
%             load(fullfile(fnout, dataFolder, 'registration'))
%     end
%     
%     % register
%     fprintf('stackRegister_MA, using shifts from previous registration\n')
%     % uses previous registration shifts (out) to re-register data quickly
%     switch analyzer
%         case 'ashley'
%             out = outs;
%     end
%     [outs, data_reg]=stackRegister_MA(data,[],[],double(out));
%     fprintf('Previous registration loaded...\n')
%     
%     % save new input
%     save(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']), 'input')
%     end
%     
% elseif ~exist(fullfile(fnout,dataFolder), 'dir') || ~exist(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']))
    % else means no previous analysis present
    % use data_avg selected above (could move here?)
    % then create new directory and save analysis
    fprintf('\nCreating new analysis!')
    
    meanrng = regIntv*(chooseInt)+(1:500);
    data_avg = mean(data(:,:,meanrng),3);
    fprintf(['\nRegister frame averaged from ' num2str(meanrng(1)) ' - ' num2str(meanrng(end)) '\n'])
    
    % register
    fprintf('stackRegister\n')
    [out, data_reg] = stackRegister(data,data_avg);
    
    % save
    fprintf('Registration complete, now saving...\n')
    mkdir(fullfile(fnout,dataFolder))
    save(fullfile(fnout, dataFolder, [mouse '_' expDate 'ret_reg_shifts.mat']), 'out', 'data_avg','meanrng')
    save(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']), 'input')
%end
clear data % depending on memory

%% test stability
% figure 2 shows the registered images to check the stability
fprintf('\nExamine registered images for stability\n')
figure(2);clf;
[n,n2] = optimizeSubplotDim(nep);
for i = 1:nep
    subplot(n,n2,i);
    imagesc(mean(data_reg(:,:,(1:500)+((i-1)*regIntv)),3));
    title([num2str(i) ': ' num2str(1+((i-1)*regIntv)) '-' num2str(500+((i-1)*regIntv))]);
end

% figure 3 shows imagesq (scaled?) of the registered data 1-10000, then prints to FOV_avg.mat
fprintf('Examine FOV\n')
figure(3);
if nframes>10000
    imagesq(mean(data_reg(:,:,1:10000),3));
else
    imagesq(mean(data_reg,3));
end
truesize;
% print to file
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOV_avg.pdf']),'-dpdf','-fillpage')

%% find activated cells
% calculate dF/F
fprintf('\nBegin image analysis...\n')

useFilt = 2;
reduceFlag = 4; % 1: 3x3 max project, 4: max project across El or Az
% max by trial type
if isfield(input, 'nScansOn')
    % nScansOn -> passive vis stim ret
    fprintf('nScansOn method - get dF/F\n')
    
    % load defining variables
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
    sz = size(data_reg);
    
    fprintf('Calculating dF/F...\n')
    
    % split into trials by reshape
    data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_fOn = mean(data_tr(:,:,nOff:(nOn+nOff),:),3);
    data_dfof = squeeze((data_fOn-data_f)./data_f);
    
    % previous code would do data_dfof for every frame
    % data_dfof = (double(data_tr)-data_f)./data_f;
    
    clear data_tr data_f data_fOn % data_fLate
    
    if useFilt
        fprintf('Filtering images and time averaging...\n')
        switch useFilt
            case 1
                % filter with a 20x20 gaussian sigma=0.7
                fprintf('useFilt=1: Gaussian smoothing filter\n')
                myfilter = fspecial('gaussian',[20 20], 0.7);
                data_dfof = imfilter(data_dfof,myfilter);
            case 2
                % filter with median filter (default 3x3)
                fprintf('useFilt=2: Median filter\n')
                data_dfof = medfilt3(data_dfof, [3 3 1]);
            case 3
                % filter with wiener adaptive filter (default 3x3)
                fprintf('useFilt=3: Wiener adaptive filter\n')
                %data_dfof = wiener2(data_dfof,[3 3]);
        end
    end
    
    fprintf('done\n')
end

% with dF/F, average by each stimulus, depending on experiment
if input.doRetStim
    % doRetStim -> retinotopy
    % requires data_dfof from above nScansOn method
    fprintf('input.doRetStim method - varying Az+El position\n')
    
    % store vectors for azimuth and elevation
    % what is celleqel2mat_padded?
    Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
    El = celleqel2mat_padded(input.tGratingElevationDeg);
    Azs = unique(Az);
    Els = unique(El);
    if min(Els,[],2)<0
        Els = fliplr(Els);
    end
    
    nStim = length(Azs)*length(Els);
    fprintf([num2str(nStim) ' unique Az+El stimuli\n'])
    
    fprintf('Averaging dF/F for each Az+El combination...\n')
    data_dfof_avg = zeros(sz(1), sz(2), nStim);
    start=1;
    Stims = [];
    %zeros(nStim, 2);
    
    for iEl = 1:length(Els)
        % runs through all unique Els to find indices of all same El
        indE = find(El == Els(iEl));
        for iAz = 1:length(Azs)
            % runs through all unique Azs to find indices of all same Az
            % then choose common El and Az indices
            % average over selected indices for all pixels and frames
            indA = find(Az == Azs(iAz));
            ind = intersect(indE,indA);
            data_dfof_avg(:,:,start) = mean(data_dfof(:,:,ind),3);
            
            % stores combination to Stims and iterates start
            Stims = [Stims; Els(iEl) Azs(iAz)];
            start = start+1;
        end
    end
    
    % plot dF/F for all stimuli
    figure(4);clf;
    for i = 1:nStim
        subplot(length(Els),length(Azs),i);
        imagesc(data_dfof_avg(:,:,i));
        clim([0 max(data_dfof_avg(:))])
        title(['El: ' num2str(Stims(i,1)) ', Az: ' num2str(Stims(i,2))])
    end
    % print to pdf
    print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOV_resp_Ret.pdf']), '-dpdf')
    
elseif input.doSizeStim
    % doSizeStim -> size tuning
    % requires data_dfof from above nScansOn method
    fprintf('input.doSizeStim method - varying grating diameter\n')
    fprintf('ERROR: this script only does retinotopy not size tuning')
    
elseif input.doTFStim && ~input.doMatrix
    % doTFStim + !doMatrix -> temporal frequency?
    % requires data_dfof from above nScansOn method
    fprintf('input.doTFStim method - varying grating TF+SF\n')
    fprintf('ERROR: this script only does retinotopy not TF tuning')
    
elseif input.doDirStim
    % doDirStim -> directional stimulation tuning curve
    % requires data_dfof from above nScansOn method
    fprintf('input.doDirStim method - varying grating direction\n')
    fprintf('ERROR: this script only does retinotopy not direction tuning')
end

% take max across stimuli
fprintf('Final step: take max across stimuli\n')
data_dfof_max = max(data_dfof_avg,[],3);
figure(10);clf;
imagesc(data_dfof_max)
clim([0 max(data_dfof_max(:))])
title('Maximum dF/F across all stimuli')
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOV_resp_max.pdf']), '-dpdf')

% save stimActFOV.mat containing: data_dfof_max, data_dfof_avg, nStim
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg', 'nStim')

if reduceFlag
    fprintf('Reducing data_dfof_avg, method: ')
    nStim_reduced = (length(Els)-2)*(length(Azs)-2); % default, unless case 4
    switch reduceFlag
        case 1
            fprintf('3x3 max project\n')
        case 2
            fprintf('3x3 pixel mean\n')
        case 3
            fprintf('3x3 pixel median\n')
        case 4
            fprintf('Max project across El,Az\n')
            nStim_reduced = length(Els)+length(Azs); % overwrite for different stim #s
    end
    data_dfof_avg_reduced = zeros(sz(1),sz(2),nStim_reduced);
    Stims_reduced = zeros(nStim_reduced,2);
    start=1;
    if reduceFlag <4 % 3x3 reduce
        for i = 1:length(Els)-2
            for j = 1:length(Azs)-2
                % find indices of 3x3 range to reduce
                reduce_El = ismember(Stims(:,1),Els(i+[0 1 2]));
                reduce_Az = ismember(Stims(:,2),Azs(j+[0 1 2]));
                reduce_ind = reduce_El.*reduce_Az;
                % reduce by median or mean or max
                switch reduceFlag
                    case 1 % max
                        data_dfof_avg_reduced(:,:,start) = max(data_dfof_avg(:,:,find(reduce_ind)),[],3);
                    case 2 % mean
                        data_dfof_avg_reduced(:,:,start) = mean(data_dfof_avg(:,:,find(reduce_ind)),3);
                    case 3 % median
                        data_dfof_avg_reduced(:,:,start) = median(data_dfof_avg(:,:,find(reduce_ind)),3);
                end
                Stims_reduced(start,:) = [Els(i+1), Azs(j+1)];
                start = start+1;
            end
        end
    else % reduceFlag >= 4, El/Az reduce (across all Az at each El, then vice versa)
        for i = 1:length(Els)
            % find indices of El to reduce by max
            reduce_ind = ismember(Stims(:,1),Els(i));
            data_dfof_avg_reduced(:,:,start) = max(data_dfof_avg(:,:,find(reduce_ind)),[],3);
            Stims_reduced(start,:) = [Els(i), NaN];
            start = start+1;
        end
        for j = 1:length(Azs)
            % find indices of Az to reduce by max
            reduce_ind = ismember(Stims(:,2),Azs(j));
            data_dfof_avg_reduced(:,:,start) = max(data_dfof_avg(:,:,find(reduce_ind)),[],3);
            Stims_reduced(start,:) = [NaN, Azs(j)];
            start = start+1;
        end
    end
end

% examine total image average over all trials
fprintf('Examine image average dF/F by trial\n')
tot_avg = squeeze(mean(mean(data_dfof,1),2));
figure;
plot(1:ntrials,tot_avg)
title('Total image average dF/F vs trial')
xlabel('Trial')
ylabel('dF/F')

fprintf('Examine image average dF/F by stimulus\n')
maxR = max(tot_avg(:));
tot_avg_sort = [tot_avg 0.5*El'/max(El)*maxR 0.5*Az'/max(Az)*maxR];
tot_avg_sort = sortrows(tot_avg_sort,[2 3]);
figure;
plot(1:ntrials,tot_avg_sort)
title('Total image average dF/F vs stimulus')
xlabel('Trial')
ylabel('dF/F')
legend('dF/F', 'El', 'Az', 'Location', 'se')

for i=1:length(Els)
    for j = 1:length(Azs)
        ind = (El==Els(i)).*(Az==Azs(j));
        respMap(i, j) = mean(tot_avg(find(ind)));
    end
end
figure;
ax=gca;
imagesc(respMap);
ax.XTickLabel = Azs;
ax.YTickLabel = Els;
title('Total image average dF/F response map')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_totalImageRF.pdf']), '-dpdf')


% scroll through stims
if 0
    %% scroll
    if reduceFlag
        data_temp = data_dfof_avg_reduced;
        Stims_temp = Stims_reduced;
        nStim_temp = nStim_reduced;
    else
        data_temp = data_dfof_avg;
        Stims_temp = Stims;
        nStim_temp = nStim;
    end
    fig=figure;
    ax1=axes('Parent',fig);
    i=1;
    while i <= nStim_temp
        %imagesc(ax,data_dfof_avg(:,:,i));
        %clim(ax,[0 max(data_dfof_avg(:))])
        imagesc(ax1,data_temp(:,:,i));
        clim(ax1,[0 max(data_temp(:))])
        title(ax1,['El: ' num2str(Stims_temp(i,1)) ', Az: ' num2str(Stims_temp(i,2))])
        
        was_a_key = waitforbuttonpress;
        if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
            i = i - 1;
            if ~i; i=1; end
        elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'x')
            break
        else
            i = i + 1;
        end
    end
    close(fig)
end

%% axon segmentation
% find hotspots

fprintf('\nBegin axon segmentation...')
min_df = 0.05;
b = 5;
mask_cell = zeros(sz(1), sz(2));
mask_all = zeros(sz(1), sz(2));
temp_data = ones(sz(1),sz(2));
temp_max = squeeze(max(max(data_dfof_avg,[],1),[],2));
[a, max_sort] = sort(temp_max,'descend');
data_dfof_avg_sort = cat(3,data_dfof_max,data_dfof_avg(:,:,max_sort));
for ia = 1:size(data_dfof_avg_sort,3)
    fprintf([ '\n img:' num2str(ia)])
    temp_data_log = ~isnan(temp_data);
    [x, y, v] = find(data_dfof_avg_sort(:,:,ia).*temp_data_log);
    [a, ind_sort] = sort(v,'descend');
    ind = find(a>min_df);
    fprintf(['- ' num2str(length(ind)) ' pix: '])
    for a = 1:length(ind)
        i = x(ind_sort(a));
        j = y(ind_sort(a));
        if i>5 & j>5 & i<sz(1)-b & j<sz(2)-b
            if ~isnan(temp_data(i-1:i+1,j-1:j+1))
                all_pix = data_dfof_avg_sort(i-1:i+1,j-1:j+1,ia);
                [max_val max_ind] = max(all_pix(:));
                if max_ind == 5
                    h = zeros(length(Els),length(Azs));
                    for iEl = 1:length(Els)
                        indE = find(El == Els(iEl));
                        for iAz = 1:length(Azs)
                            indA = find(Az == Azs(iAz));
                            ind = intersect(indE,indA);
                            [h(iEl,iAz) p] = ttest(squeeze(mean(mean(data_dfof(i-1:i+1,j-1:j+1,ind),1),2)),0,'tail','right','alpha',0.05./(length(Els)*length(Azs)));
                        end
                    end
                    if sum(h(:))>2
                        mask_cell(i,j) = 1;
                        mask_all(i-1:i+1,j-1:j+1) = ones(3,3);
                        temp_data(i-2:i+2,j-2:j+2) = NaN(5,5);
                        fprintf('.')
                    end
                end
            end
        end
    end
end      
   
figure; imagesc(mask_cell)
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_bouton_mask.pdf']), '-dpdf')

save(fullfile(fnout, dataFolder, [mouse '_' expDate '_mask_cell.mat']), 'mask_cell', 'mask_all')
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOVs.mat']), 'data_dfof_avg','data_dfof_avg_reduced', 'data_dfof_max','Azs', 'Els','Stims')

% look at masks overlaid on dfof average FOV for each stimulus
figure;
start = 1;
for i = 1:size(data_dfof_avg_sort,3)
    subplot(3,3,start);
    shade_img = imShade(data_dfof_avg_sort(:,:,i), mask_all);
    imagesc(shade_img)
    if i==1
        title('max')
    else
        title(num2str(Stims(i-1,:)))
    end
    clim([0 max(data_dfof_avg(:))]);
    colormap(gray)
    if start<9
        start = start+1;
    else
        figure;
        start = 1;
    end
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOVresp.pdf']), '-dpdf')


% clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir

%% Get time courses
nCells = sum(mask_cell(:));
fprintf(['Found ' num2str(nCells) ' boutons\n'])

fprintf('Extracting cell signal...\n')
data_tc = zeros(sz(3),nCells);
iC = 1;
for i = 1:sz(1)
    ind = find(mask_cell(i,:));
    if length(ind)>0
        for ii = 1:length(ind)
            fprintf([num2str(iC) ' '])
            j = ind(ii);
            data_tc(:,iC) = squeeze(mean(mean(data_reg(i-1:i+1,j-1:j+1,:),1),2));
            iC = 1+iC;
        end
    end
end
fprintf([num2str(nCells) ' total cells extracted\n'])
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs.mat']), 'data_tc')

fprintf('Calculating dF/F time courses...\n')
%get dF/F
tc_mat = zeros(nOn+nOff, nCells, ntrials);
fulltc_mat = zeros(nOn+nOff, 2, ntrials); % 1 for overall image, 2 for cells average
full_tc = squeeze(mean(mean(data_reg,1),2));
% reshape into trials, could do this in one line with reshape?
for itrial = 1:ntrials
    tc_mat(:,:,itrial) = data_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)),:);
    fulltc_mat(:,1,itrial) = full_tc(1+((itrial-1)*(nOn+nOff)):(itrial*(nOn+nOff)));
    fulltc_mat(:,2,itrial) = mean(tc_mat(:,:,itrial),2);
end
tc_f = mean(tc_mat(nOff/2:nOff,:,:),1);
fulltc_f = mean(fulltc_mat(nOff/2:nOff,:,:),1);
tc_dfof = (tc_mat - tc_f) ./ tc_f;
fulltc_dfof = (fulltc_mat - fulltc_f) ./ fulltc_f;
clear tc_mat tc_f fulltc_mat fulltc_f

fprintf('Time course extraction complete.\n')

%% retinotopy for these cells
% ADD IN: normalize each cell by its own max

% examine heat map of cell responses over all trials
% this allows examination of long-term effects over trials (ex: eye goop)
% responses are average of all stim on time (nOff+1:nOn+nOff) for a trial
fprintf('Examine cell responses across all trials\n')
trialRespMap = squeeze(mean(tc_dfof(nOff+1:nOn+nOff,:,:),1));
maxR = max(trialRespMap(:));
trialRespMap = [trialRespMap;El/max(El)*maxR;Az/max(Az)*maxR]; %add rows of normalized El and Az
figure;
imagesc(trialRespMap);
title('Avg dF/F resp by cell,trial')
xlabel('Trial')
ylabel('Cell')
colorbar

fprintf('Re-order by stimulus\n')
trialRespMapSort = sortrows(trialRespMap',[nCells+1 nCells+2])';
figure;
imagesc(trialRespMapSort);
title('Avg dF/F resp by cell,stimulus')
xlabel('Stim (El,Az)')
ylabel('Cell')
colorbar

%% calculate tuning mat and plot

%fprintf('\nPlotting timecourses and measuring stimOn response\n')
fprintf('\nMeasuring stimOn response\n')
tuning_mat = zeros(nStim, 2, nCells);
fulltuning_mat = zeros(nStim, 2, 2);
Ind_struct = [];
if nCells<36
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./frame_rate);
% figure;
% start = 1;
% f = 1;
for iCell = 1:nCells
%     if start >36 & start < 217
%         set(gcf, 'Position', [0 0 800 1000]);
%         print(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs' num2str(f) '.pdf']), '-dpdf')
%         start = 1;
%         f= f+1;
%         figure;
%     end
%     subplot(n, n2, start)
    for iStim = 1:nStim
        indA = find(Az == Stims(iStim,2));
        indE = find(El == Stims(iStim,1));
        ind = intersect(indE,indA);
%         if start <217
%             plot(tt', squeeze(mean(tc_dfof(:,iCell,ind),3)))
%             hold on
%         end
        tuning_mat(iStim,1,iCell) = mean(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),1),3);
        tuning_mat(iStim,2,iCell) = std(mean(tc_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
    end
%     ylim([-0.05 0.25])
%     vline(nOff)
    start = start + 1;
end
for iCell = 1:2
    for iStim = 1:nStim
        indA = find(Az == Stims(iStim,2));
        indE = find(El == Stims(iStim,1));
        ind = intersect(indE,indA);
        fulltuning_mat(iStim,1,iCell) = mean(mean(fulltc_dfof(nOff+1:nOn+nOff,iCell,ind),1),3);
        fulltuning_mat(iStim,2,iCell) = std(mean(fulltc_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
    end
end
% set(gcf, 'Position', [0 0 800 1000]);
% print(fullfile(fnout, dataFolder, [mouse '_' expDate '_TCs' num2str(f) '.pdf']), '-dpdf')

% fprintf('Plotting tuning maps\n')
% figure;
% start = 1;
% f = 1;
% for iCell = 1:nCells
%     if start >36
%         set(gcf, 'Position', [0 0 800 1000]);
%         print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning' num2str(f) '.pdf']), '-dpdf')
%         start = 1;
%         f= f+1;
%         figure;
%     end
%     subplot(n, n2, start)
%     ret_mat = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)]);
%     ret_mat = ret_mat';
%     imagesc(ret_mat)
%     colormap gray
%     %clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
%     %clim([0 chop(max(tuning_mat(:,1,iCell),[],1),2)])
%     title(num2str(chop(max(tuning_mat(:,1,iCell),[],1),2)))
%     start = start +1;
% end
% % set(gcf, 'Position', [0 0 800 1000]);
% % print(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning.mat']), 'tc_dfof', 'tuning_mat', 'Stims', 'Ind_struct')

% plot tc and ret_mat for full and cell avg
figure;
subplot(2,2,1)
ret_mat = reshape(fulltuning_mat(:,1,1), [length(Azs) length(Els)]);
ret_mat = ret_mat';
imagesc(ret_mat)
colormap gray
title(['Full image dF/F tuning map, max:' num2str(chop(max(fulltuning_mat(:,1,1),[],1),2))])
subplot(2,2,2)
ret_mat = reshape(fulltuning_mat(:,1,2), [length(Azs) length(Els)]);
ret_mat = ret_mat';
imagesc(ret_mat)
colormap gray
title(['All Cell Average dF/F tuning map, max:' num2str(chop(max(fulltuning_mat(:,1,2),[],1),2))])

% look at neuropil response
% remove top and bottom trials at each position
% look at standard deviation
% try 15-20 degree step in rough retinotopy

%% fit retinotopy data

close all

fprintf('\nBegin fitting retinotopy data...\n')

% fprintf('Plot tc_dfof for all stims of cell 10\n')
% figure;
% for iCell = 19
%     for iCond = 1:nStim
%         subplot(7,7,iCond)
%         ind_all = Ind_struct(iCond).all_trials;
%         plot(squeeze(tc_dfof(:,iCell,ind_all)))
%         ylim([-0.1 0.4])
%     end
% end

Fit_struct = [];
[AzAz, ElEl] = meshgrid(Azs,Els);
grid2.AzAz = AzAz;
grid2.ElEl = ElEl;

dAz = median(diff(Azs));
dEl = median(diff(Els));
Az_vec00 = Azs(1):(dAz/10):Azs(end);
El_vec00 = Els(1):(dEl/10):Els(end);
[AzAz00,ElEl00]=meshgrid(Az_vec00,El_vec00);
grid2.AzAz00 = AzAz00;
grid2.ElEl00 = ElEl00;
Nshuf = 500;
fprintf(['Nshuf = ' num2str(Nshuf) '\n'])
resp_dFoverF = squeeze(mean(tc_dfof(nOff:nOn+nOff,:,:),1));
base_dFoverF = squeeze(mean(tc_dfof(nOff/2:nOff,:,:),1));
p_ttest = zeros(nCells,nStim);
h_ttest = zeros(nCells,nStim);
h_all = zeros(1,nCells);

fprintf('Begin shuffling...\n')
%temp
nCells = size(resp_dFoverF,1);
% for iCell = 1:nCells; rs(iCell) = Fit_struct(iCell).True.s_.Rsq; end
% ind_cells = intersect(find(h_all),find(rs>0.7));
% nCells = length(ind_cells);
for count_shuf = 0:Nshuf
    fprintf(['count_shuf: ' num2str(count_shuf) '/' num2str(Nshuf) '\n'])
    Im_mat_USE = zeros(nCells, nStim);
    for iCond = 1:nStim
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;
            [h_ttest(:,iCond), p_ttest(:,iCond)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
    end
    
    %temp 
    for iCell = 1:nCells
        if count_shuf== 0
            if sum(h_ttest(iCell,:),2) == 0
                ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/2));
                if length(ind_p)<2
                    ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/3));
                    if length(ind_p)<3
                        ind_p = find(p_ttest(iCell,:)< 0.05./((nStim-1)/4));
                        if length(ind_p)<4
                            h_all(1,iCell) = 0;
                        else
                            h_all(1,iCell) = 1;
                        end
                    else
                        h_all(1,iCell) = 1;
                    end
                else
                    h_all(1,iCell) = 1;
                end
            else
                h_all(1,iCell) = 1;
            end
        end
        if count_shuf == 1 & iCell == 1
            x = zeros(1,nCells);
            for i = 1:nCells
                x(1,i) = Fit_struct(iCell).True.s_.Rsq;
            end
            n = length(find(x>=0.5));
            fprintf([num2str(n) '/' num2str(nCells) ' r-squared > 0.5\n'])
        end
        if count_shuf>0
            if h_all(1,iCell) == 0
                continue
            end
            if nCells > 2000 & Fit_struct(iCell).True.s_.Rsq < 0.7
                continue
            end
        end
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0
            b = reshape(a',length(Azs),length(Els));
            data = b';
            if count_shuf == 0
                PLOTIT_FIT = 0;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG_Ret_KM % modified due to error from file saving in script, saves to kevin analysis folder
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG_Ret_KM
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end
    end
    if PLOTIT_FIT & count_shuf == 0  
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_RFfits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
    end
end
fprintf('\nShuffling done, saving fit results\n')
s = whos('Fit_struct');
if s.bytes < 2300000000
    save(fullfile(fnout, dataFolder, [mouse '_' expDate '_Fit_struct.mat']), 'Fit_struct')
    fprintf('\nSaved all shuffles\n')
else 
    Fit_struct_sub = rmfield(Fit_struct,'Shuf');
    save(fullfile(fnout, dataFolder, [mouse '_' expDate '_Fit_struct_sub.mat']), 'Fit_struct_sub')
    fprintf('\nSaved only true fits\n')
end
%% assess fits
resp_ind = find(h_all); % h_all indicates responsive cell (by t-test against baseline)

fprintf('Assessing goodness of fit\n')
if Nshuf>1
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)
            eval('tmp = Fit_struct(iCell).True.s_.x;');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_50];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_50];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_10];');
            eval('tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_10];');
            % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
            fit_true_vec(iCell,:) = tmp;
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval('tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10];');
                eval('tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];');
                % A sigma_Az sigma_El Az0 El0 xi El50 Az50 El10 Az10
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end
    save(fullfile(fnout, dataFolder, [mouse '_' expDate '_fit_shuf.mat']), 'fit_shuf_vec', 'fit_true_vec')
    
    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    ind_shuf_lb = ceil(Nshuf*alpha_bound); % 0.025 percentile
    ind_shuf_ub = ceil(Nshuf*(1-alpha_bound)); % 0.975 percentile
    for iCell = 1:nCells
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp); % sort in order
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb); %lower 0.025
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub); %upper 0.975
            lbub_fits(iCell,count2,3) = mean(i); %mean
            lbub_fits(iCell,count2,5) = std(i); %stdev
        end
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:); % true (no shuffle)
    end
end

lbub_diff = lbub_fits(:,:,2)-lbub_fits(:,:,1);

goodfit_ind = [];
for iCell = 1:nCells
    if lbub_diff(iCell,4)<input.gratingAzimuthStepDeg*2
        if lbub_diff(iCell,5)<input.gratingAzimuthStepDeg*2
            goodfit_ind = [goodfit_ind iCell];
        end
    end
end
fprintf(['#Good cells = ' num2str(length(goodfit_ind)) ' (first pass)...\nNow removing RFs at ret perimeter\n'])

goodfit_ind2 = zeros(size(goodfit_ind));
for i=1:length(goodfit_ind)
    if sum(round(lbub_fits(goodfit_ind(i),4,4))==[min(Azs) max(Azs)])
        continue
    elseif sum(round(lbub_fits(goodfit_ind(i),5,4))==[min(Els) max(Els)])
        continue
    end
    goodfit_ind2(i) = goodfit_ind(i);
end
goodfit_ind2(goodfit_ind2==0) = [];
goodfit_ind = goodfit_ind2;
fprintf(['#Good cells = ' num2str(length(goodfit_ind)) ' (final)\nSaving good fits\n'])

fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_lbub_fits.mat']);
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind')

figure;
subplot(2,2,1)
for i = goodfit_ind
    plot(lbub_fits(i,4,4), lbub_fits(i,5,4), 'o')
    hold on;
end
xlim([min(Azs,[],2) max(Azs,[],2)])
ylim([min(Els,[],2) max(Els,[],2)])
axis equal
title('RF center (elim)')
subplot(2,2,2)
for i = goodfit_ind
    ellipse(lbub_fits(i,2,4), lbub_fits(i,3,4), 0, lbub_fits(i,4,4), lbub_fits(i,5,4));
    hold on;
end
axis equal
title('RF- 1 sigma (elim)')
subplot(2,2,3)
erat = (lbub_fits(goodfit_ind,3,4).^2./lbub_fits(goodfit_ind,2,4).^2);
erat(erat>1) = 1./erat(erat>1);
e = sqrt(1-erat);
hist(e)
title('eccentricity')
subplot(2,2,4)
cdiff = lbub_fits(goodfit_ind,3,4).^2-lbub_fits(goodfit_ind,2,4).^2;
cdiff(cdiff<0) = -cdiff(cdiff<0);
c = sqrt(cdiff);
hist(c)
title('linear eccentricity (elim)')
set(gcf, 'Position', [0 0 800 1000]);
fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_RFs.pdf']);
print(fn_out,'-dpdf')

figure;
subplot(2,2,1)
hist(lbub_fits(goodfit_ind,2,4))
title('Sigma azimuth (elim)')
subplot(2,2,2)
hist(lbub_fits(goodfit_ind,3,4))
title('Sigma elevation (elim)')
subplot(2,2,3)
a = lbub_fits(goodfit_ind,3,4).*lbub_fits(goodfit_ind,2,4).*pi;
hist(a)
title('Area (elim)')
subplot(2,2,4)
scatter(a, lbub_fits(goodfit_ind,1,4),'o')
xlabel('Area (elim)')
ylabel('Peak dF/F')
set(gcf, 'Position', [0 0 800 1000]);
fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_RFdists.pdf']);
print(fn_out,'-dpdf')

%% visualize retinotopic organization
% takes each of the goodfit_inds and colors masks by El+Az of RF center
mask_label = bwlabel(mask_all);
retMap_El = NaN(size(mask_cell));
retMap_Az = retMap_El;
for i=1:length(goodfit_ind)
    ind = find(mask_label == goodfit_ind(i));
    retMap_El(ind) = lbub_fits(goodfit_ind(i),5,4);
    retMap_Az(ind) = lbub_fits(goodfit_ind(i),4,4);
end

imAlpha=ones(size(retMap_El));
imAlpha(isnan(retMap_El))=0; % set all unmasked pixels to alpha=0

figure(1);clf;
colormap default
subplot(1,2,1)
imagesc(retMap_El,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by El')
h = colorbar;
ylabel(h,'El (deg)','Rotation',270.0,'VerticalAlignment','bottom')
%set(get(h,'label'),'string','El (deg)','Rotation',270.0); 
set(gca,'color',0*[1 1 1]);
subplot(1,2,2)
imagesc(retMap_Az,'AlphaData',imAlpha)
title('Retinotopy of goodfit cells by Az')
h = colorbar;
ylabel(h,'Az (deg)','Rotation',270.0,'VerticalAlignment','bottom')
set(gca,'color',0*[1 1 1]); 
set(gcf, 'Position', [100,300,1200,400])
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_retMap.pdf']), '-dpdf')

%% more plots
% right now just show raw data vs fit data at 5x5 resolution
% only cells with dF/F > 0.05

% load data
% _FOVs.mat: 'Azs', 'Els', 'Stims', 'data_dfof_avg'
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_FOVs.mat']))
% _Tuning.mat: 'tc_dfof', 'tuning_mat', 'Stims', 'Ind_struct'
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning.mat']))
% % _Fit_struct.mat: Fit_Struct?
% fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Fit_struct.mat']);
% load(fn_out)
% _lbub_fits.mat: 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind'
fn_out = fullfile(fnout, dataFolder, [mouse '_' expDate '_lbub_fits.mat']);
load(fn_out);

% use tuning_mat for dF/F response during stimOn (mean 1st column, stdev 2nd column)
% paired with Stims which is the El, Az

% examine raw data and fit data in 5x5 resolution
% include in title if good fit
% first cutoff at peak dF/F >0.05
% figure layout: 36 subplots per, using 2 per cell, so 18 cells per figure
figure;
sp = 1;
f = 1;
for iCell = 1:size(tuning_mat,3)
    if max(tuning_mat(:,1,iCell)) < 0.05
        continue
    else
        if sp >36
            set(gcf, 'Position', [0 0 800 1000]);
            print(fullfile(fnout, dataFolder, [mouse '_' expDate '_RawFit' num2str(f) '.pdf']), '-dpdf')
            sp = 1;
            f = f+1;
            figure;
        end
        subplot(6,6,sp)
        ret_raw = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)]);
        ret_raw = ret_raw';
        imagesc(ret_raw)
        colormap gray
        %clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
        %clim([0 chop(max(tuning_mat(:,1,iCell),[],1),2)])
        title(['#' num2str(iCell) ' ' num2str(chop(max(tuning_mat(:,1,iCell),[],1),2))])
        
        % fit data
        % plug in lbub_fits(i,:,4) (all 10 parameters, only 6 used, and 4 corresponds to truefit
        % and 25x2 Stims, with columns swapped ([Az El] instead of [El Az])
        % to get back 2D gaussian 
        subplot(6,6,sp+1)
        fit_mat = Gauss2D_ellipseMA(lbub_fits(iCell,:,4),[Stims(:,2) Stims(:,1)]);%get gaussian as 25x1 vector
        ret_fit = reshape(fit_mat, [length(Azs) length(Els)]);
        ret_fit = ret_fit';
        imagesc(ret_fit)
        colormap gray
        %clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
        %clim([0 chop(max(tuning_mat(:,1,iCell),[],1),2)])
        if sum(iCell==goodfit_ind)
            title(['**' num2str(chop(max(fit_mat,[],1),2))])
        else
            title(num2str(chop(max(fit_mat,[],1),2)))
        end
        
        sp = sp + 2;
    end
end
set(gcf, 'Position', [0 0 800 1000]);
print(fullfile(fnout, dataFolder, [mouse '_' expDate '_RawFit' num2str(f) '.pdf']), '-dpdf')


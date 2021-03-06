%% get path names
% date is day2
date = '200108'; 
ImgFolder = strvcat('003');
time = strvcat('1158');
mouse = 'i1316';
% set align to 1 if aligning any day to day 1. set to 0 if you are just
% working with day 1
alignToRef = 1;
% ref_date is day1
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
% Type in your own stuff here:
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jerry\2P_Imaging_Grace';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
% loading the original matfile with raw data
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
%     defining the frames that are presenting stimuli (On) and those that
%     are not (off)
     temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);
% this is the trial chopper to make sure you have 160 trials
        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end

%     defining the unregistered data
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data
% these are the 500 frames that you select to register the entire
% dataset to
data_avg = mean(data(:,:,6001:6500),3);

% this part proceeds if the data has been previously registered and needs
% to be reloaded
if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
%     this part proceeds if the data has not been previously registered.
%     Data_reg_avg will be important for the images later
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg(:,:,1:10000),3);
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% find activated cells
    
%max by trial type    
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
%         calculating dfof and formatting it as a four dim array
%         (sz1,sz1,frames,trails)
    if nOn>29
        sz = size(data_reg);
        data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
        data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_f data_df data_tr
    else
        sz = size(data_reg);
        data_tr = zeros(sz(1),sz(2), 100, ntrials-1);
        for itrial = 1:ntrials-1
            data_tr(:,:,:,itrial) = data_reg(:,:,((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)));
        end
        data_f = mean(data_tr(:,:,1:50,:),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_f data_df data_tr
    end
end

if input.doDirStim
%     defining each direction (16 in total)
    Dir = cell2mat(input.tGratingDirectionDeg);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
%     obtaining average dfof values for each direction
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
%     here is the data_dfof_max that we will use later
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    
    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i));
        clim([0 max(data_dfof_avg_all(:))])
        title(num2str(Dirs(i)))
        colormap(gray)
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3);
        end
    end
    
    figure;
    [n n2] = subplotn(nDirs/2);
    for i = 1:nStim/2
        subplot(n,n2,i)
        imagesc(data_dfof_avg_ori(:,:,i));
        clim([0 max(data_dfof_avg_ori(:))])
        title(num2str(Dirs(i)))
        axis off
    end
    subplot(n,n2,i+1)
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
    data_dfof = cat(3,data_dfof_avg_ori,max(data_dfof_avg_ori,[],3));
    
    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
end   
    
 save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')

%% cell segmentation 
if ~alignToRef
    mask_all = zeros(sz(1), sz(2));
    mask_data = data_dfof;

%     selecting active cells from an avg FOV for each direction-only done
%     on day 1 to create a cell mask that will overlap the FOVs of the
%     following days
    for iStim = 1:size(data_dfof,3)    
        mask_data_temp = mask_data(:,:,iStim);
        mask_data_temp(find(mask_all >= 1)) = 0;
        bwout = imCellEditInteractive(mask_data_temp);
        mask_2 = bwlabel(bwout);
        mask_all = mask_all+mask_2;
        close all
    end
    mask_cell = bwlabel(mask_all);
    mask_np = imCellNeuropil(mask_cell, 3, 5);
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

    clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
else
%     loading data_dfof_max and data_reg_avg from day 1 (_reg)
    load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_stimActFOV.mat']))
    load(fullfile(fnout, [day2 '_' mouse], [day2 '_' mouse '_' run_str], [day2 '_' mouse '_' run_str '_reg_shifts.mat']))
    data_dfof_max_reg = data_dfof_max;
    data_avg_reg = data_reg_avg;
%     loading data_dfof_max and data_reg_avg from day 1 (_ref)
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']))
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']))
    data_dfof_max_ref = data_dfof_max;
    data_avg_ref = data_reg_avg;
    
%     defining day 1 and day 2 data_reg_avg and taking out the brightest
%     points of the image
    reg = data_avg_reg;
    ref = data_avg_ref;
    reg(find(reg>2000)) = 0;
    reg = (reg./max(max(abs(reg))));
    ref(find(ref>2000)) = 0;
    ref = (ref./max(max(abs(ref)))); 
    
%     selecting the matching landmarks on day 1 and day 2 data_reg_avg and
%     shifting the position of day 2 data_reg_avg and data_dfof_max accordingly
    sz_target  = size(reg);
    [input_points, base_points] = cpselect(double(reg),double(ref),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    reg2ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    reg2ref_dfof = imtransform(double(data_dfof_max_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    
%     Red/green images created here. rgb_ref2ref is using data_reg_avg, and rgb_reg2ref_dfof is using
%     data_dfof_max
    rgb_reg2ref = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref_dfof = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref(:,:,1) = ref;
    rgb_reg2ref(:,:,2) = reg2ref;
    rgb_reg2ref_dfof(:,:,1) = data_dfof_max_ref;
    rgb_reg2ref_dfof(:,:,2) = reg2ref_dfof;
    figure; imagesc(rgb_reg2ref); title(['day 2 on day 1, data reg rgb overlay'])
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays']))
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_reg2ref_overlay.pdf']), '-dpdf','-bestfit')
    figure; imagesc(rgb_reg2ref_dfof); title(['day 2 on day 1, dfof rgb overlay'])
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_dfof_overlay.pdf']), '-dpdf','-bestfit')
    
%     the first figure I showed
    figure; 
    subplot(2,3,1); imagesc(ref); axis off; axis equal; title('day 1 data avg')
    subplot(2,3,2); imagesc(reg); axis off; axis equal; title('day 2 data avg')
    subplot(2,3,3); imagesc(reg2ref); axis off; axis equal; title('registered day 2 data avg')
    subplot(2,3,4); imagesc(rgb_reg2ref_dfof(:,:,1)); axis off; axis equal; title('day 1 data dfof')
    subplot(2,3,5); imagesc(data_dfof_max_reg); axis off; axis equal; title('day 2 data dfof')
    subplot(2,3,6); imagesc(rgb_reg2ref_dfof(:,:,2)); axis off; axis equal; title('registered day 2 data dfof')
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_rgb_individuals.pdf']), '-dpdf','-bestfit')

%    each frame within day 2 data_reg is shifted according to how the
%    landmarks match up
    for i = 1:nframes
        data_reg(:,:,i) = imtransform(double(data_reg(:,:,i)),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        if rem(i,50) == 0
            fprintf([num2str(i) '/n'])
        end
    end
    
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']), 'input_points', 'base_points', 'mytform', 'data_dfof_max_ref', 'ref', 'reg2ref', 'reg2ref_dfof');
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')
end

%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
%np_tc = stackGetTimeCourses(data_reg,mask_np);
clear np_tc np_tc_down
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')



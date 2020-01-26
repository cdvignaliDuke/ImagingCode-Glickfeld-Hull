%% get path names
date = '200120';
ImgFolder = strvcat('002');
time = strvcat('1445');
mouse = 'i1313';
alignToRef = 1;
ref_date = '200118';
ref_run = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
lg_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(lg_fn, mouse, date, ImgFolder(irun,:));
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    
    if isfield(input, 'cLeverUp') 
        if irun>1
            ntrials = size(input.trialOutcomeCell,2);
            for itrial = 1:ntrials
                temp(irun).cLeverDown{itrial} = temp(irun).cLeverDown{itrial}+offset;
                temp(irun).cFirstStim{itrial} = temp(irun).cFirstStim{itrial}+offset;
                temp(irun).cLeverUp{itrial} = temp(irun).cLeverUp{itrial}+offset;
                temp(irun).cStimOn{itrial} = temp(irun).cStimOn{itrial}+offset;
                temp(irun).cTargetOn{itrial} = temp(irun).cTargetOn{itrial}+offset;
            end
        end
    end
    offset = offset+nframes;
    
    
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

%% Choose register interval
nframes = 2000;
nep = floor(size(data,3)./nframes);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end

%% Register data

data_avg = mean(data(:,:,6001:6500),3);

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
% elseif doFromRef
%     ref_str = ['runs-' ref];
%     if size(ref,1)>1
%         ref_str = [ref_str '-' ref(size(ref,1),:)];
%     end
%     load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
%     [out, data_reg] = stackRegister(data,data_avg);
%     mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
%     save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
%     load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
%     load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
%     save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg,3);
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end
figure; imagesq(data_reg_avg); truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')
%% find activated cells
    
%max by trial type    
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);
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

if isfield(input,'tCyclesOn')
    tCyc = cell2mat(input.tCyclesOn);
    cStart = cell2mat(input.cFirstStim);
    cStim = cell2mat(input.cStimOn);
    cTarget = celleqel2mat_padded(input.cTargetOn);
    nTrials = length(tCyc);
    sz = size(data_reg);
    data_f = zeros(sz(1),sz(2),nTrials);
    data_base = zeros(sz(1),sz(2),nTrials);
    data_base2 = zeros(sz(1),sz(2),nTrials);
    data_targ = zeros(sz(1),sz(2),nTrials);
    for itrial = 1:nTrials
        data_f(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)-20:cStart(itrial)-1),3);
        data_base(:,:,itrial) = mean(data_reg(:,:,cStart(itrial)+9:cStart(itrial)+19),3);
        if cStim(itrial) > cStart(itrial)
            data_base2(:,:,itrial) = mean(data_reg(:,:,cStim(itrial)+9:cStim(itrial)+19),3);
        else
            data_base2(:,:,itrial) = nan(sz(1),sz(2));
        end
        if ~isnan(cTarget(itrial))
            data_targ(:,:,itrial) = mean(data_reg(:,:,cTarget(itrial)+9:cTarget(itrial)+19),3);
        else
            data_targ(:,:,itrial) = nan(sz(1),sz(2));
        end
    end
    data_base_dfof = (data_base-data_f)./data_f;
    data_base2_dfof = (data_base2-data_f)./data_f;
    data_targ_dfof = (data_targ-data_f)./data_f;
    baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
    dirs = unique(baseDir);
    ndir = length(dirs);
    targetDelta = round(celleqel2mat_padded(input.tGratingDirectionDeg),0)-baseDir;
    deltas = unique(targetDelta);
    nDelta = length(deltas);
    data_dfof_dir = zeros(sz(1),sz(2),ndir);
    data_dfof2_dir = zeros(sz(1),sz(2),ndir);
    [n n2] = subplotn(ndir);
    figure;
    for idir = 1:ndir
        ind = find(baseDir == dirs(idir));
        data_dfof_dir(:,:,idir) = mean(data_base_dfof(:,:,ind),3);
        data_dfof2_dir(:,:,idir) = nanmean(data_base2_dfof(:,:,ind),3);
        subplot(n,n2,idir)
        imagesc(data_dfof_dir(:,:,idir))
        title(dirs(idir))
    end
    data_dfof_dir_all = cat(3, data_dfof_dir, data_dfof2_dir);
    data_dfof_targ = zeros(sz(1),sz(2),nDelta);
    [n n2] = subplotn(nDelta);
    figure;
    for idir = 1:nDelta
        ind = find(targetDelta == deltas(idir));
        data_dfof_targ(:,:,idir) = nanmean(data_targ_dfof(:,:,ind),3);
        subplot(n,n2,idir)
        imagesc(data_dfof_targ(:,:,idir))
        title(deltas(idir))
    end
    data_dfof = cat(3,data_dfof_dir_all,data_dfof_targ);
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_max = max(imfilter(data_dfof,myfilter),[],3);
    figure;
    imagesc(data_dfof_max)
end

if isfield(input,'tLeftTrial')
    tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
    tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
    cStimOn = celleqel2mat_padded(input.cStimOn);
    cDecision = celleqel2mat_padded(input.cDecision);
    SIx = strcmp(input.trialOutcomeCell, 'success');
    nTrials = length(tLeftTrial);
    sz = size(data_reg);
    data_f = nan(sz(1),sz(2),nTrials);
    data_targ = nan(sz(1),sz(2),nTrials);
    data_resp = nan(sz(1),sz(2),nTrials);
    for itrial = 1:nTrials
        data_f(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)-20:cStimOn(itrial)-1),3);
        data_targ(:,:,itrial) = mean(data_reg(:,:,cStimOn(itrial)+5:cStimOn(itrial)+25),3);
        if cDecision+25<nframes
            data_resp(:,:,itrial) = mean(data_reg(:,:,cDecision(itrial)+5:cDecision(itrial)+25),3);
        end
    end
    data_targ_dfof = (data_targ-data_f)./data_f;
    data_resp_dfof = (data_resp-data_f)./data_f;
    indL = find(tLeftTrial);
    data_dfof_L = mean(data_targ_dfof(:,:,indL),3);
    indR = intersect(find(tGratingContrast==1),find(~tLeftTrial));
    data_dfof_R = mean(data_targ_dfof(:,:,indR),3);
    data_dfof_resp = mean(data_resp_dfof(:,:,find(SIx)),3);
    figure;
    subplot(2,2,1)
    imagesc(data_dfof_L)
    title('Left Stim')
    subplot(2,2,2)
    imagesc(data_dfof_R)
    title('Right Stim')
    subplot(2,2,3)
    imagesc(data_dfof_resp)
    title('Decision')
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_dFoF.pdf']), '-dpdf')
    data_dfof = cat(3, data_dfof_resp, cat(3,data_dfof_L,data_dfof_R));
    data_dfof_max = max(data_dfof,[],3);
    figure; imagesc(data_dfof_max)
end

if input.doRetStim
    Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
    El = celleqel2mat_padded(input.tGratingElevationDeg);
    Azs = unique(Az);
    Els = unique(El);
    if min(Els,[],2)<0
        Els = fliplr(Els);
    end
    nStim = length(Azs).*length(Els);
    Stims = [];
    data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
    start = 1;
    for iEl = 1:length(Els)
        ind1 = find(El == Els(iEl));
        for iAz = 1:length(Azs)
            Stims = [Stims; Els(iEl) Azs(iAz)];
            ind2 = find(Az == Azs(iAz));
            ind = intersect(ind1,ind2);
            data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
            start = start +1;
        end
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

    figure; 
    for i = 1:nStim; 
        subplot(length(Els),length(Azs),i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        colormap(gray)
        title(num2str(Stims(i,:)))
        clim([0 max(data_dfof_avg_all(:))])
    end
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')

elseif input.doSizeStim
    sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
    szs = unique(sz_mat);
    nStim = length(szs);
    Stims = szs;
    for isz = 1:nStim
        ind = find(sz_mat == szs(isz));
        data_dfof_avg(:,:,:,isz) = mean(data_dfof(:,:,:,ind),4);
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));
    clear data_dfof_avg

    figure; 
    [n n2] = subplotn(nStim);
    for i = 1:nStim; 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        title([num2str(szs(i)) ' deg'])
        clim([0 max(data_dfof_avg_all(:))])
        colormap(gray)
    end
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')
elseif input.doTFStim & ~input.doMatrix
    TF = cell2mat(input.tGratingTemporalFreqCPS);
    TFs = unique(TF);
    SF = cell2mat(input.tGratingSpatialFreqCPD);
    nStim = length(TFs);
    Stims = TFs;
    for iTF = 1:nStim
        ind = find(TF == TFs(iTF));
        data_dfof_avg(:,:,:,iTF) = mean(data_dfof(:,:,:,ind),4);
        SFs(iTF) = SF(ind(1));
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

    figure; 
    [n n2] = subplotn(nStim);
    for i = 1:nStim; 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i)); 
        title([num2str(TFs(i)./SFs(i)) ' deg/s'])
        clim([0 max(data_dfof_avg_all(:))])
        colormap(gray)
    end
elseif input.doDirStim
    Dir = cell2mat(input.tGratingDirectionDeg);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
    [n n2] = subplotn(nDirs);
    figure;
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir))
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
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
    title('Max')
    axis off
    data_dfof = cat(3,data_dfof_avg_ori,max(data_dfof_avg_ori,[],3));
end

 save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')

%% cell segmentation 
if ~alignToRef
    mask_all = zeros(sz(1), sz(2));
    %mask_data = squeeze(max(reshape(data_dfof_avg_all, [sz(1) sz(2) 2 nStim/2]),[],3));
    mask_data = data_dfof;
    % figure;
    % [n, n2] = subplotn(size(mask_data,3));
    % for iStim = 1:nStim
    %     subplot(n, n2, iStim)
    %     imagesc(mask_data(:,:,iStim))
    %     colormap gray
    % end

    for iStim = 1:size(data_dfof,3)    
        mask_data_temp = mask_data(:,:,iStim);
        mask_data_temp(find(mask_all >= 1)) = 0;
        bwout = imCellEditInteractiveLG(mask_data_temp);
        mask_2 = bwlabel(bwout);
        mask_all = mask_all+mask_2;
        close all
    end
    mask_cell = bwlabel(mask_all);


    figure; 
    [n n2] = subplotn(nStim);
    for i = 1:nStim; 
        subplot(n,n2,i); 
        shade_img = imShade(data_dfof_avg_all(:,:,i), mask_all);
        imagesc(shade_img)
        if input.doSizeStim
        title([num2str(szs(i)) ' deg'])
        elseif input.doRetStim
            title([num2str(Stims(i,:))])
        end
        clim([0 max(data_dfof_avg_all(:))])
        colormap(gray)
    end
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')


    mask_np = imCellNeuropil(mask_cell, 3, 5);
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

    clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
else
    data_dfof_max_reg = data_dfof_max;
    data_avg_reg = data_reg_avg;
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']))
    load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']))
    data_dfof_max_ref = data_dfof_max;
    data_avg_ref = data_reg_avg;
    
    reg = data_avg_reg;
    ref = data_avg_ref;
    reg(find(reg>1e5)) = 0;
    reg = (reg./max(max(abs(reg))));
    ref(find(ref>1e5)) = 0;
    ref = (ref./max(max(abs(ref)))); 
    
    figure; imagesc(reg); title('reg')
    figure; imagesc(ref); title('ref')
    
    sz_target  = size(reg);
    [input_points, base_points] = cpselect(double(reg),double(ref),'Wait', true);
    mytform = maketform('affine',input_points(1:3,:), base_points(1:3,:));
    reg2ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    reg2ref_dfof = imtransform(double(data_dfof_max_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    
    rgb_reg2ref = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref_dfof = zeros(sz_target(1), sz_target(2), 3);
    rgb_reg2ref(:,:,1) = ref;
    rgb_reg2ref(:,:,2) = reg2ref;
    rgb_reg2ref_dfof(:,:,1) = data_dfof_max_ref;
    rgb_reg2ref_dfof(:,:,2) = reg2ref_dfof;
    figure; imagesc(rgb_reg2ref)
    figure; imagesc(rgb_reg2ref_dfof)
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_registered_dfof_overlay.pdf']), '-dpdf','-bestfit')
    
    for i = 1:nframes
        data_reg(:,:,i) = imtransform(double(data_reg(:,:,i)),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        if rem(i,50) == 0
            fprintf([num2str(i) '\n'])
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

%% extract tuning
if isfield(input, 'nScansOn');
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    if nOn>29
        data_mat = zeros(nOn+nOff, nCells, ntrials);
        for itrial = 1:ntrials
            data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(nOff/2:nOff,:,:),1);
    else
        data_mat = zeros(100, nCells, ntrials-1);
        for itrial = 1:ntrials-1
            data_mat(:,:,itrial) = npSub_tc(((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(1:50,:,:),1);
    end
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f);
    
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = 50:60;
    resp_win = 70:90;
    base = squeeze(mean(data_dfof(base_win,:,:),1));
    resp = squeeze(mean(data_dfof(resp_win,:,:),1));
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right');
    no_match = find(isnan(x));
    max_dir = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
        subplot(n, n2, i)
            for idir = 1:ndir
                if nOn>29
                    ind = find(Dir == Dirs(idir));
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1));
                if h_dir(i,idir)
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_dir(i,:),2)>0
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_dir(i,:) = max_ind;
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning.pdf']),'-dpdf')

    
    nori = length(Dirs)/2;
    Ori = Dir;
    for iori = 1:nori
        ind = find(Dir == Dirs(iori+nori));
        Ori(ind) = Dirs(iori);
    end
    Oris = unique(Ori);
    h_ori = zeros(nCells, nori);
    p_ori = zeros(nCells, nori);
    ori_resp = zeros(nCells,nori);
    max_ori = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
            subplot(n, n2, i)
            for iori = 1:nori
                if nOn>29
                    ind = find(Ori == Oris(iori));
                else
                    ind = find(Ori(1:ntrials-1) == Oris(iori));
                end
                [h_ori(i,iori), p_ori(i,iori)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(nori-1));
                if h_ori(i,iori)
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_ori(i,:),2)>0
                temp_resp = ori_resp(i,:);
                temp_resp(find(h_ori(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_ori(i,:) = max_ind;
            else
                [max_val, max_ind] = max(ori_resp(i,:),[],2);
                max_ori(i,:) = max_ind;
            end
            title([num2str(Oris(max_ori(i,:))) ' deg'])
        end
    end
    
    good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]);
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','max_dir','h_dir', 'h_ori', 'max_ori','good_ind')
end

%% ori fitting
nOn = input.nScansOn;
nOff = input.nScansOff;
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
input.trialSinceReset = nTrials;

down = 10;
nframes = size(npSub_tc,1)./down;
nCells = size(npSub_tc,2);
data_tc_down = squeeze(mean(reshape(npSub_tc, [down,nframes,nCells]),1));

tuningDownSampFactor = down;
[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuningLG(data_tc_down,input,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')

%%
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
figure; 
if nCells<49
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
start = 1;
x = 0;
for ic = 1:nCells
    if start > 49
        suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    tit_str = num2str(chop(R_square(1,ic),2));
    if fitReliability(ic)<22.5
        tit_str = [tit_str '- R'];
    end
    title(tit_str)
    start = start+1;
end
suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);
edges = [0 22.5:45:180]; 
[bin ind_bin] = histc(prefOri(1,:),edges);
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for i = 1:length(bin)
    tunedCells{i} = intersect(find(ind_bin==i),ind_theta90);
end

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningInfo.mat']),...
    'prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells');

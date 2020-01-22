%% clear everything
clear all
clear global
%% get path names
date = '200106';
ImgFolder = strvcat('003');
run = strvcat('000');
time = strvcat('1251');
mouse = 'i1316';
doFromRef = 0;
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = ['Z:\All_staff\home\grace\2P_imaging\' mouse '\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_' run '.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_' run],0,nframes);
    
    
    
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:floor(nframes./(nOn+nOff)));
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
t = 3000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data
data_avg = mean(data(:,:,6001:6500),3);

% if exist(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
%     load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
%     save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
%     [outs, data_reg]=stackRegister_MA(data,[],[],out);
%     clear out outs  
% elseif doFromRef
%     ref_str = ['runs-' ref_run];
%     if size(ref_run,1)>1
%         ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
%     end
%     load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']))
%     [out, data_reg] = stackRegister(data,data_avg);
%     mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
%     save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
%     load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']))
%     load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_input.mat']))
%     save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
% else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    data_reg_avg = mean(data_reg,3);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg', 'data_reg_avg')
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
% end
clear data

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end
figure; imagesq(mean(data_reg(:,:,1:1000),3)); truesize;
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']), '-dpdf')
%% Register day x to day 1
%load reference data
if doFromRef
    ref_str = ['runs-' ref_run];
    if size(ref_run,1)>1
        ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
    end
    %reload and reregister reference day data
    refData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
    data_ref_avg = refData.data_reg_avg;
    
    sz = size(data_reg);
    Dir = cell2mat(input.tGratingDirectionDeg);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
    [n n2] = subplotn(nDirs);
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
    end
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);

    %create averages of current day and reference day
    target = data_ref_avg;
    target(find(target>1e5)) = 0;
    target = target./max(max(abs(target)));
    
    reg = data_reg_avg;
    reg(find(reg>1e5)) = 0;
    reg = reg./max(max(abs(reg)));
    
    [i_pts b_pts] = cpselect(double(reg), double(target),'Wait',true);
    sz_target = size(data_reg_avg);
    mytform = maketform('affine',i_pts(1:3,:),b_pts(1:3,:));
    img = imtransform(double(data_reg_avg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    dfof_reg2ref = imtransform(double(data_dfof_max),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    data_reg_to_ref = imtransform(double(data_reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfof_max.mat']), 'dfof_reg2ref')

    %register ref avg to reg avg 
%     [out_avg, img] = stackRegister(data_reg_avg, data_ref_avg);
%     [out, data_reg_to_ref] = stackRegister_MA(data_reg,[],[],repmat(out_avg,[size(data_reg,3) 1]));
	figure; 
    subplot(2,1,1)
    imagesc(data_ref_avg); truesize;
    title('Reference image')
    subplot(2,1,2)
    imagesc(data_reg_avg); truesize;
    title('Todays image')
    figure;
    imagesc(img); truesize;
    title('Registered image')
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_registered_FOV.pdf']), '-dpdf')
end
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
            data_tr(:,:,:,itrial) = data_reg(:,:,((itrial-1)*(nOn+nOff))+ceil(nOff/2):nOn+nOff+ceil(nOff/2)+((itrial-1)*(nOn+nOff)));
        end
        data_f = mean(data_tr(:,:,1:ceil(nOff/2),:),3);
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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_dFoF.pdf']), '-dpdf')
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
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')

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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')
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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_16Dirs.pdf']), '-dpdf')
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    data_dfof = cat(3,data_dfof_max, data_dfof_avg_all); 
    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
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
    for i = 1:nStim/2
        subplot(2,2,i)
        imagesc(data_dfof_avg_ori(:,:,i));
        clim([0 max(data_dfof_avg_ori(:))])
    end
    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')
end

save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%% axon segmentation
% max_interp = interp2(data_dfof_max);
% f1 = fspecial('average');   
% max_interp_sm = filter2(f1, max_interp);
% sz2 = size(max_interp);
% Xi = 1:2:sz2(1);
% Yi = 1:2:sz2(2);
% stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);
% stack_max_sm_long = reshape(stack_max_interp_sm,[sz(1)*sz(2) 1]);
% 
% local_max = zeros(sz(1), sz(2));
% mask = zeros(sz(1), sz(2));
% yborder = 10;
% xborder = 40;
% rg = std(reshape(stack_max_interp_sm(yborder:(sz(1)-yborder),xborder:(sz(2)-xborder)), [1 length(yborder:(sz(1)-yborder)).*length(xborder:(sz(2)-xborder))]),[],2);
% hb = zeros(size(local_max));
% ht = zeros(size(local_max));
% pb = zeros(size(local_max));
% pt = zeros(size(local_max));
% for iy = yborder:(sz(1)-yborder);
%     for ix = xborder:(sz(2)-xborder);
%         if stack_max_interp_sm(iy,ix)> (3.*rg)
%             sub = stack_max_interp_sm(iy-2:iy+2,ix-2:ix+2);
%             sub_long = reshape(sub, [1 25]);
%             [sub_long_order ind_order] = sort(sub_long);
%             if ind_order(end)==13
%                 [hb(iy,ix) pb(iy,ix)] = ttest(mean(mean(data_base(iy-1:iy+1,ix-1:ix+1,:),1),2), mean(mean(data_f(iy-1:iy+1,ix-1:ix+1,:),1),2), 'tail', 'right');
%                 [ht(iy,ix) pt(iy,ix)] = ttest(mean(mean(data_targ(iy-1:iy+1,ix-1:ix+1,:),1),2), mean(mean(data_f(iy-1:iy+1,ix-1:ix+1,:),1),2), 'tail', 'right');
%                 if hb(iy,ix) + ht(iy,ix) > 0
%                     local_max(iy,ix) = 1;
%                     mask(iy-1:iy+1,ix-1:ix+1) = 1;
%                 end
%             end
%         end
%     end
% end
% n_pix = sum(sum(local_max));
% [i, j] = find(local_max ==1); 
% cell_mask = bwlabel(mask);
% data_tc = stackGetTimeCourses(data_reg, cell_mask);
% nCells = n_pix;
% npSub_tc = data_tc;
% 
% save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'local_max', 'cell_mask', 'n_pix', 'i', 'j')
% save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'npSub_tc')
% save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
% clear data_base data_base2 data_dfof data_targ data_base_dfof data_base2_dfof data_targ_dfof data_f data_dfof_dir data_reg

% %% cell segmentation 
% 
% sz = size(data_reg);
% mask_all = zeros(sz(1), sz(2));
% %mask_data = squeeze(max(reshape(data_dfof_avg_all, [sz(1) sz(2) 2 nStim/2]),[],3));
% mask_data = data_dfof;
% % figure;
% % [n, n2] = subplotn(size(mask_data,3));
% % for iStim = 1:nStim
% %     subplot(n, n2, iStim)
% %     imagesc(mask_data(:,:,iStim))
% %     colormap gray
% % end
% 
% 
% for iStim = 1:size(data_dfof,3)    
%     mask_data_temp = mask_data(:,:,iStim);
%     mask_data_temp(find(mask_all >= 1)) = 0;
%     bwout = imCellEditInteractive(mask_data_temp);
%     mask_2 = bwlabel(bwout);
%     mask_all = mask_all+mask_2;
%     close all
% end
% mask_cell = bwlabel(mask_all);
% 
% mask_data_temp = data_dfof_max;
% mask_data_temp(find(mask_cell >= 1)) = 0;
% figure; imagesc(mask_data_temp)
% print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_data.pdf']), '-dpdf')
% 
% 
% figure; 
% [n n2] = subplotn(nStim);
% for i = 1:nStim; 
%     subplot(n,n2,i); 
%     shade_img = imShade(data_dfof_avg_all(:,:,i), mask_all);
%     imagesc(shade_img)
%     if input.doSizeStim
%     title([num2str(szs(i)) ' deg'])
%     elseif input.doRetStim
%         title([num2str(Stims(i,:))])
%     end
%     clim([0 max(data_dfof_avg_all(:))])
%     colormap(gray)
% end
%     print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')
% 
% mask_np = imCellNeuropil(mask_cell, 3, 5);
% save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
% 
% clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

%% neuropil subtraction
if doFromRef
    ref_str = ['runs-' ref_run];
    if size(ref_run,1)>1
        ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
    end
refMaskData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
mask_cell =  refMaskData.mask_cell;
data_reg = data_reg_to_ref;
end

mask_np = imCellNeuropil(mask_cell, 3, 5);
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
[max_skew, ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down data_reg_to_ref

save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')


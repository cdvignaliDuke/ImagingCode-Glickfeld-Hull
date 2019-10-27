%% get path names
date = '190906';
ImgFolder = strvcat('002');
time = strvcat('1353');
mouse = 'i1306';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = ['Z:\All_staff\home\grace\2P_imaging\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    
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

if exist(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
elseif doFromRef
    ref_str = ['runs-' ref];
    if size(ref,1)>1
        ref_str = [ref_str '-' ref(size(ref,1),:)];
    end
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_mask_cell.mat']))
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ref_str], [date '_' mouse '_' ref_str '_trialData.mat']))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.mat']))
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
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    data_dfof = cat(3,data_dfof_max, data_dfof_avg_all); 
%     figure; 
%     Stims = Dirs;
%     nStim = length(Dirs);
%     [n n2] = subplotn(nDirs);
%     data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
%     for i = 1:nStim 
%         subplot(n,n2,i); 
%         imagesc(data_dfof_avg_all(:,:,i));
%         clim([0 max(data_dfof_avg_all(:))])
%         title(num2str(Dirs(i)))
%         colormap(gray)
%         if i<(nDirs/2)+1
%             data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3);
%         end
%     end
%     for i = 1:nStim/2
%         subplot(2,2,i)
%         imagesc(data_dfof_avg_ori(:,:,i));
%         clim([0 max(data_dfof_avg_ori(:))])
%     end
%     figure;
%     imagesc(max(data_dfof_avg_ori,[],3))
end

 save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
%% axon segmentation
max_interp = interp2(data_dfof_max);
f1 = fspecial('average');   
max_interp_sm = filter2(f1, max_interp);
sz2 = size(max_interp);
Xi = 1:2:sz2(1);
Yi = 1:2:sz2(2);
stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);
stack_max_sm_long = reshape(stack_max_interp_sm,[sz(1)*sz(2) 1]);

local_max = zeros(sz(1), sz(2));
mask = zeros(sz(1), sz(2));
yborder = 10;
xborder = 40;
rg = std(reshape(stack_max_interp_sm(yborder:(sz(1)-yborder),xborder:(sz(2)-xborder)), [1 length(yborder:(sz(1)-yborder)).*length(xborder:(sz(2)-xborder))]),[],2);
hb = zeros(size(local_max));
ht = zeros(size(local_max));
pb = zeros(size(local_max));
pt = zeros(size(local_max));
for iy = yborder:(sz(1)-yborder);
    for ix = xborder:(sz(2)-xborder);
        if stack_max_interp_sm(iy,ix)> (3.*rg)
            sub = stack_max_interp_sm(iy-2:iy+2,ix-2:ix+2);
            sub_long = reshape(sub, [1 25]);
            [sub_long_order ind_order] = sort(sub_long);
            if ind_order(end)==13
                [hb(iy,ix) pb(iy,ix)] = ttest(mean(mean(data_base(iy-1:iy+1,ix-1:ix+1,:),1),2), mean(mean(data_f(iy-1:iy+1,ix-1:ix+1,:),1),2), 'tail', 'right');
                [ht(iy,ix) pt(iy,ix)] = ttest(mean(mean(data_targ(iy-1:iy+1,ix-1:ix+1,:),1),2), mean(mean(data_f(iy-1:iy+1,ix-1:ix+1,:),1),2), 'tail', 'right');
                if hb(iy,ix) + ht(iy,ix) > 0
                    local_max(iy,ix) = 1;
                    mask(iy-1:iy+1,ix-1:ix+1) = 1;
                end
            end
        end
    end
end
n_pix = sum(sum(local_max));
[i, j] = find(local_max ==1); 
cell_mask = bwlabel(mask);
data_tc = stackGetTimeCourses(data_reg, cell_mask);
nCells = n_pix;
npSub_tc = data_tc;

save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'local_max', 'cell_mask', 'n_pix', 'i', 'j')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'npSub_tc')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
clear data_base data_base2 data_dfof data_targ data_base_dfof data_base2_dfof data_targ_dfof data_f data_dfof_dir data_reg
%% cell segmentation 

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
    bwout = imCellEditInteractive(mask_data_temp);
    mask_2 = bwlabel(bwout);
    mask_all = mask_all+mask_2;
    close all
end
mask_cell = bwlabel(mask_all);

mask_data_temp = data_dfof_max;
mask_data_temp(find(mask_cell >= 1)) = 0;
figure; imagesc(mask_data_temp)


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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')


mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 

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

save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')

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
    h = zeros(nCells, ndir);
    p = zeros(nCells, ndir);
    base = squeeze(mean(data_dfof(40:50,:,:),1));
    resp = squeeze(mean(data_dfof(60:70,:,:),1));
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right');
    max_dir = zeros(nCells,1);
    figure;
    start = 1;
    for i = 1:nCells
        if start >16
            start = 1;
            suptitle(['Cells' num2str(i-16) '-' num2str(i-1)])
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning_Cells' num2str(i-17) '-' num2str(i-1) '.pdf']),'-dpdf')
            figure;    
        end
            subplot(4, 4, start)
        for idir = 1:ndir
            if nOn>29
                ind = find(Dir == Dirs(idir));
            else
                ind = find(Dir(1:ntrials-1) == Dirs(idir));
            end
            [h(i,idir), p(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1));
            if h(i,idir)
                errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
            else
                errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
            end
            dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2);
            hold on
        end
        if sum(h(i,:),2)>0
            temp_resp = dir_resp(i,:);
            temp_resp(find(h(i,:)==0)) = NaN;
            [max_val max_ind] = max(temp_resp,[],2);
            max_dir(i,:) = max_ind;
        else
            [max_val max_ind] = max(dir_resp(i,:),[],2);
            max_dir(i,:) = max_ind;
        end
        title([num2str(Dirs(max_dir(i,:))) ' deg'])
        start = start+1;
    end
    good_ind = intersect(find(x), find(sum(h,2)>0));
    suptitle(['Cells' num2str(i-16) '-' num2str(i-1)])
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning_Cells' num2str(i-16) '-' num2str(i-1) '.pdf']),'-dpdf')
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','max_dir','good_ind')
end

%% Cycle analysis
if isfield(input, 'nFramesOn');
    tCyc = cell2mat(input.tCyclesOn);
    if length(unique(tCyc)) == 1
        nOn = input.nFramesOn;
        nOff = input.nFramesOff;
        cStart = cell2mat(input.cLeverDown);
        nTrials = length(tCyc);
        sz = size(npSub_tc);
        data_trial = zeros((nOn+nOff)*(tCyc(1)+1)+100,sz(2),nTrials);
        for itrial = 1:nTrials
            data_trial(:,:,itrial) = npSub_tc(cStart(itrial)-20:79+cStart(itrial)+((nOn+nOff)*(tCyc(itrial)+1)),:);
        end
        data_trial_f = mean(data_trial(1:20,:,:),1);
        data_trial_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_trial_f),data_trial_f);
        baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
        dirs = unique(baseDir);
        ndir = length(dirs);
        data_dfof_dir = zeros(size(data_trial,1),sz(2),ndir);
        nCells = sz(data_trial,2);
        h = zeros(ndir,nCells);
        p = zeros(ndir,nCells);
        
        for idir = 1:ndir
            ind = find(baseDir == dirs(idir));
            data_trial_dir(:,:,idir) = mean(data_trial_dfof(:,:,ind),3);
            for iCell = 1:nCells
                [h(idir,iCell), p(idir,iCell)] = ttest(mean(data_trial_dfof(1:20,iCell,ind),1),mean(data_trial_dfof(30:40,iCell,ind),1),'tail','left','alpha',0.05./(ndir-1));
            end
        end
        good_ind = find(sum(h,1)>0);
            
        figure;
        for i = 1:length(good_ind)
            [n n2] = subplotn(length(good_ind));
            subplot(n,n2,i); 
            iC = good_ind(i);
            plot(squeeze(data_trial_dir(:,iC,:))); 
            vline([21:nOn+nOff:21+(nOn+nOff)*(tCyc(1)-1)]); 
            vline(21+(nOn+nOff)*(tCyc(1)), 'k'); 
            title(num2str(iC))
        end
        
        figure; 
        [n n2] = subplotn(length(good_ind));
        for i = 1:length(good_ind)
            subplot(n,n2,i); 
            iC = good_ind(i);
            for idir = 1:ndir
                ind = find(baseDir == dirs(idir));
                if h(idir,iC)
                    errorbar(dirs(idir),squeeze(mean(data_trial_dir(31:35,iC,idir),1)),squeeze(std(mean(data_trial_dfof(31:35,iC,ind),1),[],3)./sqrt(length(ind))), 'or');
                else
                    errorbar(dirs(idir),squeeze(mean(data_trial_dir(31:35,iC,idir),1)),squeeze(std(mean(data_trial_dfof(31:35,iC,ind),1),[],3)./sqrt(length(ind))), 'ok');
                end
                hold on
            end
            ylim([0 max(mean(data_trial_dir(31:35,iC,:),1),[],3)*1.5])
            xlim([-dirs(2) dirs(end)+dirs(2)])
            title(num2str([osi(iC,1) chop(osi(iC,2),2)]))
        end
        
        
        figure;
        [n n2] = subplotn(length(good_ind));
        data_dir_avg = zeros(nCells, ndir, 6);
        max_dir = zeros(nCells,1);
        min_dir = zeros(nCells,1);
        for i = 1:length(good_ind)
            subplot(n,n2,i); 
            iC = good_ind(i);
            for idir = 1:ndir
                ind = find(baseDir == dirs(idir));
                data_dir_avg(iC,idir,1) = mean(data_trial_dir(31:35,iC,idir),1);
                data_dir_avg(iC,idir,2) = mean(data_trial_dir(42:46,iC,idir),1)-mean(data_trial_dir(37:39,iC,idir),1);
                data_dir_avg(iC,idir,3) = mean(data_trial_dir(53:57,iC,idir),1)-mean(data_trial_dir(48:50,iC,idir),1);
                data_dir_avg(iC,idir,4) = mean(data_trial_dir(64:68,iC,idir),1)-mean(data_trial_dir(59:61,iC,idir),1);
                data_dir_avg(iC,idir,5) = mean(data_trial_dir(75:79,iC,idir),1)-mean(data_trial_dir(70:72,iC,idir),1);
                data_dir_avg(iC,idir,6) = mean(data_trial_dir(86:90,iC,idir),1)-mean(data_trial_dir(81:83,iC,idir),1);
            end
            h_ind = find(h(:,iC));
            [max_val max_ind] = max(data_dir_avg(iC,h_ind,1),[],2);
            max_dir(iC,:) = h_ind(max_ind);
            plot(squeeze(data_trial_dir(11:35,iC,max_dir(iC,:))));
            hold on
            plot(squeeze(data_trial_dir(21+(4*(nOn+nOff))-10:21+14+(4*(nOn+nOff)),iC,max_dir(iC,:))-mean(data_trial_dir(21+(4*(nOn+nOff))-10:21+(4*(nOn+nOff)),iC,max_dir(iC,:)),1)));
            ind_maxdir = find(baseDir == dirs(max_dir(iC,:)));
            tarDir = celleqel2mat_padded(input.tGratingDirectionDeg);
            tars = unique(tarDir(ind_maxdir));
            ntar = length(tars);
%             for itar = 1:ntar
%                 ind_tar = intersect(find(tarDir == tars(itar)),ind_base);
%                 plot(squeeze(mean(data_trial_dfof(21+(5*(nOn+nOff))-10:21+14+(5*(nOn+nOff)),iC,ind_tar),3)- mean(mean(data_trial_dfof(21+(5*(nOn+nOff))-10:21+(5*(nOn+nOff)),iC,ind_tar),3))));
%             end
%             if i== 1
%             legend('base1','base5',num2str(dirs(max_dir)-tars(1)),num2str(dirs(max_dir)-tars(2)));
%             end
%             ylim([-0.2 0.8])
            if dirs(max_dir(iC,:))<90
                min_dir(iC,:) = find(dirs == dirs(max_dir(iC,:))+90);
            else
                min_dir(iC,:) = find(dirs == dirs(max_dir(iC,:))-90);
            end
            plot(squeeze(data_trial_dir(11:35,iC,min_dir(iC,:))));
            ind_mindir = find(baseDir == dirs(min_dir(iC,:)));
            min_tars = unique(tarDir(ind_mindir));
            ind_mintar = intersect(find(tarDir == min_tars(2)),ind_mindir);
            plot(squeeze(mean(data_trial_dfof(21+(5*(nOn+nOff))-10:21+14+(5*(nOn+nOff)),iC,ind_mintar),3)- mean(mean(data_trial_dfof(21+(5*(nOn+nOff))-10:21+(5*(nOn+nOff)),iC,ind_mintar),3))));
            if i== 1
                legend('base1','base5','opp1','opp2target');
            end
            %ylim([-0.2 0.8])
            title(num2str(chop(osi(iC,:),2)))
        end
        
        data_dir_temp = data_dir_avg;
        data_dir_temp(find(data_dir_avg<0)) = 0;   
        osi = nans(nCells,1);
        for iC = good_ind
            osi(iC,:) = (data_dir_temp(iC,max_dir(iC),1)-data_dir_temp(iC,min_dir(iC),1))./(data_dir_temp(iC,max_dir(iC),1)+data_dir_temp(iC,min_dir(iC),1));
        end

        data_dir_norm = bsxfun(@rdivide, data_dir_avg, max(data_dir_avg(:,:,1),[],2));
        data_dir_pp = data_dir_avg(:,:,2)./data_dir_avg(:,:,1);
        data_dir_ad = data_dir_avg(:,:,5)./data_dir_avg(:,:,1);
        data_dir_oad = circshift(data_dir_avg(:,:,6),3,2)./data_dir_avg(:,:,1);
        
        
        
        figure;
        temp = nan(length(good_ind),1);
        for iC = good_ind
            temp(iC,:) = data_dir_pp(iC,(osi(iC,1)/30)+1);
        end
        subplot(3,2,1)
        hist(temp)
        xlabel('Paired-pulse ratio')
        xlim([-0.5 2])
        subplot(3,2,2)
        scatter(temp,osi(:,2))
        ylabel('OSI')
        xlabel('Paired-pulse adaptation')
        xlim([-0.5 2])
        
        temp = nan(length(good_ind),1);
        for iC = good_ind
            temp(iC,:) = data_dir_ad(iC,(osi(iC,1)/30)+1);
        end
        subplot(3,2,3)
        hist(temp)
        xlabel('Steady-state adaptation')
        xlim([-0.5 2])
        subplot(3,2,4)
        scatter(temp,osi(:,2))
        ylabel('OSI')
        xlabel('Steady-state adaptation')
        xlim([-0.5 2])
        
        temp = nan(length(good_ind),1);
        for iC = good_ind
            temp(iC,:) = data_dir_oad(iC,(osi(iC,1)/30)+1);
        end
        subplot(3,2,5)
        hist(temp)
        xlabel('Orientation selective adaptation')
        xlim([-0.5 2])
        subplot(3,2,6)
        scatter(temp,osi(:,2))
        ylabel('OSI')
        xlabel('Orientation selective adaptation')
        xlim([-0.5 2])
    else
        if iscell(input.nFramesOn)
            nOn = input.nFramesOn{1};
        else
            nOn = input.nFramesOn;
        end
        cStart = cell2mat(input.cFirstStim);
        nTrials = length(tCyc);
        nCells = size(npSub_tc,2);
        maxCyc = max(tCyc,[],2);
        data_trial = nan(40,nCells,maxCyc+1,nTrials);
        tFramesOff = nan(nTrials,maxCyc);
        SIx = strcmp(input.trialOutcomeCell, 'success');
        MIx = strcmp(input.trialOutcomeCell, 'ignore');
        FIx = strcmp(input.trialOutcomeCell, 'failure');
        nCyc = tCyc;
        nCyc([find(MIx) find(SIx)]) = tCyc([find(MIx) find(SIx)])+1;
        for itrial = 1:nTrials
            if isfield(input, 'tFramesOff')
                if length(input.tFramesOff{itrial}>0)
                    tempFramesOff = input.tFramesOff{itrial};
                else
                    tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
                    input.tFramesOff{itrial} = tempFramesOff;
                end
            else
                if iscell(input.nFramesOff)
                    tempFramesOff = input.nFramesOff{itrial}.*(ones(1,tCyc(itrial)));
                else
                    tempFramesOff = input.nFramesOff.*(ones(1,tCyc(itrial)));
                end
            end
                
            tFramesOff(itrial,1:tCyc(itrial)) = tempFramesOff(1:tCyc(itrial));
            for icyc = 1:nCyc(itrial)
                if icyc > 1
                    cyc_add = ((icyc-1)*nOn)+sum(tempFramesOff(1:icyc-1))+1;
                else
                    cyc_add = 0;
                end
                data_trial(:,:,icyc,itrial) = npSub_tc(cStart(itrial)-20+cyc_add:cStart(itrial)+19+cyc_add,:);
            end
        end
        data_f = mean(data_trial(1:20,:,1,:),1);
        data_dfof = bsxfun(@rdivide,bsxfun(@minus,data_trial,data_f),data_f);
        baseDir = celleqel2mat_padded(input.tBaseGratingDirectionDeg);
        dirs = unique(baseDir);
        ndir = length(dirs);
        data_dfof_dir = zeros(size(data_trial,1),nCells,maxCyc+1,ndir);
        tGratingDir = round(cell2mat(input.tGratingDirectionDeg),0);
        targetDelta = tGratingDir-baseDir;
        deltas = unique(targetDelta);
        nDelta = length(deltas);
        offs = unique(tFramesOff(:,1));
        noff = length(offs);
        frameRateHz = input.frameRateHz;
        
        base_win = [21:23];
        resp_win = [28:30];

        save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']),'baseDir', 'dirs', 'ndir', 'tFramesOff', 'offs', 'noff', 'tGratingDir', 'targetDelta', 'deltas', 'nDelta', 'tCyc', 'nCyc', 'maxCyc', 'nCells', 'frameRateHz', 'nTrials', 'SIx', 'MIx', 'FIx', 'cTarget', 'cStart')
        %find good cells- responsive to all base or at least one base direction
        [x,y] = ttest(squeeze(mean(data_dfof(base_win,:,1,:),1))',squeeze(mean(data_dfof(resp_win,:,1,:),1))','tail','left','alpha',0.05);
        h = zeros(ndir,nCells);
        p = zeros(ndir,nCells);
        

        if ndir > 1
            for idir = 1:ndir
                ind = find(baseDir == dirs(idir));
                data_dfof_dir = nan(40, nCells, maxCyc, ndir);
                for icyc = 1:maxCyc
                    ind2 = intersect(ind,find(nCyc >= icyc));
                    data_dfof_dir(:,:,icyc,idir)= nanmean(nanmean(data_dfof(:,:,icyc,ind2),3),4);
                end
                for iCell = 1:nCells
                    [h(idir,iCell), p(idir,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,1,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,1,ind),1)),'tail','left','alpha',0.05./(ndir-1));
                end
            end
            good_ind = unique([find(x) find(sum(h,1)>0)]);
        else
            data_dfof_dir = nan(40, nCells, maxCyc);
            ind_cyc = zeros(1,maxCyc);
            for icyc = 1:maxCyc
                if icyc == 1
                    %ind = intersect(find(SIx),find(nCyc>icyc));
                    ind = find(tCyc>icyc);
                    ind_cyc(1,icyc)= length(ind);
                else
                    %ind = intersect(find(SIx), find(nCyc >= icyc));
                    ind = find(tCyc >= icyc);
                    ind_cyc(1,icyc)= length(ind);
                end
                data_dfof_dir(:,:,icyc)= squeeze(nanmean(data_dfof(:,:,icyc,ind),4));
            end
            good_ind = find(x);
        end
        
        %plot ori tuning of all good cells
        if ndir > 1
            figure; 
            [n n2] = subplotn(length(good_ind));
            data_dir_avg = zeros(nCells,ndir);
            max_dir = NaN(nCells,1);
            for i = 1:length(good_ind)
                subplot(n,n2,i); 
                iC = good_ind(i);
                for idir = 1:ndir
                    ind = find(baseDir == dirs(idir));
                    data_dir_avg(iC,idir) = squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1));
                    if h(idir,iC)
                        errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'or');
                    else
                        errorbar(dirs(idir),squeeze(mean(data_dfof_dir(resp_win,iC,1,idir),1)),squeeze(std(mean(data_dfof(resp_win,iC,1,ind),1),[],4)./sqrt(length(ind))), 'ok');
                    end
                    hold on
                    h_ind = find(h(:,iC));
                    if length(h_ind)>0
                        [max_val max_ind] = max(data_dir_avg(iC,h_ind),[],2);
                        max_dir(iC,:) = h_ind(max_ind);
                    else
                        [max_val max_ind] = max(data_dir_avg(iC,:),[],2);
                        max_dir(iC,:) = max_ind;
                    end
                end
                title(['Cell ' num2str(iC) '- Pref ' num2str(dirs(max_dir(iC,:)))])
                ylim([0 max(mean(data_dfof_dir(31:35,iC,:),1),[],3)*1.5])
                xlim([-dirs(2) dirs(end)+dirs(2)])
            end
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
        elseif ~doFromRef
            max_dir = ones(nCells,1);
        end
        %response by cycle- avg all cells
        tt = [-20:19].*frameRateHz;
        data_dfof_base = zeros(40,nCells,maxCyc);
        if ~doFromRef
            for i = 1:length(good_ind)
                iC = good_ind(i);
                data_dfof_base(:,iC,:) = (data_dfof_dir(:,iC,:,max_dir(iC,:)));
            end
        else
            %data_dfof_base = data_dfof_dir(:,:,1:maxCyc+1);
        end
        
        figure;
        subplot(2,1,1)
        for icyc = 1:maxCyc
            plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,icyc),nanmean(data_dfof_base(base_win,good_ind,icyc),1)),2))
            hold on
        end
        ylim([-0.02 0.15])
        xlabel('Time (ms)')
        legend([repmat('cycle ',[maxCyc 1]) num2str([1:maxCyc]')])
        legend('Location','Northwest')
        subplot(2,1,2)
        data_dfof_base_diff = zeros(maxCyc,length(good_ind));
        for icyc = 1:maxCyc
            data_dfof_base_diff(icyc,:) = nanmean(data_dfof_base(resp_win,good_ind,icyc),1)- nanmean(data_dfof_base(base_win,good_ind,icyc),1);
            errorbar(icyc,nanmean(data_dfof_base_diff(icyc,:),2),nanstd(data_dfof_base_diff(icyc,:),[],2)./sqrt(length(good_ind)),'o')
            hold on
        end
        ylim([0 0.15])
        ylabel('dF/F')
        xlabel('cycle #')
        set(gca,'XTick',1:maxCyc)
        suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells- ' num2str(ind_cyc)]) 
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc.pdf']),'-dpdf','-bestfit')
       
        figure
        [n n2] = subplotn(length(good_ind));
        for iCell = 1:length(good_ind)
            iC = good_ind(iCell);
            subplot(n,n2,iCell)
            for icyc = 1:2
                plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,iC,icyc),nanmean(data_dfof_base(base_win,iC,icyc),1)),2))
                hold on
            end
        end
        suptitle([mouse ' ' date])
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_2Cyc_byCell.pdf']),'-dpdf','-bestfit')
        
        off_str{1} = 'base';
        for ioff = 1:noff
            off_str{ioff+1} = num2str(offs(ioff).*frameRateHz);
        end
        %response by off frames, avg all cycles all cells
        data_dfof_off = zeros(40,nCells,noff);
        for i = 1:length(good_ind)
            iC = good_ind(i);
            for ioff = 1:noff
                temp = [];
                for icyc = 1:maxCyc-1
                    if ~doFromRef
                        ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
                    else
                        ind = find(tFramesOff(:,icyc) == offs(ioff));
                    end
                    temp = [temp squeeze(data_dfof(:,iC,icyc+1,ind))];
                end
                data_dfof_off(:,iC,ioff) = nanmean(temp,2);
            end
        end
        figure;
        subplot(2,1,1)
        plot(tt,nanmean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
        hold on
        for ioff = 1:noff
            plot(tt,nanmean(bsxfun(@minus,data_dfof_off(:,good_ind,ioff), nanmean(data_dfof_off(base_win,good_ind,ioff),1)),2));
        end
        legend(off_str)
        legend('Location','Northwest')
        ylim([-0.02 0.15])
        xlabel('Time (ms)')
        subplot(2,1,2)
        errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
        hold on
        data_dfof_base_off = zeros(noff,length(good_ind));
        for ioff = 1:noff
            temp = nanmean(data_dfof_off(28:31,good_ind,ioff),1)- nanmean(data_dfof_off(20:23,good_ind,ioff),1);
            errorbar((offs(ioff)*frameRateHz)./1000,mean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
            hold on
        end
        set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
        xlim([0.05 10])
        ylim([0 0.15])
        ylabel('dF/F')
        xlabel('Off Interval (ms)')
        suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt.pdf']),'-dpdf', '-bestfit')
        
        %response by off frames by cell, avg all cycles
        figure;
        [n n2] = subplotn(length(good_ind));
        for iCell = 1:length(good_ind)
            iC = good_ind(iCell);
            subplot(n, n2, iCell)
            plot(tt,bsxfun(@minus,data_dfof_base(:,iC,1),nanmean(data_dfof_base(base_win,iC,1),1)))
            hold on
            for ioff = 1:noff
                plot(tt,bsxfun(@minus,data_dfof_off(:,iC,ioff), nanmean(data_dfof_off(base_win,iC,ioff),1)));
            end
        end
        suptitle([mouse ' ' date]) 
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byInt_byCell.pdf']),'-dpdf', '-bestfit')
        
        
        %response by off frames by cycle, avg all cells
        data_dfof_off_cyc = zeros(40,nCells,noff,maxCyc);
        ind_n = zeros(maxCyc,noff);
        for i = 1:length(good_ind)
            iC = good_ind(i);
            for ioff = 1:noff
                for icyc = 1:maxCyc
                    if ~doFromRef
                        ind = intersect(find(tCyc >= icyc+1),intersect(find(baseDir == dirs(max_dir(iC,:))), find(tFramesOff(:,icyc) == offs(ioff))));
                        ind_n(icyc,ioff) = length(ind);
                    else
                        ind = find(tFramesOff(:,icyc) == offs(ioff));
                    end
                    data_dfof_off_cyc(:,iC,ioff,icyc) = nanmean(data_dfof(:,iC,icyc+1,ind),4);
                end
            end
        end
        figure;
        [n n2] = subplotn(double(maxCyc-1));
        for icyc = 1:maxCyc-1
            subplot(n,n2,double(icyc))
            plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_ind,1),nanmean(data_dfof_base(base_win,good_ind,1),1)),2))
            hold on
            for ioff = 1:noff
                plot(tt,mean(bsxfun(@minus,data_dfof_off_cyc(:,good_ind,ioff,icyc),mean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1)),2));
            end
            legend(off_str,'Location','Northwest')
            title(['Cycle ' num2str(icyc+1) '- ' num2str(ind_n(icyc,:)) ' trials'])
            ylim([-0.02 0.15])
            xlabel('Time (ms)')
        end
        suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt.pdf']),'-dpdf', '-bestfit')
        
        figure;
        for icyc = 1:(maxCyc-1)
            subplot(n,n2,double(icyc))
            errorbar(8,nanmean(data_dfof_base_diff(1,:),2),nanstd(data_dfof_base_diff(1,:),[],2)./sqrt(length(good_ind)),'o')
            hold on
            for ioff = 1:noff
                temp = bsxfun(@minus, nanmean(data_dfof_off_cyc(resp_win,good_ind,ioff,icyc),1), nanmean(data_dfof_off_cyc(base_win,good_ind,ioff,icyc),1));
                errorbar((offs(ioff)*frameRateHz)./1000,nanmean(temp,2),nanstd(temp,[],2)./sqrt(length(good_ind)),'o')
                hold on
            end
            set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
            xlim([0.05 10])
            ylim([-0.02 0.15])
            ylabel('dF/F')
            xlabel('Off Interval (s)')
            title(['Cycle ' num2str(icyc+1)])
        end
        suptitle([mouse ' ' date '- ' num2str(length(good_ind)) ' cells']) 
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byCyc_byInt_quant.pdf']),'-dpdf', '-bestfit')
               
        %find cells responsive to target or at least one target direction
        data_dfof_targ = nan(40,nCells,nTrials);
        tFramesOff_targ = nan(1,nTrials);
        for itrial = 1:nTrials
            if ~isnan(cTarget(itrial))
                icyc = nCyc(itrial);
                data_dfof_targ(:,:,itrial) = data_dfof(:,:,icyc,itrial);
                tFramesOff_targ(1,itrial) = tFramesOff(itrial,icyc-1);
            end
        end
        [x_targ y_targ] = ttest(squeeze(nanmean(data_dfof_targ(base_win,:,:),1))',squeeze(nanmean(data_dfof_targ(resp_win,:,:),1))','tail','left','alpha',0.05);
        
        if ndir > 1
            h_targ = zeros(ndir,nCells);
            p_targ = zeros(ndir,nCells);
            for idir = 1:ndir
                ind1 = intersect(find(baseDir == dirs(idir)),find(targetDelta== deltas(1)));
                shift = idir-2;
                if shift<1
                    shift = shift+input.baseGratingDirectionStepN;
                end
                ind2 = intersect(find(baseDir == dirs(shift)),find(targetDelta== deltas(2)));
                ind = [ind1 ind2];
                for iCell = 1:nCells
                    [h_targ(idir,iCell), p_targ(idir,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof_targ(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
                end
            end
        else
            h_targ = zeros(nDelta,nCells);
            p_targ = zeros(nDelta,nCells);
            for itarg = 1:nDelta
                ind = find(targetDelta== deltas(itarg));
                for iCell = 1:nCells
                    [h_targ(itarg,iCell), p_targ(itarg,iCell)] = ttest(squeeze(nanmean(data_dfof_targ(base_win,iCell,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
                end
            end
        end
        good_targ_ind = intersect(good_ind, unique([find(sum(h_targ,1)>0) find(x_targ)]));
        good_resp_ind = unique([find(sum(h_targ,1)>0) find(x_targ) good_ind]);
        
        
        save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']), 'data_dfof', 'max_dir','good_ind', 'good_targ_ind', 'good_resp_ind')

        if ndir > 1
            %figure for each cell with base, target and target by change
            for iCell = 1:nCells
                figure;
                for idir = 1:ndir
                    ind = find(baseDir == dirs(idir));
                    subplot(2,2,1)
                    plot(mean(data_dfof(:,iCell,1,ind),4))
                    hold on
                    subplot(2,2,2)
                    ind1 = intersect(ind,find(targetDelta== deltas(1)));
                    shift = idir-2;
                    if shift<1
                        shift = shift+input.baseGratingDirectionStepN;
                    end
                    ind2 = intersect(find(baseDir == dirs(shift)),find(targetDelta== deltas(2)));
                    ind_t = [ind1 ind2];
                    plot(mean(data_dfof(:,iCell,6,ind_t),4)-mean(mean(data_dfof(base_win,iCell,6,ind_t),1),4))
                    hold on
                    for itarg = 1:length(deltas)
                        ind3 = intersect(ind, find(targetDelta == deltas(itarg)));
                        subplot(2,2,2+itarg)
                        plot(mean(data_dfof(:,iCell,6,ind3),4)-mean(mean(data_dfof(base_win,iCell,6,ind3),1),4))
                        hold on
                    end
                    subplot(2,2,1)
                    title(['Base- ' num2str(dirs(find(h(:,iCell))))])
                    subplot(2,2,2)
                    title(['Target- ' num2str(dirs(find(h_targ(:,iCell))))])
                    subplot(2,2,3)
                    title(['Target- ' num2str(deltas(1)) ' deg change'])
                    subplot(2,2,4)
                    title(['Target- ' num2str(deltas(2)) ' deg change'])
                end
                base = num2str(length(find(good_ind == iCell)));
                if base == '1'
                    pref = num2str(dirs(max_dir(iCell,:)));
                    targ =  num2str(length(find(good_targ_ind == iCell)));
                else
                    pref = 'NaN';
                    targ = 'NaN';
                end
                suptitle(['Cell ' num2str(iCell) '- Pref is ' pref '- Base is ' base '- Target is ' targ])
                print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'AllCells',[date '_' mouse '_' run_str '_avgRespCell' num2str(iCell) '.pdf']),'-dpdf')
                close all
            end

            data_dfof_targ = zeros(40,nCells,ndir,nDelta);
            for idir = 1:ndir
                ind_dir = find(baseDir == dirs(idir));
                for itarg = 1:nDelta
                    ind_targ = intersect(ind_dir, find(targetDelta == deltas(itarg)));
                    if length(ind_targ)>4
                        data_dfof_targ(:,:,idir,itarg) = squeeze(mean(data_dfof(:,:,end,ind_targ),4));
                    else
                        data_dfof_targ(:,:,idir,itarg) = NaN(40,1);
                    end
                end
            end

            %average target response for absolute preferred stim
            data_dfof_delta = zeros(40,sz(2),nDelta);
            for i = 1:length(good_targ_ind)
                iC = good_targ_ind(i);
                idir = max_dir(iC,:);
                shift1 = idir-1;
                if shift1<1
                    shift1 = shift1+ndir;
                end
                ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));    
                data_dfof_delta(:,iC,1) = mean(data_dfof(:,iC,6,ind1),4);
                shift2 = idir-3;
                if shift2<1
                    shift2 = shift2+ndir;
                end
                ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));    
                data_dfof_delta(:,iC,2) = mean(data_dfof(:,iC,6,ind2),4);
            end
            figure;
            subplot(2,1,1)
            plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_targ_ind,1),mean(data_dfof_base(base_win,good_targ_ind,1),1)),2))
            hold on
            for itarg = 1:nDelta
                plot(tt, mean(bsxfun(@minus,data_dfof_delta(:,good_targ_ind,itarg),mean(data_dfof_delta(base_win,good_targ_ind,itarg),1)),2))
                hold on
            end
            legend('base', [num2str(deltas(1)) ' deg target'], [num2str(deltas(2)) ' deg target'])
            legend('Location','Northwest')
            xlabel('Time (ms)')
            ylabel('dF/F')
            subplot(2,1,2)
            temp = mean(data_dfof_base(resp_win,good_targ_ind,1),1)-mean(data_dfof_base(base_win,good_targ_ind,1),1);
            errorbar(0,mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
            hold on
            for itarg = 1:nDelta
                temp  = mean(data_dfof_delta(resp_win,good_targ_ind,itarg),1)-mean(data_dfof_delta(base_win,good_targ_ind,itarg),1);
                errorbar(deltas(itarg), mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
                hold on
            end
            xlabel('Degree change')
            ylabel('dF/F')
            ylim([0 0.1])
            xlim([-10 100])
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta.pdf']),'-dpdf')

            %tuning curves for baseline
            data_dfof_base_tuning = zeros(40,sz(2),ndir,nCyc);
            for i = 1:length(good_ind)
                iC = good_ind(i);
                for idir = 1:ndir
                    dir_temp = max_dir(iC,:)+(idir-1);
                    if dir_temp>ndir
                        dir_temp = rem(dir_temp,ndir);
                    end
                    ind = find(baseDir == dirs(dir_temp));
                    for icyc = 1:nCyc
                        data_dfof_base_tuning(:,iC,idir,icyc) = mean(data_dfof(:,iC,icyc,ind),4);
                    end
                end
            end
            %for base responsive cells
            figure;
            for icyc = 1:nCyc
                for idir = 1:ndir
                    subplot(2,double(nCyc),double(icyc))
                    plot(tt,mean(bsxfun(@minus,data_dfof_base_tuning(:,good_ind,idir,icyc),mean(data_dfof_base_tuning(base_win,good_ind,idir,icyc),1)),2))
                    hold on
                    subplot(2,double(nCyc),double(icyc+nCyc))
                    temp = bsxfun(@minus,mean(data_dfof_base_tuning(resp_win,good_ind,idir,icyc),1),mean(data_dfof_base_tuning(base_win,good_ind,idir,icyc),1));
                    errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_ind)),'o');
                    hold on
                end
                ylim([0 0.1])
                xlim([-30 180])
                xlabel('Ori (deg)')
                subplot(2,double(nCyc),double(icyc))
                ylim([-0.02 0.1])
                title(['Cycle # ' num2str(icyc)])
                xlabel('Time (ms)')
            end
            suptitle('Baseline responsive cells')
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_BaseOriTuning_goodCells.pdf']),'-dpdf')

            %for target responsive cells
            figure;
            for icyc = 1:nCyc
                for idir = 1:ndir
                    subplot(2,double(nCyc),double(icyc))
                    plot(tt,mean(bsxfun(@minus,data_dfof_base_tuning(:,good_targ_ind,idir,icyc),mean(data_dfof_base_tuning(base_win,good_targ_ind,idir,icyc),1)),2))
                    hold on
                    subplot(2,double(nCyc),double(icyc+nCyc))
                    temp = bsxfun(@minus,mean(data_dfof_base_tuning(resp_win,good_targ_ind,idir,icyc),1),mean(data_dfof_base_tuning(base_win,good_targ_ind,idir,icyc),1));
                    errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
                    hold on
                end
                ylim([0 0.1])
                xlim([-30 180])
                xlabel('Ori (deg)')
                subplot(2,double(nCyc),double(icyc))
                ylim([-0.02 0.1])
                title(['Cycle # ' num2str(icyc)])
                xlabel('Time (ms)')
            end
            suptitle('Target responsive cells')
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_BaseOriTuning_goodTargCells.pdf']),'-dpdf')

            %tuning curves for targets
            data_dfof_targ_tuning = zeros(40,sz(2),ndir,nDelta+4);
            for i = 1:length(good_targ_ind)
                iC = good_targ_ind(i);
                for idir = 1:ndir
                    dir_temp = max_dir(iC,:)+(idir-1);
                    if dir_temp>ndir
                        dir_temp = rem(dir_temp,ndir);
                    end
                    ind = find(baseDir == dirs(dir_temp));
                    data_dfof_targ_tuning(:,iC,idir,1) = mean(data_dfof(:,iC,1,ind),4);
                    data_dfof_targ_tuning(:,iC,idir,2) = mean(mean(data_dfof(:,iC,3:5,ind),3),4);
                    shift1 = max_dir(iC,:)-1+(idir-1);
                    if shift1<1
                        shift1 = shift1+ndir;
                    elseif shift1>ndir
                        shift1 = rem(shift1,ndir);
                    end
                    ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));    
                    data_dfof_targ_tuning(:,iC,idir,3) = mean(data_dfof(:,iC,6,ind1),4);
                    shift2 = max_dir(iC,:)-1+(idir-3);
                    if shift2<1
                        shift2 = shift2+ndir;
                    elseif shift2>ndir
                        shift2 = rem(shift2,ndir);
                    end
                    ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));    
                    data_dfof_targ_tuning(:,iC,idir,4) = mean(data_dfof(:,iC,6,ind2),4);
                end
            end
            figure;
            for i = 1:4
                for idir = 1:ndir
                    subplot(2,4,i)
                    plot(tt,mean(bsxfun(@minus,data_dfof_targ_tuning(:,good_targ_ind,idir,i),mean(data_dfof_targ_tuning(base_win,good_targ_ind,idir,i),1)),2))
                    hold on
                    subplot(2,4,i+4)
                    temp = bsxfun(@minus,mean(data_dfof_targ_tuning(resp_win,good_targ_ind,idir,i),1),mean(data_dfof_targ_tuning(base_win,good_targ_ind,idir,i),1));
                    errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
                    hold on
                end
                ylim([0 0.1])
                ylabel('dF/F')
                xlim([-30 180])
                xlabel('Ori (deg)')
                subplot(2,4,i)
                ylim([-0.02 0.1])
                ylabel('dF/F')
                xlabel('Time (ms)')
            end
            subplot(2,4,1)
            title('Base 1')
            subplot(2,4,2)
            title('Base 3-5')
            subplot(2,4,3)
            title(['Target- ' num2str(deltas(1)) ' deg change'])
            subplot(2,4,4)
            title(['Target- ' num2str(deltas(2)) ' deg change'])
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_OriTuning.pdf']),'-dpdf')

            %tuning curves for targets by interval
            data_dfof_targ_tuning_off = zeros(40,sz(2),ndir,nDelta,noff);
            for i = 1:length(good_targ_ind)
                iC = good_targ_ind(i);
                for idir = 1:ndir
                    shift1 = max_dir(iC,:)-1+(idir-1);
                    if shift1<1
                        shift1 = shift1+ndir;
                    elseif shift1>ndir
                        shift1 = rem(shift1,ndir);
                    end
                    shift2 = max_dir(iC,:)-1+(idir-3);
                    if shift2<1
                        shift2 = shift2+ndir;
                    elseif shift2>ndir
                        shift2 = rem(shift2,ndir);
                    end
                    for ioff = 1:noff
                        ind1 = intersect(find(tFramesOff(:,5) == offs(ioff)),intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1))));    
                        data_dfof_targ_tuning_off(:,iC,idir,1,ioff) = mean(data_dfof(:,iC,6,ind1),4);
                        ind2 = intersect(find(tFramesOff(:,5) == offs(ioff)),intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2))));    
                        data_dfof_targ_tuning_off(:,iC,idir,2,ioff) = mean(data_dfof(:,iC,6,ind2),4);
                    end
                end
            end

            figure; 
            for itarg = 1:nDelta
                for ioff = 1:noff
                    subplot(2,3,ioff +((itarg-1)*noff))
                    for idir = 1:ndir
                        temp = bsxfun(@minus,mean(data_dfof_targ_tuning_off(resp_win,good_targ_ind,idir,itarg,ioff),1),mean(data_dfof_targ_tuning_off(base_win,good_targ_ind,idir,itarg,ioff),1));
                        errorbar(dirs(idir),mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o');
                        hold on
                    end
                    xlim([-30 180])
                    xlabel('Ori (deg)')
                    ylim([-0.04 0.12])
                    ylabel('dF/F')
                    hline(0, '--k')
                    title([num2str(offs(ioff)*frameRateHz) ' ms; ' num2str(deltas(itarg)) ' deg change'])
                end
            end
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_OriTuning_byInt.pdf']),'-dpdf')

            %target response (for absolute pref stim) by interval
            data_dfof_delta_off = zeros(40,sz(2),nDelta,noff);
            for i = 1:length(good_targ_ind)
                iC = good_targ_ind(i);
                idir = max_dir(iC,:);
                shift1 = idir-1;
                if shift1<1
                    shift1 = shift1+ndir;
                end
                ind1 = intersect(find(baseDir == dirs(shift1)),find(targetDelta == deltas(1)));
                for ioff = 1:noff
                    ind1_off = intersect(ind1, find(tFramesOff(:,5) == offs(ioff)));
                    if length(ind1_off)>4
                        data_dfof_delta_off(:,iC,1,ioff) = mean(data_dfof(:,iC,6,ind1_off),4);
                    else
                        data_dfof_delta_off(:,iC,1,ioff) = NaN(40,1);
                    end
                end
                shift2 = idir-3;
                if shift2<1
                    shift2 = shift2+ndir;
                end
                ind2 = intersect(find(baseDir == dirs(shift2)),find(targetDelta == deltas(2)));
                for ioff = 1:noff
                    ind2_off = intersect(ind2, find(tFramesOff(:,5) == offs(ioff)));
                    if length(ind2_off)>4
                        data_dfof_delta_off(:,iC,2,ioff) = mean(data_dfof(:,iC,6,ind2_off),4);
                    else
                        data_dfof_delta_off(:,iC,2,ioff) = NaN(40,1);
                    end
                end
            end
            figure;
            for itarg = 1:nDelta
                subplot(2,2,itarg)
                plot(tt,mean(bsxfun(@minus,data_dfof_base(:,good_targ_ind,1),mean(data_dfof_base(base_win,good_targ_ind,1),1)),2))
                hold on
                for ioff = 1:noff
                    plot(tt, nanmean(bsxfun(@minus, data_dfof_delta_off(:,good_targ_ind,itarg,ioff),nanmean(data_dfof_delta_off(base_win,good_targ_ind,itarg,ioff),1)),2))
                    hold on
                end
                xlabel('Time (ms)')
                title([num2str(deltas(itarg)) ' deg target'])
                subplot(2,2,itarg+2)
                temp = mean(data_dfof_base(resp_win,good_targ_ind,1),1)-mean(data_dfof_base(base_win,good_targ_ind,1),1);
                errorbar(8,mean(temp,2),std(temp,[],2)./sqrt(length(good_targ_ind)),'o')
                hold on
                for ioff = 1:noff
                    temp  = bsxfun(@minus, nanmean(data_dfof_delta_off(resp_win,good_targ_ind,itarg,ioff),1),nanmean(data_dfof_delta_off(base_win,good_targ_ind,itarg,ioff),1));
                    errorbar((offs(ioff)*frameRateHz)./1000, nanmean(temp,2),nanstd(temp,[],2)./sqrt(length(~isnan(temp))),'o')
                    hold on
                end
                set(gca,'xscale','log','XTick', [0.1 0.2 0.4 0.8, 1.6 3.2 6.4 12.8])
                xlim([0.05 10])
                ylim([-0.02 0.1])
                ylabel('dF/F')
                xlabel('Off Interval (s)')
            end
            subplot(2,2,2)
            legend('base', num2str(offs(1).*frameRateHz),num2str(offs(2).*frameRateHz),num2str(offs(3).*frameRateHz))
            legend('Location','Northwest')
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byInt.pdf']),'-dpdf')
        else
            figure;
            [n n2] = subplotn(length(good_resp_ind));
            for iCell = 1:length(good_resp_ind)
                iC = good_resp_ind(iCell);
                subplot(n,n2,iCell)
                plot(tt,data_dfof_base(:,iC,1))
                hold on
                for itarg = 1:nDelta
                    ind = find(targetDelta == deltas(itarg));
                    plot(tt,nanmean(data_dfof_targ(:,iC,ind),3)-nanmean(nanmean(data_dfof_targ(base_win,iC,ind),1),3));
                    hold on;
                end
            end
            suptitle([mouse ' ' date]) 
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byCell.pdf']),'-dpdf','-bestfit')

            targ_only_ind = unique([find(sum(h_targ,1)>0) find(x_targ)]);
            indn = zeros(nDelta, noff);
            figure;
            for itarg = 1:nDelta
                for ioff = 1:noff
                    ind = intersect(find(targetDelta == deltas(itarg)), find(tFramesOff_targ== offs(ioff)));
                    indn(itarg,ioff) = length(ind);
                    subplot(2,nDelta,itarg)
                    plot(tt, nanmean(nanmean(bsxfun(@minus,data_dfof_targ(:,targ_only_ind,ind), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3),2))
                    hold on
                    subplot(2,nDelta,itarg+nDelta)
                    temp = nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,targ_only_ind,ind),1), nanmean(data_dfof_targ(base_win,targ_only_ind,ind),1)),3);
                    errorbar(offs(ioff)*frameRateHz, nanmean(temp,2), nanstd(temp,[],2)./sqrt(length(targ_only_ind)),'o')
                    hold on
                end
                ylim([0 0.2])
                ylabel('dF/F')
                xlabel('Off interval (ms)')
                subplot(2,nDelta,itarg)
                ylim([-0.05 0.2])
                ylabel('dF/F')
                xlabel('Time (ms)')
                title([num2str(deltas(itarg)) ' deg change- ' num2str(indn(itarg,:))])
            end
            suptitle([mouse ' ' date '- ' num2str(length(targ_only_ind)) ' cells']) 
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgResp_byTargDelta_byInt_targInd.pdf']),'-dpdf','-bestfit')
                
            figure; 
            [n n2] = subplotn(length(good_resp_ind));
            dirs = [0:30:150];
            for iCell = 1:length(good_resp_ind)
                iC = good_resp_ind(iCell);
                subplot(n, n2, iCell)
                errorbar(dirs(1),mean(bsxfun(@minus,nanmean(data_dfof(resp_win,iC,1,:),1),nanmean(data_dfof(base_win,iC,1,:),1)),4), nanstd(bsxfun(@minus,nanmean(data_dfof(resp_win,iC,1,:),1),nanmean(data_dfof(base_win,iC,1,:),1)),[],4)./sqrt(nTrials), 'o')
                hold on
                %errorbar(dirs(1),nanmean(bsxfun(@minus,mean(data_dfof_base(resp_win,iC,5,:),1),mean(data_dfof(base_win,iC,5,:),1)),4), std(bsxfun(@minus,mean(data_dfof(resp_win,iC,5,:),1),mean(data_dfof(base_win,iC,5,:),1)),[],4)./sqrt(nTrials), 'o')
                for itarg = 1:nDelta
                    ind = find(targetDelta == deltas(itarg));
                    errorbar(deltas(itarg),nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),3), nanstd(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),[],3)./sqrt(length(ind)), 'o')
                end
                xlim([-30 100])
                if length(unique(max_dir)) > 1
                    title(num2str(dirs(max_dir(iC,:))))
                end
            end
            
            baseV90 = zeros(2,nCells);
            ind = find(targetDelta == deltas(itarg));
            for iCell = 1:length(good_targ_ind)
                iC = good_targ_ind(iCell);
                baseV90(:,iC) = [mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,1,:),1),mean(data_dfof(base_win,iC,1,:),1)),4) mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,6,ind),1),mean(data_dfof(base_win,iC,6,ind),1)),4)];
            end
            
            
            for i = 1:length(good_targ_ind)
                figure;
                iC = good_targ_ind(i);
                for ioff = 1:noff
                    subplot(2, noff, ioff)
                    ind = find(tFramesOff(:,5) == offs(ioff));
                    plot(tt,mean(mean(bsxfun(@minus,data_dfof(:,iC,5,ind),mean(data_dfof(base_win,iC,5,ind),1)),4),2))
                    hold on
                    subplot(2, noff, ioff+noff)
                    errorbar(dirs(1),mean(mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,5,ind),1),mean(data_dfof(base_win,iC,5,ind),1)),4),2), std(mean(bsxfun(@minus,mean(data_dfof(resp_win,:,5,:),1),mean(data_dfof(base_win,:,5,:),1)),4),[],2)./sqrt(nCells), 'o')
                    hold on
                    for itarg = 1:nDelta
                        subplot(2, noff, ioff)
                        ind = intersect(find(tFramesOff(:,5) == offs(ioff)), find(targetDelta == deltas(itarg)));
                        plot(tt,mean(mean(bsxfun(@minus,data_dfof(:,iC,6,ind),mean(data_dfof(base_win,iC,6,ind),1)),4),2))
                        hold on
                        subplot(2, noff, ioff+noff)
                        errorbar(deltas(itarg),mean(mean(bsxfun(@minus,mean(data_dfof(resp_win,iC,6,ind),1),mean(data_dfof(base_win,iC,6,ind),1)),4),2), std(mean(bsxfun(@minus,mean(data_dfof(resp_win,:,6,ind),1),mean(data_dfof(base_win,:,6,ind),1)),4),[],2)./sqrt(nCells), 'o')
                        hold on
                    end
                    ylim([-0.02 0.3])
                    xlim([-30 100])
                    subplot(2, noff, ioff)
                    title([num2str(offs(ioff)*frameRateHz) ' ms interval'])
                    ylim([-0.02 0.3])
                end
                suptitle(['Cell #' num2str(iC) '; Pref- ' num2str(dirs(max_dir(iC,:))) ' deg'])
            end
            
            figure;
            delt_col = strvcat('k', 'b');
            [n n2] = subplotn(length(good_resp_ind));
            for iCell = 1:length(good_resp_ind)
                iC = good_resp_ind(iCell);
                subplot(n,n2,iCell)
                for ioff = 1:noff
                    for idelta = 1:nDelta
                        ind = intersect(find(tFramesOff_targ == offs(ioff)), find(targetDelta == deltas(idelta)));
                        errorbar(offs(ioff).*frameRateHz,nanmean(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),3),nanstd(bsxfun(@minus,nanmean(data_dfof_targ(resp_win,iC,ind),1),nanmean(data_dfof_targ(base_win,iC,ind),1)),[],3), ['o' delt_col(idelta)])
                        hold on
                        ylim([-0.05 0.2])
                    end
                end
            end
            
            
        end    
        ind_n = zeros(ndir,nDelta,noff);
        for idir = 1:ndir
            ind_dir = find(baseDir == dirs(idir));
            for itarg = 1:nDelta
                ind_targ = intersect(ind_dir, find(targetDelta == deltas(itarg)));
                for ioff = 1:noff
                    ind_off = intersect(ind_targ,find(tFramesOff_targ == offs(ioff)));
                    ind_n(idir,itarg,ioff) = length(ind_off);
                end
            end
        end
    end
end
clear data_mat data_f data_df

if input.doAnnulusSizeStim
    nStim = input.AnnulusgratingDiameterStepN;
    Stims = unique(celleqel2mat_padded(input.tAnnulusGratingDiameterDeg));
end

figure;
if nCells<37
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./frame_rate);
tuning_mat = zeros(nStim,2,nCells);
tc_mat = zeros(nOn+nOff,nStim,nCells);
start = 1;
f = 1;
cmap = flipud(gray(nStim+1));
Ind_struct = [];
for iCell = 1:nCells
    if start >36
        for i = 1:36
            subplot(n,n2,i)
            ylim([min(min(min(tc_mat(:,:,1:start-1),[],1),[],2),[],3) max(max(max(tc_mat(:,:,1:start-1),[],1),[],2),[],3)])
            vline(nOff)
        end
        suptitle(['Stims: ' num2str(chop(Stims,2))])
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        figure;
        f = 1+f;
        start = 1;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        if input.doDirStim
            ind = find(Dir == Dirs(iStim));
        elseif input.doTFStim
            ind = find(TF == TFs(iStim));
        elseif input.doSizeStim
            ind = find(sz_mat == Stims(iStim));
        elseif input.doRetStim
            ind1 = find(Az == Stims(iStim,2));
            ind2 = find(El == Stims(iStim,1));
            ind = intersect(ind1,ind2);
        end
        plot(tt',squeeze(mean(data_dfof(:,iCell,ind),3)), 'Color', cmap(iStim+1,:))
        hold on
        tc_mat(:,iStim,iCell) = squeeze(mean(data_dfof(:,iCell,ind),3));
        tuning_mat(iStim,1,iCell) = mean(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
        tuning_mat(iStim,2,iCell) = std(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
    end
    start = start+1;
end
for i = 1:start-1
    subplot(n,n2,i)
    ylim([min(min(min(tc_mat,[],1),[],2),[],3) max(max(max(tc_mat,[],1),[],2),[],3)])
end
if size(Stims,2) == 1
suptitle(['Stims: ' num2str(chop(Stims,2))])
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        figure;
        f = 1+f;
        start = 1;
    end
    subplot(n, n2, start)
    if input.doRetStim
        ret_mat = reshape(tuning_mat(:,1,iCell), [length(Els) length(Azs)]);
        if min(Els,[],2) <0
            ret_mat = flipud(ret_mat');
        else
            ret_mat = ret_mat';
        end
        imagesc(ret_mat)
        colormap gray
        clim([0 0.2])
        title(num2str(chop(max(tuning_mat(:,1,iCell),[],1),2)))        
    else    
        errorbar(Stims, tuning_mat(:,1,iCell), tuning_mat(:,2,iCell), '-ok');
        %ylim([-.02 0.2]) 
    end
    start = start + 1;
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning'  num2str(f) '.pdf']), '-dpdf')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TuningSummary.mat']), 'tuning_mat', 'tc_mat', 'Ind_struct')

%% 2AFC attention

tLeftTrial = celleqel2mat_padded(input.tLeftTrial);
nside = 2;
tGratingContrast = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tGratingContrast);
ncon = length(cons);
cStimOn = celleqel2mat_padded(input.cStimOn);
cDecision = celleqel2mat_padded(input.cDecision);
nTrials = length(tLeftTrial);
data_stim = nan(50,nCells,nTrials);
data_choice = nan(50,nCells,nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial)+29 < nframes
        data_stim(:,:,itrial) = npSub_tc(cStimOn(itrial)-20:cStimOn(itrial)+29,:);
    end
    if cDecision(itrial)
        if cDecision(itrial)+29 < nframes
            data_choice(:,:,itrial) = npSub_tc(cDecision(itrial)-20:cDecision(itrial)+29,:);
        end
    end
end
dataf = mean(data_stim(1:20,:,:),1);
data_stim_dfof = bsxfun(@rdivide, bsxfun(@minus, data_stim, dataf), dataf);
data_choice_dfof = bsxfun(@rdivide, bsxfun(@minus, data_choice, dataf), dataf);

ind = intersect(find(~tLeftTrial),find(tGratingContrast ==1));
figure; plot(mean(mean(data_stim_dfof(:,:,ind),2),3));

SIx = strcmp(input.trialOutcomeCell, 'success');
MIx = strcmp(input.trialOutcomeCell, 'ignore');
FIx = strcmp(input.trialOutcomeCell, 'incorrect');

h = zeros(nside, nCells);
p = zeros(nside, nCells);

for iS = 1:nside
    indS = find(tLeftTrial == iS-1);
    for icon = 1:ncon
        indC = intersect(indS, find(tGratingContrast == cons(icon)));
        if icon == ncon
            [h(iS,:) p(iS,:)] = ttest(squeeze(mean(data_stim_dfof(25:29,:,indC),1))', squeeze(mean(data_stim_dfof(19:23,:,indC),1))', 'tail', 'right');
        end
    end
end
good_ind = find(sum(h,1)>0);

tt = (-19:30)*frame_rate;


figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
       indS = intersect(find(SIx), find(tLeftTrial == iS-1));
       plot(tt,nanmean(data_stim_dfof(:,iCell,indS),3));
       hold on;
        if find(good_ind == iCell)
            good_str = ' resp';
        else
            good_str = ' not resp';
        end
        title(['Cell # ' num2str(iCell) ' is' good_str])
    end
end
suptitle([mouse ' ' date '- stimAlign- blue is right; red is left'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

figure;
[n n2] = subplotn(nCells);
for iCell = 1:nCells
    subplot(n, n2, iCell)
    for iS = 1:nside
       indS = intersect(find(SIx), find(tLeftTrial == iS-1));
       plot(tt,nanmean(data_choice_dfof(:,iCell,indS),3));
       hold on;
        if find(good_ind == iCell)
            good_str = ' resp';
        else
            good_str = ' not resp';
        end
        title(['Cell # ' num2str(iCell) ' is' good_str])
    end
end
suptitle([mouse ' ' date '- choiceAlign- blue is right; red is left'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byCell.pdf']), '-dpdf', '-bestfit')

figure; 
for icon = 1:ncon
    subplot(2,2,1)
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRH(1,icon) = length(ind1); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind1),2),3))
    hold on
    subplot(2,2,2)
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRM(1,icon) = length(ind2); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind2),2),3))
    hold on
    subplot(2,2,3)
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLH(1,icon) = length(ind3); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind3),2),3))
    hold on
    subplot(2,2,4)
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLM(1,icon) = length(ind4); 
    plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind4),2),3))
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
ylim([-.05 .2])
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
ylim([-.05 .2])
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
ylim([-.05 .2])
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
ylim([-.05 .2])
suptitle([mouse ' ' date '- stimAlign'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byOutcome.pdf']), '-dpdf', '-bestfit')

figure; 
for icon = 1:ncon
    subplot(2,2,1)
    ind1 = intersect(intersect(find(SIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRH(1,icon) = length(ind1); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind1),2),3))
    hold on
    subplot(2,2,2)
    ind2 = intersect(intersect(find(FIx), find(~tLeftTrial)), find(tGratingContrast == cons(icon)));
    indRM(1,icon) = length(ind2); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind2),2),3))
    hold on
    subplot(2,2,3)
    ind3 = intersect(intersect(find(SIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLH(1,icon) = length(ind3); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind3),2),3))
    hold on
    subplot(2,2,4)
    ind4 = intersect(intersect(find(FIx), find(tLeftTrial)), find(tGratingContrast == cons(icon)));
    indLM(1,icon) = length(ind4); 
    plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind4),2),3))
    hold on
end
subplot(2,2,1)
title(['Right Hit- ' num2str(indRH)])
ylim([-.05 .2])
subplot(2,2,2)
title(['Right Miss- ' num2str(indRM)])
ylim([-.05 .2])
subplot(2,2,3)
title(['Left Hit- ' num2str(indLH)])
ylim([-.05 .2])
subplot(2,2,4)
title(['Left Miss- ' num2str(indLM)])
ylim([-.05 .2])
suptitle([mouse ' ' date '- choiceAlign'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byOutcome.pdf']), '-dpdf')


tt = (-19:30)*frame_rate;
if input.doRandProb
    figure;
    tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
    nprob = length(input.ProbList);
    probs = cell2mat(input.ProbList);
    indR = zeros(nprob,ncon);
    indL = zeros(nprob,ncon);
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            indR(iprob,icon) = length(ind1); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            indL(iprob,icon) = length(ind2); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')

    figure;
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(SIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Hit- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- choiceAlign'])
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byProb.pdf']), '-dpdf', '-bestfit')

    figure;
    tProbLeft = celleqel2mat_padded(input.tStimProbAvgLeft);
    nprob = length(input.ProbList);
    probs = cell2mat(input.ProbList);
    indR = zeros(nprob,ncon);
    indL = zeros(nprob,ncon);
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            indR(iprob,icon) = length(ind1); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));           
            indL(iprob,icon) = length(ind2); 
            plot(tt, nanmean(nanmean(data_stim_dfof(:,good_ind,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- stimAlign'])
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimAlignResp_bySide_byProb_Miss.pdf']), '-dpdf', '-bestfit')

    figure;
    for iprob = 1:nprob
        for icon = 1:ncon
            subplot(2,nprob,iprob)
            ind1 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(~tLeftTrial), find(tGratingContrast == cons(icon)))));
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind1),2),3))
            hold on
            if icon == ncon
                title(['Right Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indR(iprob,:))])
                ylim([-.05 .2])
            end
            subplot(2,nprob,iprob+nprob)
            ind2 = intersect(find(tProbLeft == probs(iprob)), intersect(find(FIx), intersect(find(tLeftTrial), find(tGratingContrast == cons(icon)))));            
            plot(tt, nanmean(nanmean(data_choice_dfof(:,good_ind,ind2),2),3))
            hold on
            if icon == ncon
                title(['Left Miss- ' num2str((1-probs(iprob))*100) '% Right- ' num2str(indL(iprob,:))])
                ylim([-.05 .2])
            end
        end
    end
    suptitle([mouse ' ' date '- choiceAlign'])
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_choiceAlignResp_bySide_byProb_Miss.pdf']), '-dpdf', '-bestfit')

end

%% Size tuning quantification
m = @(pars,Stims) DiffofGauss_LG(pars,Stims);

%BaseA CA SA C0 S0 sigmaC sigmaS
s.lb =  [.001 .001 log2(min(Stims,[],2)) log2(0) log2(0) .001];
s.ub =  [1 1 log2(max(Stims,[],2)) log2(1000) log2(1000) 1];

clear temp;
index = 1;

names = {'x2','resnorm'};
args = cell(5);
x = log2(Stims);
Nshuf = 500;
resp_dFoverF = squeeze(mean(data_dfof(nOff:nOn+nOff,:,:),1));
base_dFoverF = squeeze(mean(data_dfof(nOff/2:nOff,:,:),1));
p_ttest = zeros(nCells,nStim);
h_ttest = zeros(nCells,nStim);
for count_shuf = 0:Nshuf
    Im_mat_USE = zeros(nCells, nStim);
    for iCond = 1:nStim        
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;
            [h_ttest(:,iCond) p_ttest(:,iCond)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
    end
    start = 1;
    for iCell = 1:nCells;
        if count_shuf == 0
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
    end
        
    figure;
    [n n2] = subplotn(nCells);
    for iCell = 1:nCells
    y = squeeze(tuning_mat(:,1,iCell));
    [max_amp max_ind] = max(y);
    if max_ind == length(x)
        ind2 = max_ind;
    else
        ind2 = max_ind+1;
    end

    s.x0 =  [2*max_amp (max_amp-y(end)) x(max_ind) 2 4 min(y)]; 
    options = optimset('Display', 'off');
    [x2,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB]  = lsqcurvefit(m,s.x0,x',y,s.lb,s.ub,options);
    temp = cell2struct({x2,Resnorm},names,2);
    s.k2_plot_oversamp0 = m(temp.x2,min(x,[],2):.1:max(x,[],2));

    subplot(n, n2, iCell)
    scatter(x',y);
    hold on; plot(min(x,[],2):.1:max(x,[],2), s.k2_plot_oversamp0);
    end
end

%% Contrast and Size tuning
if input.doConStim & input.doSizeStim & input.doMatrix
    
    if nCells<37
        [n, n2] = subplotn(nCells);
    else
        [n, n2] = subplotn(36);
    end

    con = cell2mat(input.tGratingContrast);
    cons = unique(con);
    nStim = length(cons).*length(szs);
    Stims = [];
    for isz = 1:length(szs)
        for icon = 1:length(cons)
            Stims = [Stims; szs(isz) cons(icon)];
        end
    end
    tuning_mat = zeros(nStim,nCells);
    figure;
    start = 1;
    f = 1;
    for iCell = 1:nCells
        if start >36
            print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeByConTuning' num2str(f) '.pdf']), '-dpdf')
            figure;
            f = 1+f;
            start = 1;
        end
        for iStim = 1:nStim
            ind1 = find(sz_mat == Stims(iStim,1));
            ind2 = find(con == Stims(iStim,2));
            ind = intersect(ind1,ind2);
            tuning_mat(iStim,1,iCell) = mean(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
            tuning_mat(iStim,2,iCell) = std(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
            Ind_struct(iStim).all_trials = ind;
        end
        size_con_mat = reshape(tuning_mat(:,1,iCell), [length(cons) length(szs)]);
        size_con_err = reshape(tuning_mat(:,2,iCell), [length(cons) length(szs)]);
        subplot(n,n2,start)
        for icon = 1:length(cons)
            errorbar(szs, size_con_mat(icon,:), size_con_err(icon,:), '-o', 'MarkerEdgeColor', repmat(cons(icon), [1 3]), 'Color', repmat(cons(icon), [1 3]))
            hold on
        end
        start = start +1; 
    end
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeByConTuning' num2str(f) '.pdf']), '-dpdf')
end

save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims')

%% retinotopy for these cells

ret_run_str = ['runs-' ret_ImgFolder(1,:)];
if ret_nrun>1
    ret_run_str = [ret_run_str '-' ret_ImgFolder(ret_nrun,:)];
end

data = [];
for irun = 1:ret_nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ret_ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ret_ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ret_ImgFolder(irun,:) '_000_000'],0,nframes);
    
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end
clear data_temp

%load average FOV from experiment
if exist(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str]))
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str]))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_reg_shifts.mat']), 'out', 'data_avg')
end
clear data

%load mWorks file
clear input temp
for irun = 1:ret_nrun
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' ret_time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
end
input = concatenateDataBlocks(temp);

nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2);

data_reg = data_reg(:,:,1:ntrials*(nOn+nOff));
sz = size(data_reg);
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr

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
data_dfof_avg_ret = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

%load masks from experiment
load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

figure; 
mask_all = mask_cell;
mask_all(find(mask_all>0)) = 1;
for i = 1:nStim; 
    subplot(length(Els),length(Azs),i);
    shade_img = imShade(data_dfof_avg_ret(:,:,i), mask_all);
    imagesc(shade_img)
    title(num2str(Stims(i,:)))
    clim([0 max(data_dfof_avg_all(:))]);
    colormap(gray)
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_FOVresp.pdf']), '-dpdf')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_FOVs.mat']), 'data_dfof_avg_ret', 'Azs', 'Els','Stims')
    
%extract timecourses and subtract neuropil
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(sz(3)./down, nCells);
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
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_reg data_reg_down 

%get dF/F
nCells = size(npSub_tc,2);
data_mat = zeros(nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
end
data_f = mean(data_mat(nOff/2:nOff,:,:),1);
data_df = bsxfun(@minus, data_mat, data_f);
data_dfof = bsxfun(@rdivide, data_df, data_f);
clear data_mat data_f data_df


tuning_mat = zeros(nStim, 2, nCells);
Ind_struct = [];
if nCells<37
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./frame_rate);
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind1 = find(Az == Stims(iStim,2));
        ind2 = find(El == Stims(iStim,1));
        ind = intersect(ind1,ind2);
        plot(tt', squeeze(mean(data_dfof(:,iCell,ind),3)))
        hold on
        tuning_mat(iStim,1,iCell) = mean(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
        tuning_mat(iStim,2,iCell) = std(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
        Ind_struct(iStim).all_trials = ind;
    end
    ylim([-0.05 0.25])
    vline(nOff)
    start = start + 1;
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    ret_mat = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)]);
    ret_mat = ret_mat';
    imagesc(ret_mat)
    colormap gray
    clim([0 max(max(tuning_mat(:,1,:),[],1),[],3)])
    title(num2str(chop(max(tuning_mat(:,1,iCell),[],1),2)))  
    start = start +1;
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims', 'Ind_struct')

%% fit retinotopy data 

% figure; 
% for iCell = 6
%     for iCond = 1:nStim
%         subplot(5,7,iCond)
%         ind_all = Ind_struct(iCond).all_trials;
%         plot(squeeze(data_dfof(:,iCell,ind_all)))
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
resp_dFoverF = squeeze(mean(data_dfof(nOff:nOn+nOff,:,:),1));
base_dFoverF = squeeze(mean(data_dfof(nOff/2:nOff,:,:),1));
p_ttest = zeros(nCells,nStim);
h_ttest = zeros(nCells,nStim);
h_all = zeros(1,nCells);
for count_shuf = 0:Nshuf
    fprintf('.')
    Im_mat_USE = zeros(nCells, nStim);
    for iCond = 1:nStim        
        ind_all = Ind_struct(iCond).all_trials;
        if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
            ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
        else
            ind_all_1 = ind_all;
            [h_ttest(:,iCond) p_ttest(:,iCond)] = ttest(resp_dFoverF(:,ind_all), base_dFoverF(:,ind_all), 'tail', 'right', 'dim', 2, 'alpha', 0.05./(nStim-1));
        end
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
    end

    start = 1;
    for iCell = 1:nCells;
        if count_shuf == 0
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
        if count_shuf>0
            if h_all(1,iCell) == 0
                continue
            end
        end
        a = Im_mat_USE(iCell,:);
        if max(a,[],2) > 0
            b = reshape(a',length(Azs),length(Els));
            data = b';
            if count_shuf == 0
                PLOTIT_FIT = 1;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG_Ret
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG_Ret
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end               
    end
    if count_shuf == 0
        fn_out = fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFfits' num2str(ifig) '.pdf']);   
        print(fn_out,'-dpdf')
    end
end


fn_out = fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Fit_struct.mat']);   
save(fn_out, 'Fit_struct')

resp_ind = find(h_all);

 if Nshuf>1;
    for iCell = 1:nCells
        if ~isempty(Fit_struct(iCell).True)                
            eval(['tmp = Fit_struct(iCell).True.s_.x;']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Elhicut_10];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.Azhicut_10];']);
            fit_true_vec(iCell,:) = tmp;
        end
    end
    
    fit_shuf_vec = NaN(nCells,10,Nshuf);
    for count_shuf = 1:Nshuf
        for iCell = 1:nCells
            if ~isempty(Fit_struct(iCell).Shuf)
                eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Elhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.Azhicut_10];']);
                %fit is: %A sigma_Az sigma_El Az0 El0 xi
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end

    Npars = size(fit_shuf_vec,2);
    lbub_fits = NaN(nCells,Npars,5);
    alpha_bound = .025;
    for iCell = 1:nCells
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp);
            ind_shuf_lb = ceil(Nshuf*alpha_bound);
            ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb);
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub);
            lbub_fits(iCell,count2,3) = mean(i); 
            lbub_fits(iCell,count2,5) = std(i);
        end
        %now take means from truedata fit:
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
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

fn_out = fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_lbub_fits.mat']);   
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind', 'resp_ind')


figure
subplot(2,2,1)
for i = goodfit_ind
    plot(lbub_fits(i,4,4), lbub_fits(i,5,4), 'o')
    hold on;
    xlim([min(Azs,[],2) max(Azs,[],2)])
    ylim([min(Els,[],2) max(Els,[],2)])
end
axis equal
title('RF center')
subplot(2,2,2)
for i = goodfit_ind
    ellipse(lbub_fits(i,2,4), lbub_fits(i,3,4), 0, lbub_fits(i,4,4), lbub_fits(i,5,4));
    hold on;
end
axis equal
title('RF- 1 sigma')
fn_out = fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFs.pdf']);   
print(fn_out,'-dpdf')

figure;
subplot(2,2,1)
hist(lbub_fits(goodfit_ind,2,4))
title('Sigma azimuth')
subplot(2,2,2)
hist(lbub_fits(goodfit_ind,3,4))
title('Sigma elevation')
subplot(2,2,3)
a = lbub_fits(goodfit_ind,3,4).*lbub_fits(goodfit_ind,2,4).*pi;
hist(a)
title('Area')
subplot(2,2,4)
scatter(a, lbub_fits(goodfit_ind,1,4),'o')
xlabel('Area')
ylabel('Peak dF/F')
fn_out = fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFdists.pdf']);   
print(fn_out,'-dpdf')

%% annulus for these cells

nrun = size(ann_ImgFolder,1);

ann_run_str = ['runs-' ann_ImgFolder(1,:)];
if nrun>1
    ann_run_str = [ann_run_str '-' ann_ImgFolder(nrun,:)];
end

data = [];
for irun = 1:nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ann_ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ann_ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ann_ImgFolder(irun,:) '_000_000'],0,nframes);
    
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
end
clear data_temp

%load average FOV from experiment
if exist(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str]))
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str]))
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_reg_shifts.mat']), 'out', 'data_avg')
end
clear data

%load mWorks file
clear input temp
for irun = 1:nrun
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' ann_time(irun,:) '.mat'];
    load(fName);
    temp(irun) = input;
end
input = concatenateDataBlocks(temp);

nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2);

data_reg = data_reg(:,:,1:ntrials*(nOn+nOff));
sz = size(data_reg);
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f); 
clear data_f data_df data_tr
    
Ann = double(celleqel2mat_padded(input.tAnnulusGratingDiameterDeg))-double(input.greymaskDiameterDeg);
Anns = unique(Ann);
nStim = length(Anns);
Stims = Anns;
data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, nStim);
for iStim = 1:nStim
    ind = find(Ann == Anns(iStim));
    data_dfof_avg(:,:,:,iStim) = mean(data_dfof(:,:,:,ind),4);
end
clear data_dfof
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));

figure; 
for i = 1:nStim; 
    subplot(2,2,i); 
    imagesc(data_dfof_avg_all(:,:,i));
    clim([0 .3])
    title(Anns(i));
    colormap(gray)
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_FOVresp.pdf']), '-dpdf')

%load masks from experiment
load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

%extract timecourses and subtract neuropil
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(sz(3)./down, nCells);
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
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
clear data_reg data_reg_down 

%get dF/F
nOn = input.nScansOn;
nOff = input.nScansOff;
nCells = size(npSub_tc,2);
data_mat = zeros(nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
end
data_f = mean(data_mat(nOff/2:nOff,:,:),1);
data_df = bsxfun(@minus, data_mat, data_f);
data_dfof = bsxfun(@rdivide, data_df, data_f);
clear data_mat data_f data_df

tuning_mat = zeros(nStim, 2, nCells);
if nCells<37
    [n, n2] = subplotn(nCells);
else
    [n, n2] = subplotn(36);
end
tt= (1-nOff:nOn)*(1000./frame_rate);
figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    for iStim = 1:nStim
        ind = find(Ann == Anns(iStim));
        plot(tt', squeeze(mean(data_dfof(:,iCell,ind),3)))
        hold on
        tuning_mat(iStim,1,iCell) = mean(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),3),1);
        tuning_mat(iStim,2,iCell) = std(mean(data_dfof(nOff+1:nOn+nOff,iCell,ind),1),[],3)./sqrt(length(ind));
    end
    ylim([-0.05 0.25])
    vline(nOff)
    start = start + 1;
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    errorbar(Anns, tuning_mat(:,1,iCell),tuning_mat(:,2,iCell), '-o');
    start = start +1;
end
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims')

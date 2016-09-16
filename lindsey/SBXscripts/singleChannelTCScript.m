%% get path names
date = '160714';
ImgFolder = strvcat('002', '003', '004');
time = strvcat('1712', '1726', '1742');
mouse = 'i924';
nrun = size(ImgFolder,1);
frame_rate = 15;

%% load and register

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

data = [];
clear temp
trial_n = [];
for irun = 1:nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
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
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    data_avg = mean(data(:,:,1:100),3);
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
end
clear data

%% find activated cells
    
%max by trial type    
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    ntrials = size(input.tGratingDirectionDeg,2);

    sz = size(data_reg);
    data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(data_tr), data_f); 
    data_dfof = bsxfun(@rdivide,data_df, data_f); 
    clear data_f data_df data_tr
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
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')

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
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_resp.pdf']), '-dpdf')
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
    data_dfof_avg = zeros(sz(1),sz(2),nOn+nOff,length(Dirs));
    for idir = 1:length(Dirs)
        ind = find(Dir == Dirs(idir));
        data_dfof_avg(:,:,:,idir) = mean(data_dfof(:,:,:,ind),4);
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = squeeze(std(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),[],3));
    data_dfof_avg_avg = squeeze(mean(mean(imfilter(data_dfof_avg(:,:,nOff/2:nOff,:),myfilter),3),4));

    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i));
        clim([0 max(data_dfof_avg_all(:))])
        title(num2str(nDirs))
        colormap(gray)
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3);
        end
    end
end

%% segmentation 

mask_all = zeros(sz(1), sz(2));
%mask_data = squeeze(max(reshape(data_dfof_avg_all, [sz(1) sz(2) 2 nStim/2]),[],3));
mask_data = data_dfof_avg_all;
% figure;
% [n, n2] = subplotn(size(mask_data,3));
% for iStim = 1:nStim
%     subplot(n, n2, iStim)
%     imagesc(mask_data(:,:,iStim))
%     colormap gray
% end

for iStim = 1:nStim    
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_all > 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_2 = bwlabel(bwout);
    mask_all = mask_all+mask_2;
    close all
end


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
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')

mask_cell = bwlabel(mask_all);
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_avg_all', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2
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
clear data_reg data_reg_down

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%% extract tuning
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
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning'  num2str(f) '.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TuningSummary.mat']), 'tuning_mat', 'tc_mat', 'Ind_struct')

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
            print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeByConTuning' num2str(f) '.pdf']), '-dpdf')
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
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_SizeByConTuning' num2str(f) '.pdf']), '-dpdf')
end

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims')

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
if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str]))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str]))
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_reg_shifts.mat']), 'out', 'data_avg')
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
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_FOVresp.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_FOVs.mat']), 'data_dfof_avg_ret', 'Azs', 'Els','Stims')
    
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
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
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
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims', 'Ind_struct')

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
        fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFfits' num2str(ifig) '.pdf']);   
        print(fn_out,'-dpdf')
    end
end


fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_Fit_struct.mat']);   
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

fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_lbub_fits.mat']);   
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
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFs.pdf']);   
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
fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_run_str], [date '_' mouse '_' ret_run_str '_RFdists.pdf']);   
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
if exist(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str]))
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_reg_shifts.mat']))
    [outs, data_reg]=stackRegister_MA(data,[],[],out);
    clear out outs
else
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str]))
    save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_reg_shifts.mat']), 'out', 'data_avg')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_FOVresp.pdf']), '-dpdf')

%load masks from experiment
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

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
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
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
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')
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
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_TCs' num2str(f) '.pdf']), '-dpdf')

figure;
start = 1;
f = 1;
for iCell = 1:nCells
    if start >36
        print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
        start = 1;
        f= f+1;
        figure;
    end
    subplot(n, n2, start)
    errorbar(Anns, tuning_mat(:,1,iCell),tuning_mat(:,2,iCell), '-o');
    start = start +1;
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning' num2str(f) '.pdf']), '-dpdf')
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ann_run_str], [date '_' mouse '_' ann_run_str '_Tuning.mat']), 'data_dfof', 'tuning_mat', 'Stims')

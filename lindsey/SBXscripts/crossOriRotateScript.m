clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRotate_ExptList';
iexp = 6; 
rc = behavConstsAV;
eval(ds)

frame_rate = 15;

%%
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
%LG_base = '\\CRASH.dhe.duke.edu\data\home\lindsey';

fprintf(['2P imaging rotation analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load and register
tic
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\ashley\data\AW68\two-photon imaging\' date '\' ImgFolder(irun,:)];
    %CD = [LG_base '\Data\2P_images\' mouse '-KM' '\' date '_' mouse '\' ImgFolder(irun,:)];
    %CD = ['\\CRASH.dhe.duke.edu\data\home\kevin\Data\2P\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder{irun} '_000_000.mat'];
    load(imgMatFile);
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    
    temp(irun) = input;
    nframes = [temp(irun).counterValues{end}(end) info.config.frames];
    
    fprintf(['Reading run ' num2str(irun) '- ' num2str(min(nframes)) ' frames \r\n'])
    data_temp = sbxread(imgMatFile(1,1:11),0,min(nframes));
    if size(data_temp,1)== 2
        data_temp = data_temp(1,:,:,:);
    end
    
    if isfield(input, 'cStimOneOn') 
        if irun>1
            ntrials = size(input.cStimOneOn,2);
            for itrial = 1:ntrials
                temp(irun).cStimOneOn{itrial} = temp(irun).cStimOneOn{itrial}+offset;
                temp(irun).cStimOneOff{itrial} = temp(irun).cStimOneOff{itrial}+offset;
            end
        end
    end
    offset = offset+min(nframes);
        
    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateStructuresLG(temp);
clear data_temp
clear temp
toc

%% Choose register interval
nep = floor(size(data,3)./10000);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); colormap gray; clim([0 3000]); end
movegui('center')
data_avg = mean(data(:,:,30001:30500),3);
%% Register data
if exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    [outs, data_reg] = stackRegister_MA(data,[],[],out);
else
    [out, data_reg] = stackRegister(data,data_avg);
    mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_avg')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data out

%% test stability
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3)); title([num2str(1+((i-1)*10000)) '-' num2str(500+((i-1)*10000))]); end
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_byFrame.pdf']),'-dpdf', '-bestfit')
movegui('center')
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf', '-bestfit')
movegui('center')
i = 1;
sz = size(data_reg);
rg = zeros(sz(1),sz(2),3);
first = mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3);
rg(:,:,1) = first./max(first(:));
i = nep; 
last = mean(data_reg(:,:,1+((i-1)*10000):500+((i-1)*10000)),3);
rg(:,:,2) = last./max(last(:));
figure; image(rg)
movegui('center')
%% find active cells
sz = size(data_reg);
cStimOn = celleqel2mat_padded(input.cStimOneOn);
down = 5;

nTrials = length(cStimOn);
ind_nomask = find(celleqel2mat_padded(input.tMaskOneGratingContrast)==0);
ind_mask = find(celleqel2mat_padded(input.tMaskOneGratingContrast)>0);
nOn = unique(celleqel2mat_padded(input.nStimOneFramesOn));
prewin_frames = 30;
data_trial_ind = zeros(sz(1),sz(2),(nOn+prewin_frames)/down,length(ind_nomask));
for i = 1:length(ind_nomask)
    data_trial_ind(:,:,:,i) = stackGroupProject(data_reg(:,:,cStimOn(ind_nomask(i))-prewin_frames:cStimOn(ind_nomask(i))+nOn-1),down);
end
data_trial_f = mean(data_trial_ind(:,:,1:prewin_frames/down,:),3);
data_trial_df = bsxfun(@minus,data_trial_ind,data_trial_f);
data_trial_dfof = bsxfun(@rdivide, data_trial_df,data_trial_f);
data_nomask_max = max(mean(data_trial_dfof,4),[],3);
figure; movegui('center');
imagesc(data_nomask_max)
clear data_trial_ind data_trial_f data_trial_dfof

maskPhas_all = celleqel2mat_padded(input.tMaskOneGratingPhaseDeg);
maskPhas = unique(maskPhas_all);
nMaskPhas = length(maskPhas);
if isfield(input,'tDoMaskRotate')
    tMaskRot_all = celleqel2mat_padded(input.tDoMaskRotate);
else
    tMaskRot_all = zeros(size(maskPhas_all));
end
tMaskRot = unique(tMaskRot_all);
nR = length(tMaskRot);
data_mask_max = zeros(sz(1),sz(2),nMaskPhas,nR);
for ir = 1:nR
    ind_r = intersect(ind_mask,find(tMaskRot_all == tMaskRot(ir)));
    figure; movegui('center');
    for ip = 1:nMaskPhas
        ind_p = intersect(ind_r, find(maskPhas_all == maskPhas(ip)));
        data_trial_ind = zeros(sz(1),sz(2),(nOn+prewin_frames)/down,length(ind_p));
        for i = 1:length(ind_p)
            data_trial_ind(:,:,:,i) = stackGroupProject(data_reg(:,:,cStimOn(ind_p(i))-prewin_frames:cStimOn(ind_p(i))+nOn),down);
        end
        data_trial_f = mean(data_trial_ind(:,:,1:prewin_frames/down,:),3);
        data_trial_df = bsxfun(@minus,data_trial_ind,data_trial_f);
        data_trial_dfof = bsxfun(@rdivide, data_trial_df,data_trial_f);
        data_mask_max(:,:,ip,ir) = max(mean(data_trial_dfof,4),[],3);
        subplot(2,2,ip)
        imagesc(data_mask_max(:,:,ip,ir))
    end
end
clear data_trial_ind data_trial_f data_trial_dfof data_trial_df
data_mask_max = reshape(data_mask_max,[sz(1) sz(2) nMaskPhas.*nR]);
data_dfof = cat(3,data_mask_max, data_nomask_max);

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']), 'cStimOn', 'maskPhas_all', 'maskPhas', 'nMaskPhas', 'tMaskRot_all', 'tMaskRot', 'nR', 'frame_rate', 'nOn', 'ind_mask', 'ind_nomask', 'nTrials')

%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
    mask_data_temp = mask_data(:,:,end+1-iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractiveLG(mask_data_temp);
    mask_all = mask_all+bwout;
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)

%% neuropil subtraction
mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof', 'mask_cell', 'mask_np')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_data_temp mask_exp data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 


down = 5;
sz = size(data_reg);

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
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

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc', 'nCells', 'sz')

clear data_tc data_tc_down np_tc np_tc_down mask_np mask_cell

%%
nCells = size(npSub_tc,2);
data_trial_nomask = zeros(nOn+(2*frame_rate),nCells,length(ind_nomask));
for i = 1:length(ind_nomask)
    if size(npSub_tc,1)>cStimOn(ind_nomask(i))+nOn+frame_rate-1
        data_trial_nomask(:,:,i) = npSub_tc(cStimOn(ind_nomask(i))-frame_rate:cStimOn(ind_nomask(i))+nOn+frame_rate-1,:);
    else
        data_trial_nomask(:,:,i) = [];
    end
end
data_trial_nomask_f = mean(data_trial_nomask(1:frame_rate,:,:),1);
data_trial_nomask_df = bsxfun(@minus, data_trial_nomask, data_trial_nomask_f);
data_trial_nomask_dfof = bsxfun(@rdivide, data_trial_nomask_df, data_trial_nomask_f);
data_nomask_dfof = smoothdata(mean(data_trial_nomask_dfof,3),1);
[max_val max_ind] = max(data_nomask_dfof(frame_rate+1:end,:),[],1);
[val ind] = sort(max_ind,'ascend');
figure;
movegui('center')
imagesc(data_nomask_dfof(frame_rate+1:end,ind)')
ylabel('Cells')
xlabel('Time')
title('No mask')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sortTC_nomask.pdf']),'-dpdf', '-bestfit')


data_mask_dfof = zeros(nOn+(2*frame_rate),nCells,nMaskPhas,nR);
for ir = 1:nR
    ind_r = intersect(ind_mask,find(tMaskRot_all == tMaskRot(ir)));
    for ip = 1:nMaskPhas
        ind_p = intersect(ind_r, find(maskPhas_all == maskPhas(ip)));
        data_trial_mask = zeros(nOn+(2*frame_rate),nCells,length(ind_p));
        for i = 1:length(ind_p)
            if size(npSub_tc,1)>cStimOn(ind_p(i))+nOn+frame_rate-1
                data_trial_mask(:,:,i) = npSub_tc(cStimOn(ind_p(i))-frame_rate:cStimOn(ind_p(i))+nOn+frame_rate-1,:);
            else
                data_trial_mask(:,:,i) = [];
            end
        end
        data_trial_mask_f = mean(data_trial_mask(1:frame_rate,:,:),1);
        data_trial_mask_df = bsxfun(@minus, data_trial_mask, data_trial_mask_f);
        data_trial_mask_dfof = bsxfun(@rdivide, data_trial_mask_df, data_trial_mask_f);
        data_mask_dfof(:,:,ip,ir) = smoothdata(mean(data_trial_mask_dfof,3),1);
    end
end
mask_str = {'stationary','rotate'};
for ir = 1:nR
figure;
movegui('center')
for ip = 1:nMaskPhas
    subplot(2,2,ip)
    imagesc(data_mask_dfof(frame_rate+1:end,ind,ip,ir)')
    ylabel('Cells')
    xlabel('Time')
    title(num2str(maskPhas(ip)))
end
suptitle(['Mask- ' mask_str(ir)])
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sortTC_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')
end

save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respAvg.mat']),'data_nomask_dfof', 'data_mask_dfof')
%%
yB = cell(nR,nMaskPhas);
yA = cell(nR,1);
for ir = 1:nR
    X = data_nomask_dfof(frame_rate+1:end,:);
    nframes = size(X,1);
    for ip = 1:nMaskPhas
        X = [X; data_mask_dfof(frame_rate+1:end,:,ip,ir)];
    end

    figure;
    [Y,e]=cmdscale(pdist(X,'cos'));
    plot(e(1:20),'k.-','markersize',25)
    box off
    xlabel('Rank')
    ylabel(['Eigenvalue'])
    title(['Mask- ' mask_str(ir)])
    movegui('center')
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_eigenvals_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')


    figure;
    u = 1:nframes;   % indices for masked and unmasked conditions
    yA{ir}= Y(u,:);
    m = u;
    for i = 1:nMaskPhas
        subplot(2,2,i)
        plot3(yA{ir}(:,1),yA{ir}(:,2),yA{ir}(:,3),'b-','LineWidth',2)
        hold on
        m = m+nframes;
        yB{ir,i} = Y(m,:);
        plot3(yB{ir,i}(:,1),yB{ir,i}(:,2),yB{ir,i}(:,3),'r-','LineWidth',2)
        title(num2str(maskPhas(i)))
        view([7.244 11.396])
        hold off
        xlim([-0.5 0.5])
        ylim([-0.5 0.5])
        zlim([-0.5 0.5])
    end
    movegui('center')
    suptitle(['Mask- ' mask_str(ir)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allTrajectories_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')

    figure;
    plot3(yA{ir}(:,1),yA{ir}(:,2),yA{ir}(:,3),'LineWidth',2)
    legstr = cell(1,nMaskPhas+1);
    legstr{1} = 'Grating';
    hold on
    for i = 1:nMaskPhas
        plot3(yB{ir,i}(:,1),yB{ir,i}(:,2),yB{ir,i}(:,3),'LineWidth',2)
        hold on
        legstr{i+1} = num2str(maskPhas(i));
    end
    view([7.244 11.396])
    hold off
    legend(legstr)
    movegui('center')
    title(['Mask- ' mask_str(ir)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allTrajectoriesOverlay_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')

    figure;
    plot(yA{ir}(:,1),yA{ir}(:,2),'LineWidth',2)
    legstr = cell(1,nMaskPhas+1);
    legstr{1} = 'Grating';
    hold on
    for i = 1:nMaskPhas
        plot(yB{ir,i}(:,1),yB{ir,i}(:,2),'LineWidth',2)
        hold on
        legstr{i+1} = num2str(maskPhas(i));
    end
    legend(legstr)
    movegui('center')
    title(['Mask- ' mask_str(ir)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allTrajectoriesOverlay2D_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')

    ndim = 6;            % select embedding dimension
    figure
    for i = 1:nMaskPhas
        yAc = yA{ir}(:,1:ndim);
        yBc = yB{ir,i}(:,1:ndim);

        yAn = bsxfun(@rdivide,yAc,sqrt(sum(yAc.^2,2))); % normalize responses
        yBn = bsxfun(@rdivide,yBc,sqrt(sum(yBc.^2,2)));

        yAh = [yAn ones(size(yAc,1),1)]; % put in homogeneous coordinates 
        yBh = [yBn ones(size(yBc,1),1)];

        A = pinv(yAh)*yBh;              % estimate affine transform
        yBhp = yAh*A;                   % estimate prediction
        subplot(2,2,i)
        plot3(yAh(:,1),yAh(:,2),yAh(:,3),'b-','LineWidth',2)
        hold on
        plot3(yBh(:,1),yBh(:,2),yBh(:,3),'r-','LineWidth',2)
        plot3(yBhp(:,1),yBhp(:,2),yBhp(:,3),'g-','LineWidth',2)
        axis equal
        title(num2str(maskPhas(i)))
        view([-149.29 12.17])
    end
    movegui('center')
    suptitle(['Mask- ' mask_str(ir)])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_allTrajectoriesAffine_maskRotate' num2str(ir-1) '.pdf']),'-dpdf', '-bestfit')
end
save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cosSpace.mat']),'yA', 'yB')

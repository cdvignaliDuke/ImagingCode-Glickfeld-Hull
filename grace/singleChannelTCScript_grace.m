%% clear everything
clear all
clear global
%% get path names
date = '200108';
ImgFolder = strvcat('003');
run = strvcat('000');
time = strvcat('1158');
mouse = 'i1316';
doFromRef = 1;
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
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
figure; imagesq(mean(data_reg(:,:,1:10000),3)); truesize;
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']), '-dpdf')

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

if input.doDirStim
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
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')
   
    figure;
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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_16Stim.pdf']), '-dpdf')

    
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
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')
end

%% Register day x to day 1
%load reference data
zoom = '2.0';
mod = '200um';
if doFromRef
    ref_str = ['runs-' ref_run];
    if size(ref_run,1)>1
        ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
    end
    %reload and reregister reference day data
%     refData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_reg_shifts.mat']));
%     data_ref_avg = refData.data_reg_avg;
    refGFP = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_FOV_Check'], [ref_date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.mat']));
    data_ref_avg = refGFP.data_reg_avg;
    regGFP = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_FOV_Check'], [date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.mat']));
    data_reg_avg = regGFP.data_reg_avg;
    
    %create averages of current day and reference day
    ref = data_ref_avg;
    ref(find(ref>1e5)) = 0;
    ref = ref./max(max(abs(ref)));
    
    reg = data_reg_avg;
    reg(find(reg>1e5)) = 0;
    reg = reg./max(max(abs(reg)));
    
    [i_pts b_pts] = cpselect(double(reg), double(ref),'Wait',true);
    sz_target = size(reg);
    mytform = maketform('affine',i_pts(1:3,:),b_pts(1:3,:));
    img = imtransform(double(data_reg_avg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    dfof_reg2ref = imtransform(double(data_dfof_max),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    data_reg_to_ref = imtransform(double(reg),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RegData.mat']), 'dfof_reg2ref', 'reg', 'ref', 'data_reg_to_ref')

    %register ref avg to reg avg 
%     [out_avg, img] = stackRegister(data_reg_avg, data_ref_avg);
%     [out, data_reg_to_ref] = stackRegister_MA(data_reg,[],[],repmat(out_avg,[size(data_reg,3) 1]));
	figure; 
    subplot(2,1,1)
    imagesc(ref); truesize;
    title('Reference image')
    subplot(2,1,2)
    imagesc(reg); truesize;
    title('Todays image')
    figure;
    imagesc(img); truesize;
    title('Registered image')
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_registered_FOV.pdf']), '-dpdf')
end

%% cell segmentation 

if ~doFromRef
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
        print(fullfile('Z:\All_staff\home\grace\Analysis\2P',[date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')
        
    mask_np = imCellNeuropil(mask_cell, 3, 5);
    save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')
    
    clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
end

%% neuropil subtraction

if doFromRef
    ref_str = ['runs-' ref_run];
    if size(ref_run,1)>1
        ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
    end
refMaskData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
mask_cell =  refMaskData.mask_cell;
    for i = 1:nframes
        data_reg(:,:,i) = imtransform(double(data_reg(:,:,i)),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
        if rem(i,50) == 0
            fprintf([num2str(i) '/n'])
        end
    end
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


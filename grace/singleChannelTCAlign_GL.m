clear all
clear global
%% get path names
date = '200120';
ImgFolder = strvcat('003');
time = strvcat('1214');
mouse = 'i1312';
alignToRef = 1;
ref_date = '200118';
ref_run = strvcat('002');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(ref_run, size(ref_run,1));
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\\2P_images';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '\' ImgFolder(irun,:)]);
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
    
    offset = offset+nframes;

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

data_avg = mean(data(:,:,8001:8500),3);

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
    data_reg_avg = mean(data_reg(:,:,1:10000),3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data

%% test stability/pixel correlation
% figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end
figure; imagesq(data_reg_avg); truesize;
% writetiff(data_reg_avg, ['Z:\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_FOV_Check\' date '_' mouse '_run' run '_zoom' zoom '_' mod '_avgFOVgreen.tiff'])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')

load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']));
data_reg_3hz = stackGroupProject(data_reg,5);
% nFrames = size(data_reg_3hz,3);
% reg = data_reg_avg;
% sz_target = size(reg);
%     for i = 1:nFrames
%         data_reg_3hz(:,:,i) = imtransform(double(data_reg_3hz(:,:,i)),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);
%         if rem(i,50) == 0
%             fprintf([num2str(i) '/n'])
%         end
%     end
pix = getPixelCorrelationImage(data_reg_3hz);
pix(isnan(pix))=0;

% transforms
fitGeoTAf = fitgeotrans(input_points(:,:), base_points(:,:),'affine'); 
r2rFGTA = imwarp(double(reg),fitGeoTAf, 'OutputView', imref2d(size(reg)));
r2rFGTA_dfof = imwarp(double(data_dfof_max),fitGeoTAf, 'OutputView', imref2d(size(reg)));

fitGeoTNr = fitgeotrans(input_points(:,:), base_points(:,:),'nonreflectivesimilarity');
r2rFGTN = imwarp(double(reg),fitGeoTNr, 'OutputView', imref2d(size(reg)));
r2rFGTN_dfof = imwarp(double(data_dfof_max),fitGeoTNr, 'OutputView', imref2d(size(reg))); 
   
[optimizer, metric] = imregconfig('Monomodal');
imRegisTForm = imregtform(double(reg), double(ref), 'similarity', optimizer, metric);
r2rIMRT = imwarp(double(reg),imRegisTForm, 'OutputView', imref2d(size(reg))); 
% r2rIMRT_dfof = imwarp(double(data_dfof_max),imRegisTForm, 'OutputView', imref2d(size(reg))); 

% transforming pixel image
pix_fgta = imwarp(double(pix),fitGeoTAf, 'OutputView', imref2d(size(reg))); 
pix_fgtn = imwarp(double(pix),fitGeoTNr, 'OutputView', imref2d(size(reg))); 
pix_imrt = imwarp(double(pix),imRegisTForm, 'OutputView', imref2d(size(reg))); 

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pixel.mat']),'pix_fgta','pix')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']), 'input_points', 'base_points', 'data_dfof_max_ref', 'ref','reg','fitGeoTAf','r2rFGTA','r2rFGTA_dfof');
load(fullfile(fnout, [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']))
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_cell', 'mask_np')

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
%     print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_16Stim.pdf']), '-dpdf')

    
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
%     print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')

    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
    print(fullfile('Z:\All_staff\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofMax.pdf']), '-dpdf')
end

 load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all', 'nStim')

%% histogram data_reg_avg
figure;histogram(data_reg_avg(:),10,'BinLimits',[0,2000])
figure;histogram(data_dfof_max(:),10,'BinLimits',[0,2000])
%% cell segmentation 
if ~alignToRef
    mask_exp = zeros(sz(1),sz(2));
    mask_all = zeros(sz(1),sz(2));
    mask_data = data_dfof;

for iStim = 1:size(data_dfof,3)
  mask_data_temp = mask_data(:,:,iStim);
  mask_data_temp(find(mask_exp >= 1)) = 0;
  bwout = imCellEditInteractive(mask_data_temp);
  mask_all = mask_all+bwout;
  mask_exp = imCellBuffer(mask_all,3)+mask_all;
  close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell)


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
    
%     data_reg_avg_thresh = data_reg_avg;
%     data_reg_avg_thresh(find(data_reg_avg>600)) = 0;
    
    reg = data_avg_reg;
    ref = data_avg_ref;
%     reg(find(reg>6000)) = 0;
    reg = (reg./max(max(abs(reg))));
%     reg = imadjust(reg);
%     ref(find(ref>6000)) = 0;
    ref = (ref./max(max(abs(ref)))); 
%     ref = imadjust(ref);
    
%     figure; 
%     subplot(2,1,1);imagesc(reg); axis off; axis equal; title('data avg FOV day 2')
%     subplot(2,1,2);imagesc(ref); axis off; axis equal; title('data avg FOV day 1')
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_data_avg_FOVs.pdf']), '-dpdf','-bestfit')

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
    figure; imagesc(rgb_reg2ref); title(['day 2 on day 1, data reg rgb overlay'])
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays']))
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_reg2ref_overlay.pdf']), '-dpdf','-bestfit')
    figure; imagesc(rgb_reg2ref_dfof); title(['day 2 on day 1, dfof rgb overlay'])
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_registered_dfof_overlay.pdf']), '-dpdf','-bestfit')
    
    figure; 
    subplot(2,3,1); imagesc(ref); axis off; axis equal; title('day 1 data avg')
    subplot(2,3,2); imagesc(reg); axis off; axis equal; title('day 2 data avg')
    subplot(2,3,3); imagesc(reg2ref); axis off; axis equal; title('registered day 2 data avg')
    subplot(2,3,4); imagesc(rgb_reg2ref_dfof(:,:,1)); axis off; axis equal; title('day 1 data dfof')
    subplot(2,3,5); imagesc(data_dfof_max_reg); axis off; axis equal; title('day 2 data dfof')
    subplot(2,3,6); imagesc(rgb_reg2ref_dfof(:,:,2)); axis off; axis equal; title('registered day 2 data dfof')
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_rgb_individuals.pdf']), '-dpdf','-bestfit')

    figure;
    subplot(2,1,1); imagesc(data_dfof_max_reg); axis off; axis equal; title('day 2 before registration')
    subplot(2,1,2); imagesc(reg2ref_dfof); axis off; axis equal; title('day 2 after registration')
    print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], ['RGB Overlays'], [date '_' mouse '_' run_str '_before_and_after_reg.pdf']), '-dpdf','-bestfit')


    
    for i = 1:nframes
        data_reg(:,:,i) = imwarp(double(data_reg(:,:,i)),fitGeoTAf, 'OutputView', imref2d(size(reg)));
        if rem(i,50) == 0
            fprintf([num2str(i) '/n'])
        end
    end
    
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_transform.mat']), 'input_points', 'base_points', 'mytform', 'data_dfof_max_ref', 'ref', 'reg2ref', 'reg2ref_dfof','fitGeoTAf','r2rFGTA','r2rFGTA_dfof','fitGeoTNr','r2rFGTN','r2rFGTN_dfof','imRegisTForm','r2rIMRT','r2rIMRT_dfof');
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
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun  
    CD = fullfile(gl_fn, [mouse '\' date '\' ImgFolder(irun,:)]);
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
end
input = concatenateDataBlocks(temp);

 Dir = cell2mat(input.tGratingDirectionDeg);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
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
%     save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','max_dir','h_dir', 'h_ori', 'max_ori','good_ind')
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
    getOriTuningLG(data_tc_down(:,cells_all),input,tuningDownSampFactor);
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
for icell = 1:length(cells_all)
    ic = cells_all(icell);
    if start > 49
        suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,icell-(x.*49))
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
suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(length(cells_all)) '- well-fit'])
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

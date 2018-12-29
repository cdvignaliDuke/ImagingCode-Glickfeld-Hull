clear all
CRP_expt_list_Crus
id = 3;
jake_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Analysis\Cue_reward_pairing_analysis\2P';
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
nexp = size(expt(id).date,1);

for iexp = 1:nexp
    mouse = expt(id).mouse(iexp,:);
    date = expt(id).date(iexp,:);
    run = expt(id).run(iexp,:);
    time = expt(id).time(iexp,:);
    fprintf([date ' ' mouse '\n'])
    img_fn = [date '_' mouse];
    if ~exist(fullfile(lg_out,img_fn))
        mkdir(fullfile(lg_out, img_fn))
    end
    %% load and register
    if str2num(mouse)>800
        mouse_name = mouse(2:3);
    else
        mouse_name = mouse;
    end
    
    cd(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Data\2P_imaging', [date '_img' mouse_name],['img' mouse_name]))
    load(['img' mouse_name '_000_' run '.mat'])
    if exist(fullfile(lg_out, img_fn, [img_fn '_input.mat']))
        load(fullfile(lg_out, img_fn, [img_fn '_input.mat']))
    else
        load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-i' mouse '-' date '-' time '.mat'])
        mworks = input;
        clear input;
    end
    clear ttl_log
    nf= mworks.counterValues{end}(end);
    if mod(nf,10000)==2
        nf = nf-2;
    elseif mod(nf,10000)==1
        nf = nf-1;
    end
    
    if expt(id).ttl(iexp)
        load(['img' mouse_name '_000_' run '_realtime.mat'])
        fprintf('Loading data \n')
        data = sbxread(['img' mouse_name '_000_' run],0,nf);
        data = squeeze(data);
        ttl_trans = find(abs(diff(ttl_log)));
        n = length(ttl_trans);
        for i = 1:n
            ttl_log(ttl_trans(i)-4:ttl_trans(i)+4,:) = 0;
        end
        ttl_ind = find(ttl_log(1:nf,:) ==0);
        data(:,:,ttl_ind) = [];
    else
        if nf > 60000
            nf = 60000;
        end
        fprintf('Loading data \n')
        data = sbxread(['img' mouse_name '_000_' run],0,nf);
        data = squeeze(data);
    end
    if exist(fullfile(lg_out, img_fn, [img_fn '_reg.mat']))
        load(fullfile(lg_out, img_fn, [img_fn '_reg.mat']))
    else
       figure;
       nplot = floor(size(data,3)./5000);
       [n n2] = subplotn(nplot);
       for i = 1:nplot
           subplot(n,n2,i)
           imagesc(mean(data(:,:,1+(i-1).*5000:500+(i-1).*5000),3))
           title(num2str(i));
       end
       prompt = 'Choose a registration frame \n';
       x = input(prompt);
       img_ref = mean(data(:,:,1+(x-1).*5000:500+(x-1).*5000),3);
       save(fullfile(lg_out, img_fn, [img_fn '_reg.mat']), 'img_ref');
    end
    %%    
    if size(data,3) >= 60000
        data = data(:,:,1:2:end);
    end
    fprintf(['Registering ' num2str(size(data,3)) ' frames'])
    [outs img_reg]= stackRegister(data, img_ref);
    
    clear data
    [npw, nph, nt] = size(img_reg);
%     nframes(1,iexp) = nt;
%     nx = floor(nt/4);
%     figure;
%     for i = 1:4
%         subplot(2,2,i)
%         imagesc(mean(img_reg(:,:,(i-1)*nx + 1: (i-1)*nx + 500),3))
%         title([num2str((i-1)*nx + 1) ':' num2str((i-1)*nx + 500)]);
%     end
%     suptitle([date ' ' mouse ' stability'])
%     print(fullfile(lg_out, img_fn, [img_fn '_stability.pdf']),'-dpdf','-fillpage')
    %% PCA
    nPCA = 200; %100 for old datasets, 500 for newer
    nf = size(img_reg,3);
    [mixedsig, mixedfilters, CovEvals, ~, ~, ...
        ~] = CellsortPCA_2P(img_reg,[1 nf], nPCA,[], fullfile(lg_out,img_fn), img_fn, []);
    PCuse = 1:nPCA;
    %% ICA
    mu = 0.98; % spatial temporal mix
    nIC = 150;  %number of IC, 400- img90 100- img91
    termtol = 0.00001; % termination tolerance
    maxrounds = 1000; % #of iteration 

    %run ICA and save outputs
    [ica_sig, mixedfilters, ica_A, numiter] = CellsortICA_2P(mixedsig, ...
        mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);
    icasig = permute(mixedfilters,[2,3,1]);
    save(fullfile(lg_out, img_fn, [img_fn '_ICA.mat']), 'ica_sig', 'icasig');

    nIC = size(icasig,3);
    icasig = stackFilter(icasig); % filter IC with gaussian filter
    
    %% select ICs
    if expt(id).name == 'Simp'
        if id == 1 & iexp == 4
            cluster_threshold = 95;
        elseif id == 2 & iexp == 4
            cluster_threshold = 96;
        elseif id == 2 & iexp == 7
            cluster_threshold = 96.5;
        elseif id == 3 & iexp == 1
            cluster_threshold = 96.5;
        elseif id == 3 & iexp == 4
            cluster_threshold = 96;
        elseif id == 4 & iexp == 4
            cluster_threshold = 95.5;
        else
            cluster_threshold = 97;%percentile for thresholding 97- img90 97- img91   94.5- img044 97.5- img060 96- img064 97- img063 97- img065 XXXXXXXXXXXX
        end
    elseif expt(id).name == 'Crus'
        cluster_threshold = 97;
    end
    threshold = 0.8; % correlation threshold
    mask_cell = zeros(size(icasig));
    sm_logical = zeros(npw,nph);

    for ic = 1:nIC
        icasig(:,:,ic) = imclearborder(icasig(:,:,ic)); % remove ICs touching borders
        % thresholding to get each IC
        sm_logical((icasig(:,:,ic)>mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
        sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
        sm_logical = logical(sm_logical);
        mask_cell(:,:,ic) = bwlabel(sm_logical); % label multiple IC on single frame
    end
    
    % separate multiple ICs and remove small ones and deal with overlapping 
    mask_final = processMask(mask_cell);
    mask_raw = reshape(mask_final, npw, nph);
    figure; imagesc(mask_raw); truesize;
    
    % combine highly correlated ICs to one 
    [ ~, mask3D, ~] = finalMask(img_reg(:,:,1:5:end), mask_final, threshold, fullfile(lg_out,img_fn));
    
    pause
    
    %% Plotting TCs
    %reload data
    nf= mworks.counterValues{end}(end); 
    if mod(nf,10000)==2
        nf = nf-2;
    elseif mod(nf,10000)==1
        nf = nf-1;
    end
    fprintf('Loading data \n')
    data = sbxread(['img' mouse_name '_000_' run],0,nf);
    data = squeeze(data);
    [outs img_reg]= stackRegister(data, img_ref);
    clear data
    save(fullfile(lg_out, img_fn, [img_fn '_reg.mat']), 'img_ref', 'outs');
    nmask = size(mask3D,3);
    FrameRate = 30;
    fprintf('Getting TCs \n')
    tc_avg = getTC(img_reg, mask3D, nmask);

    saveData = 1;
    reg_sum = sum(img_reg,3);
    %plot TCs with and without the frames without laser power
    plotTC(tc_avg, mask3D, reg_sum, 1:size(tc_avg,2), FrameRate, fullfile(lg_out,img_fn), saveData);
    %plotTC(tc_avg(laser_on_ind,:), mask3D, reg_sum, 1:size(tc_avg,2),
    %FrameRate, fullfile(lg_out,img_fn), 0); laser frames?
    mask_flat = plotMask(mask3D, saveData, fullfile(lg_out,img_fn),1);
    title([date ' ' mouse])
    print(fullfile(lg_out,img_fn, [img_fn '_finalMask.pdf']),'-dpdf')
    data_corr = corrcoef(tc_avg);
    figure; fig = imagesc(data_corr); colorbar
    title([date ' ' mouse])
    print(fullfile(lg_out,img_fn, [img_fn '_data_corr.pdf']),'-dpdf')

    
    save(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']), 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
    
    %% make movie
    fprintf('Making movie \n')
    sz = size(img_reg);    
    cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
    nTrials = length(cTargetOn);
    prewin_frames = round(1500./mworks.frameRateHz);
    postwin_frames = round(3000./mworks.frameRateHz);
    reg_align = nan(sz(1),sz(2),prewin_frames+postwin_frames,nTrials);
    for itrial = 1:nTrials
        if cTargetOn(itrial)+postwin_frames<size(img_reg,3)
            reg_align(:,:,:,itrial) = img_reg(:,:,cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1);
        end
    end
    reg_align_avg = nanmean(reg_align,4);
    writetiff(reg_align_avg(:,:,1:2:end), fullfile('\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper', expt(id).name, 'CC_movies', ['Day' num2str(id)], [img_fn '_cueAlign_rereg.tif']))

end
        
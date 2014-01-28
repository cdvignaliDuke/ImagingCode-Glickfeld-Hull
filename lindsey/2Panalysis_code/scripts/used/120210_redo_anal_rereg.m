P = 2;
nON = 12;
nOFF = 12;
nPlanes = 1;
begin = 1;
TFSFetc = [1:2];
pre_win = [7 12];
post_win = [13 24];
Nshuf = 500;
SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;

dSF = median(diff(log2(SF_vec0)));
dTF = median(diff(log2(TF_vec0)));
SF_vec00 = log2(SF_vec0(1)):(dSF/10):log2(SF_vec0(end));
TF_vec00 = log2(TF_vec0(1)):(dTF/10):log2(TF_vec0(end));
[sfsf00,tftf00]=meshgrid(SF_vec00,TF_vec00);
grid2.sfsf00 = sfsf00;
grid2.tftf00 = tftf00;

% redo_exp.date = {'110920' '110817' '110820' '110819' '110817'};
% redo_exp.mouse = {'DR9' 'AC44' 'AC45' 'Y26' 'AC44'};
% redo_exp.userun = {[1:4] [4:6] [1:4] [3:5] [1:3]};
% redo_exp.zoom = {1.5 1 1 1.5 1};
% redo_exp.count_prot = {1 2 1 2 1};
% redo_exp.dirs = {2 2 2 2 2};
% redo_exp.run = {0 0 0 0 0};
% redo_exp.blanks = {1 1 1 1 1};
% redo_exp.area = {'PM' 'AL' 'LM' 'LM' 'LM'};

redo_exp.date = {'111125' '110916'};
redo_exp.mouse = {'Y28' 'Y27'};
redo_exp.userun = {[1:2] [3:5]};
redo_exp.zoom = {1 1.5};
redo_exp.count_prot = {1 2};
redo_exp.dirs = {1 2};
redo_exp.run = {0 0};
redo_exp.blanks = {1 1};
redo_exp.area = {'PM' 'V1'};

for iexp = 1:2
    mouse = char(redo_exp.mouse{iexp});
    date = char(redo_exp.date{iexp});
    userun = redo_exp.userun{iexp};
    count_protocol = redo_exp.count_prot{iexp};
    run = redo_exp.dirs{iexp};
    blanks = redo_exp.blanks{iexp};
    dirs = redo_exp.dirs{iexp};
    zoom = redo_exp.zoom{iexp};
    
    if dirs ==1
        nCond = 25;
    elseif dirs ==2
        nCond = 50;
    end

    base = 'G:\users\lindsey\analysisLG\active mice';
    outDir = fullfile(base, mouse, date);

    fn_Big_seq = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat']);
    if exist(fn_Big_seq,'file') == 0
        eval(['PARAMS_' date '_' mouse])
        resort_seq_only
    else
        load(fn_Big_seq)
    end
        
    stack_sort
    siz = size(stack_sorted);
    stack_sorted_np = zeros(size(stack_sorted));
    avg = mean(mean(mean(stack_sorted,3),2),1);
    for iframe = 1:siz(3)
        np_av = mean(mean(stack_sorted(:,:,iframe),2),1);
        stack_sorted_np(:,:,iframe) = stack_sorted(:,:,iframe)-np_av+avg;
    end

    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted_npsub.tif']);
    writetiff(stack_sorted_np, fn_stack);

    clear stack_sorted

    fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn);

    siz = size(stack_sorted_np);
    stack_dF_np = zeros(siz(1), siz(2), nCond+1);
    start = 0;
    for iCond = 1:nCond+1;
        nRep = length(Big_Seqposition(iCond).ind);
        rep_dF = zeros(siz(1), siz(2),nRep);
        for iRep = 1:nRep
            rep_base = mean(stack_sorted_np(:,:,start+pre_win(1):start+pre_win(2)),3);
            rep_resp = mean(stack_sorted_np(:,:,start+post_win(1):start+post_win(2)),3);
            rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
            start = ((nOFF+nON)/nPlanes)+start;
        end
        stack_dF_np(:,:,iCond) = mean(rep_dF,3);
    end
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
    writetiff(stack_dF_np, fn_stack);

    nReps = sum(stim_reps(1,:));
    if dirs == 2
        add = 1;
        stim_reps_dir = zeros(1,26);
        for iCond = 1:nCond/2
        stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
        add = add+2;
        end
        stim_reps_dir(1,end) = stim_reps(1,end);
        stim_reps = stim_reps_dir;
    end

    roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
    start = 1;
    rep = 1;
    for iCond = 1:length(stim_reps); 
        nRep = stim_reps(iCond);
        for iRep = 1:nRep;
            roi_stim(:,:,:, rep) = stack_sorted_np(:,:,start:start-1+nON+nOFF);
            start = start+nON+nOFF;
            rep = rep+1;
        end
    end
    resp_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
    resp_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));

    clear stack_sorted_np
    clear roi_stim

    alphaB = .05./(25);
    Info_ttest_mat = zeros(siz(1),siz(2),25);

    b= 5;

    f1 = fspecial('average');
    siz = size(resp_on);
    resp_on_long = reshape(resp_on, [siz(1) siz(2)*siz(3)]);
    resp_off_long = reshape(resp_off, [siz(1) siz(2)*siz(3)]);

    clear resp_on
    clear resp_off

    resp_on_sm = reshape(filter2(f1,resp_on_long),[siz(1) siz(2) siz(3)]);
    resp_off_sm = reshape(filter2(f1,resp_off_long),[siz(1) siz(2) siz(3)]);

    clear resp_on_long
    clear resp_off_long

    for iy = b+1:240-b
        fprintf([num2str(iy) ' '])
        for ix = b+1:256-b
            start = 1;
            p_ttestB = zeros(1,1,25);
            for iCond = 1:(length(stim_reps)-1)
                nRep = stim_reps(1,iCond);
                [h_ttestB1,p_ttestB1] = ttest(resp_off_sm(iy,ix,start:start-1+nRep),resp_on_sm(iy,ix,start:start-1+nRep),alphaB,'left');
                p_ttestB(1,1,iCond) = p_ttestB1;
                start = start+nRep;
            end
        Info_ttest_mat(iy,ix,:) = p_ttestB;
        end
    end

    clear resp_on_sm
    clear resp_off_sm

    siz = size(Info_ttest_mat);
    Info_ttest_mat_long = reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]);
    ttest_smooth = reshape(filter2(f1,Info_ttest_mat_long), [siz(1) siz(2) siz(3)]);
    ttest_mask = min(ttest_smooth,[],3) < alphaB;

    ttest_mask(1:5,1:end) = 0;
    ttest_mask(1:end, 1:5) = 0;
    ttest_mask(1:end, 251:end) = 0;
    ttest_mask(235:end,1:end) = 0;

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    save(fn_out, 'ttest_mask', 'Info_ttest_mat');
    
        
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
    stack_sorted = readtiff(fn_out);
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn_out);

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
    stack_dF = readtiff(fn_out);

    siz = size(stack_sorted);

    nReps = sum(stim_reps(1,:));

    roi_stim = zeros(siz(1), siz(2), nOFF+nON,nReps);
    start = 1;
    rep = 1;
    for iRep = 1:nReps;
        roi_stim(:,:,:, rep) = stack_sorted(:,:,start:start-1+nON+nOFF);
        start = start+nON+nOFF;
        rep = rep+1;
    end
    stim_off = squeeze(mean(roi_stim(:,:,pre_win(1):pre_win(2),:),3));
    stim_on = squeeze(mean(roi_stim(:,:,post_win(1):post_win(2),:),3));

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
    save(fn_out, 'stim_on', 'stim_off');

    ttest_mask_long = reshape(ttest_mask, [siz(1)*siz(2), 1]);
    ind = find(ttest_mask_long == 1);
    ttest_mask_all(:,:,iexp) = ttest_mask;

    stack_sorted_long = reshape(stack_sorted, [siz(1)*siz(2) siz(3)]);
    clear stack_sorted
    stack_avg = mean(stack_sorted_long(ind,:),1);
    clear stack_sorted_long;

    if dirs ==1
        nCond = 25;
    elseif dirs ==2
        nCond = 50;
    end

    siz = size(stack_dF);
    stack_dF_long = reshape(stack_dF, [siz(1)*siz(2) siz(3)]);
    stack_dF_avg = mean(stack_dF_long(ind,:),1); 
    [peak_dF peak_cond] = max(stack_dF_avg,[],2);
    nRep = stim_reps(1,peak_cond);
    start_rep = sum(stim_reps(1,1:(peak_cond-1)));

    ncyc = size(stack_avg,2)/(nOFF+nON);
    stdev= zeros(nRep, 1);
    start = 1;
    for iRep = start_rep+1:nRep+start_rep;
        stdev(start,:) = std(stack_avg(:,pre_win(1)+((iRep-1)*(nON+nOFF)):(pre_win(2)-1+((iRep-1)*(nON+nOFF)))));
        start = 1+start;
    end
    noise = mean(stdev,1);

    clear stack_avg

    resp_on_long = reshape(stim_on, [siz(1)*siz(2) ncyc]);
    resp_off_long = reshape(stim_off, [siz(1)*siz(2) ncyc]);

    clear stim on
    clear stimoff

    resp_dF_peak = resp_on_long(ind,start_rep+1:nRep+start_rep)- resp_off_long(ind,start_rep+1:nRep+start_rep);
    resp_dF_avg = mean(mean(resp_dF_peak,1),2);
    SNR = resp_dF_avg./noise;

    clear resp_on_long
    clear resp_off_long

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_SNR.mat']);
    save(fn_out, 'SNR', 'resp_dF_avg', 'noise');
    
    fn_reps= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
    fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
    stack_dF = readtiff(fn_stack); 

    siz = size(stack_dF);
    stack_max= max(stack_dF,[],3);
    max_interp = interp2(stack_max);
    f1 = fspecial('average');   
    max_interp_sm = filter2(f1, max_interp);
    siz2 = size(max_interp);
    Xi = 1:2:siz2(1);
    Yi = 1:2:siz2(2);
    stack_max_interp_sm = interp2(max_interp_sm, Yi', Xi);
    stack_max_sm_long = reshape(stack_max_interp_sm,[siz(1)*siz(2) 1]);


    local_max = zeros(siz(1), siz(2));
    border = 5;
    for iy = border:(siz(1)-border);
        for ix = border:(siz(2)-border);            
            sub = stack_max_interp_sm(iy-1:iy+1,ix-1:ix+1);
            sub_long = reshape(sub, [1 9]);
            [sub_long_order ind_order] = sort(sub_long);
            if ind_order(end)==5
                local_max(iy,ix) = 1;
            end
        end
    end
    local_max_long = reshape(local_max, [siz(1)*siz(2) 1]);
    ind_local_max = find(local_max_long==1);

    fn_ttest= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_ttest_mask.mat']);
    load(fn_ttest);
    siz = size(Info_ttest_mat);
    Info_ttest_mat_long = interp2(reshape(Info_ttest_mat, [siz(1) siz(2)*siz(3)]));
    Info_ttest_smooth = filter2(f1,Info_ttest_mat_long);
    siz_interp = size(Info_ttest_smooth);
    Xi = 1:2:siz_interp(1);
    Yi = 1:2:siz_interp(2);
    ttest_smooth_siz = interp2(Info_ttest_smooth, Yi', Xi);
    ttest_smooth = min(reshape(ttest_smooth_siz, [siz(1) siz(2) siz(3)]),[],3);
    ttest_long = reshape(ttest_smooth, [siz(1)*siz(2) 1]);

    ind_highP = find(ttest_long(ind_local_max,:)>=(.05./25));
    local_max_long(ind_local_max(ind_highP,:),:) = 0;

    local_max = reshape(local_max_long, [siz(1) siz(2)]);
    n_pix = sum(sum(local_max));
    [i, j] = find(local_max ==1);

    fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
    save(fn_local, 'local_max', 'n_pix', 'i', 'j');

    fn_stim = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
    load(fn_stim);

    nRep = stim_reps(1,1);
    stack_var=var(stim_on,[],3);
    stack_mean=mean(stim_on,3);
    stack_mean_long = reshape(stack_mean, [siz(1)*siz(2) 1]);
    stack_var_long = reshape(stack_var, [siz(1)*siz(2) 1]);
    b=robustfit(stack_mean_long,stack_var_long);

    nPhoton=round((stack_mean_long./b(2)).*nRep.*(nON/2.667));
    fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_photon.mat']);
    save(fn_stim, 'b', 'nPhoton');

    resp_dF = zeros(n_pix, size(stim_on,3));
    resp_dFoverF = zeros(n_pix, size(stim_on,3));
    for ipix = 1:n_pix
        sub_y = [i(ipix)-1:i(ipix)+1]; 
        sub_x = [j(ipix)-1:j(ipix)+1];
        roi_on = squeeze(mean(mean(stim_on(sub_y,sub_x,:),2),1));
        roi_off = squeeze(mean(mean(stim_off(sub_y,sub_x,:),2),1));
        resp_dF(ipix,:) = (roi_on-roi_off)';
        resp_dFoverF(ipix,:) = ((roi_on-roi_off)./roi_off)';
    end
    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
    save(fn_out, 'resp_dF', 'resp_dFoverF');
    clear stim_off
    clear stim_on
    clear roi_stim

    nReps = sum(stim_reps(1,:));
    if dirs == 2
        add = 1;
        stim_reps_dir = zeros(1,26);
        for iCond = 1:nCond/2
        stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
        add = add+2;
        end
        stim_reps_dir(1,end) = stim_reps(1,end);
        stim_reps = stim_reps_dir;
    end

    Ind_struct = [];
    start = 1;
    for iCond = 1:25
        nRep = stim_reps(iCond);
        Ind_struct(iCond).all_trials = [start:start-1+nRep];
        start = start+nRep;
    end

    Fit_struct = [];
    for count_shuf = 0:Nshuf
        fprintf('.')
        Im_mat_USE = zeros(size(resp_dFoverF,1), 25);
        Im_mat_std = zeros(size(resp_dFoverF,1), 25);
        dF_mat = zeros(size(resp_dFoverF,1), 25);
        for iCond = 1:25        
            ind_all = Ind_struct(iCond).all_trials;
            if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
            else
                ind_all_1 = ind_all;        
            end
            Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
            dF_mat(:,iCond) = mean(resp_dF(:,ind_all_1),2);
        end

        start = 1;
        for iCell = 1:n_pix;
            a = Im_mat_USE(iCell,:);
            if max(a,[],2) > 0     
                b = reshape(a',length(SF_vec0),length(TF_vec0));
                %b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
                data = b';
                ind0 = find(data<0);
                data(ind0) = NaN;
                if count_shuf == 0
                    PLOTIT_FIT = 0;
                    SAVEALLDATA = 1;
                    Fit_2Dellipse_LG
                    eval(['Fit_struct(iCell).True.s_',' = s;']);
                else
                    SAVEALLDATA = 0;
                    PLOTIT_FIT = 0;
                    Fit_2Dellipse_LG
                    eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
                end
            end               
        end
    end

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
    save(fn_out, 'Fit_struct')


    if Nshuf>1;
        for iCell = 1:n_pix
            if ~isempty(Fit_struct(iCell).True)                
                eval(['tmp = Fit_struct(iCell).True.s_.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_10];']);
                fit_true_vec(iCell,:) = tmp;
            end
        end

        for count_shuf = 1:Nshuf
            for iCell = 1:n_pix
                if ~isempty(Fit_struct(iCell).Shuf)
                    eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_50];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_10];']);
                    eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_10];']);
                    %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                    fit_shuf_vec(iCell,:,count_shuf) = tmp;
                end
            end
        end

        Npars = size(fit_shuf_vec,2);
        lbub_fits = zeros(n_pix,Npars,5);
        alpha_bound = .025;
        for iCell = 1:n_pix
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
    for iCell = 1:n_pix
        if lbub_diff(iCell,4)<2 
            if lbub_diff(iCell,5)<2
                goodfit_ind = [goodfit_ind iCell];
            end
        end
    end

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);   
    save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind')
    
    fn_stim = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
    load(fn_stim);

    fn_reps = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn_reps);

    fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
    load(fn_resp);
    i=[];
    j=[];
    local_max = [];
    fn_local = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
    load(fn_local);

    bouton_map = zeros(size(local_max));
    for ipix = 1:n_pix
        sub_y = i(ipix)-1:i(ipix)+1; 
        sub_x = j(ipix)-1:j(ipix)+1;
        bouton_map(sub_y, sub_x) = ipix;
    end

    siz = size(bouton_map); 
    npmasks = zeros(siz(1),siz(2),n_pix);
    mask_dil = zeros(siz(1),siz(2),n_pix);
    for ipix = 1:n_pix
        mask_dil(i(ipix)-2:i(ipix)+2, j(ipix)-2:j(ipix)+2,ipix) = 1;
    end
    all_mask_dil = sum(mask_dil,3);

    for ipix = 1:n_pix
        np_plus = zeros(siz(1),siz(2));
        np_plus(i(ipix)-4:i(ipix)+4, j(ipix)-4:j(ipix)+4) = 1;
        npmask = np_plus-mask_dil(:,:,ipix);
        npmask(find(all_mask_dil>0)) = 0;
        npmask(1:5,:)=0;
        npmask(236:240,:)=0;
        npmask(:,1:5)=0;
        npmask(:,252:256)=0;
        npmasks(:,:,ipix) = npmask;
    end

    npmask_sum = sum(npmasks,3);
    npmask_all = zeros(size(npmask_sum));
    npmask_all(find(npmask_sum>0))=1;
    ind_all = find(npmask_all==1);

    siz = size(stim_on); 
    stim_on_long = reshape(stim_on, siz(1)*siz(2), siz(3));
    stim_off_long = reshape(stim_off, siz(1)*siz(2), siz(3));
    stim_on_npmask_all = mean(stim_on_long(ind_all,:),1);    
    stim_off_npmask_all = mean(stim_off_long(ind_all,:),1);
    nponly_dFoverF = (stim_on_npmask_all-stim_off_npmask_all)./stim_off_npmask_all;
    np_all_avg = mean([stim_on_npmask_all stim_off_npmask_all],2);

    resp_dF_np = zeros(n_pix, siz(3));
    resp_dFoverF_np = zeros(n_pix, siz(3));
    roi_on_np = zeros(n_pix, sum(stim_reps));
    roi_off_np = zeros(n_pix, sum(stim_reps));
    for ipix = 1:n_pix
        sub_y = [i(ipix)-1:i(ipix)+1]; 
        sub_x = [j(ipix)-1:j(ipix)+1];
        stim_on_long = reshape(stim_on, siz(1)*siz(2),siz(3));
        stim_off_long = reshape(stim_off, siz(1)*siz(2),siz(3));
        ind = find(npmasks(:,:,ipix)==1);
        stim_on_npmask = mean(stim_on_long(ind,:),1);    
        stim_off_npmask = mean(stim_off_long(ind,:),1);
        np_avg = mean([stim_on_npmask stim_off_npmask],2);
        np_on = stim_on_npmask_all*(np_avg/np_all_avg);
        np_off = stim_off_npmask_all*(np_avg/np_all_avg);
        roi_on_np(ipix,:) = bsxfun(@minus, squeeze(mean(mean(stim_on(sub_y,sub_x,:),1),2)), np_on') +np_avg;
        roi_off_np(ipix,:) = bsxfun(@minus, squeeze(mean(mean(stim_off(sub_y,sub_x,:),1),2)), np_off') +np_avg;
        resp_dF_np(ipix,:) = (roi_on_np(ipix,:)-roi_off_np(ipix,:))';
        resp_dFoverF_np(ipix,:) = ((roi_on_np(ipix,:)-roi_off_np(ipix,:))./roi_off_np(ipix,:))';
    end

    nReps = sum(stim_reps(1,:));
    if dirs == 2
        add = 1;
        stim_reps_dir = zeros(1,26);
        for iCond = 1:nCond/2
        stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
        add = add+2;
        end
        stim_reps_dir(1,end) = stim_reps(1,end);
        stim_reps = stim_reps_dir;
    end        

    Ind_struct = [];
    start = 1;
    for iCond = 1:25
        nRep = stim_reps(iCond);
        Ind_struct(iCond).all_trials = [start:start-1+nRep];
        start = start+nRep;
    end

    Im_mat_nponly = zeros(1, 25);
    Im_mat_np = zeros(n_pix, 25);
    Im_mat_USE = zeros(n_pix, 25);
    for iCond = 1:25        
        ind_all = Ind_struct(iCond).all_trials;
        Im_mat_nponly(:,iCond) = mean(nponly_dFoverF(:,ind_all),2);
        Im_mat_np(:,iCond) = mean(resp_dFoverF_np(:,ind_all),2);
        Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all),2);
    end

    fn_neuropil = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_neuropil.mat']);
    save(fn_neuropil, 'npmasks', 'npmask_all', 'resp_dFoverF_np', 'nponly_dFoverF', 'Im_mat_nponly', 'Im_mat_np', 'Im_mat_USE');

end


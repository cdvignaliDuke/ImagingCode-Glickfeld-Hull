clear all
areas = ['RL'; 'AM'];
inj = 'V1';
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
for iArea = 1:2
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};

        if dirs ==1
            nCond = 25;
        elseif dirs ==2
            nCond = 50;
        end

        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
     
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
        alpha = 0.05/nCond;
        ind_highP = find(ttest_long(ind_local_max,:)>=(alpha));
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
    end
end
clear all
areas = ['PM'; 'LM'; 'AL'];
inj = 'V1';
P = 2;
nPlanes = 1;
begin = 1;

for iArea = 2;
    matrix = 'DIR8';
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
        SFs = exp_list.SF_mat{iexp};
        TFs = exp_list.TF_mat{iexp};
        pair_run = exp_list.paired_mat{iexp};
        nframes = exp_list.nframes_mat{iexp};
        
        nON = nframes;
        nOFF = nframes;
        if nframes == 12;
            pre_win = [7 12];
            post_win = [13 24];
        elseif nframes == 15;
            pre_win = [8 15];
            post_win = [16 30];
        end
        
        nSFTF = length(SFs).*length(TFs);
        nCond = 8.*nSFTF;
        
        ORI_vec0 = [0 45 90 135 180 225 270 315]*pi/180;
        
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

        ind_highP = find(ttest_long(ind_local_max,:)>=(.05./25));
        local_max_long(ind_local_max(ind_highP,:),:) = 0;
        
        local_max = reshape(local_max_long, [siz(1) siz(2)]);
        n_pix = sum(sum(local_max));
        [i, j] = find(local_max ==1);
        
        fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        save(fn_local, 'local_max', 'n_pix', 'i', 'j');
        
        fn_stack = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
        stack_sorted = readtiff(fn_stack);
        
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
    end
end

clear all
areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);
    nON = 12;
    nOFF = 12;
    nPlanes = 1;

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    for iexp = 1:nexp;

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
end



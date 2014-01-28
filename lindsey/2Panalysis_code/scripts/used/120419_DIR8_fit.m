clear all
areas = ['PM'; 'LM'; 'AL'];
inj = 'V1';
P = 2;
nPlanes = 1;
begin = 1;
Nshuf = 500;

for iArea = 1:3;
    matrix = 'DIR16';
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
        dirs = exp_list.dirs_mat{iexp};
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
        SF_vec0 = SFs;
        TF_vec0 = TFs;
        if nSFTF == 4;
            sfsf_grid = [SFs SFs];
            tftf_grid = [TFs(1) TFs(1) TFs(2) TFs(2)];
        elseif nSFTF<4;
            sfsf_grid = repmat(SFs, 1, length(TFs));
            tftf_grid = repmat(TFs, 1, length(SFs));
        end
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_reps= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
        fn_local= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);
        [i, j] = find(local_max ==1);
        
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_out);
        
        Ind_struct = [];
        start = 1;
        for iCond = 1:nCond
            nRep = stim_reps(iCond);
            Ind_struct(iCond).all_trials = [start:start-1+nRep];
            start = start+nRep;
        end

        Fit_struct = [];
        for count_shuf = 0:Nshuf
            fprintf('.')
            Im_mat_USE = zeros(size(resp_dFoverF,1), nCond);
            Im_mat_std = zeros(size(resp_dFoverF,1), nCond);
            dF_mat = zeros(size(resp_dFoverF,1), nCond);
            for iCond = 1:nCond       
                ind_all = Ind_struct(iCond).all_trials;
                if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                    ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
                else
                    ind_all_1 = ind_all;        
                end
                Im_mat_USE(:,iCond) = mean(resp_dFoverF(:,ind_all_1),2);
                Im_mat_std(:,iCond) = std(resp_dFoverF(:,ind_all_1),[],2);
                dF_mat(:,iCond) = mean(resp_dF(:,ind_all_1),2);
            end

            start = 1;
            for iCell = 1:n_pix;
                a = Im_mat_USE(iCell,:);
                b = Im_mat_std(iCell,:);
                if max(a,[],2) > 0     
                    data = reshape(a',nSFTF, length(ORI_vec0))';
                    data_std = reshape(b',nSFTF, length(ORI_vec0))';
                    if count_shuf == 0
                        SAVEALLDATA = 1;
                        Fit_Ori_LG
                        eval(['Fit_struct(iCell).True.s_',' = s;']);
                    else
                        SAVEALLDATA = 0;
                        Fit_Ori_LG
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
                    fit_true_vec(iCell,:) = tmp;
                end
            end

            for count_shuf = 1:Nshuf
                for iCell = 1:n_pix
                    if ~isempty(Fit_struct(iCell).Shuf)
                        eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                        %fit is: OI (p-n)/p; Ori vector tuning; OI2 (p-n)/(p+n); Ori peak angle; Ori vector angle;
                        % DI (p-n)/p; Dir vector tuning; DI2 (p-n)/(p+n); Dir peak angle; Dir vector angle
                        fit_shuf_vec(iCell,:,count_shuf) = tmp;
                    end
                end
            end

            Npars = size(fit_shuf_vec,2);
            lbub_fits = zeros(n_pix,Npars,6);
            alpha_bound = .025;
            for iCell = 1:n_pix
                lbub_fits(iCell,5,6) = circ_confmean(squeeze(fit_shuf_vec(iCell,5,:)));
                lbub_fits(iCell,10,6) = circ_confmean(squeeze(fit_shuf_vec(iCell,10,:)));
            end
           
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

        goodfit_ind = [];
        goodfit_ind_tuned = [];
        goodfit_ind_untuned = [];
        for iCell = 1:n_pix
            if lbub_fits(iCell,5,6)<0.05 
                goodfit_ind = [goodfit_ind iCell];
                if Fit_struct(iCell).True.s_.x(3)>0.33
                    goodfit_ind_tuned = [goodfit_ind iCell];
                end
            elseif lbub_fits(iCell,5,6)>0.05 
                if mean(Fit_struct(iCell).True.s_.data3,2)>0.05
                    if Fit_struct(iCell).True.s_.x(3)<0.33
                        if lbub_fits(iCell,3,2)<0.5
                            goodfit_ind = [goodfit_ind iCell];
                            goodfit_ind_untuned = [goodfit_ind_untuned iCell];
                        end
                    end                       
                end
            end
        end

        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);   
        save(fn_out, 'lbub_fits', 'goodfit_ind', 'goodfit_ind_tuned','goodfit_ind_untuned')
        end
    end
end
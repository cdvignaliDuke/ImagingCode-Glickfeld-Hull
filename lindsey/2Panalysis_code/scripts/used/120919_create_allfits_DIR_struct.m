clear all
areas = ['PM'; 'LM'; 'AL'];
P = 2;
matrix = 'DIR';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
all_fits = [];
for iArea = 1:3;
    image = areas(iArea,:);
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    mouse_list = [];
    mice = exp_list.mouse_mat;
    all_fits_dir(iArea).nexp = nexp;
    all_fits_dir(iArea).name = image;
    all_fits_dir(iArea).x_bestSFTF = [];
    all_fits_dir(iArea).x_bestSFTF_orituned = [];
    all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd = [];
    all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd = [];
    all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd = [];
    all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd = [];
    all_fits_dir(iArea).x_2Hz_pt16cpd = [];
    all_fits_dir(iArea).x_2Hz_pt04cpd = [];
    all_fits_dir(iArea).x_8Hz_pt16cpd = [];
    all_fits_dir(iArea).x_8Hz_pt04cpd = [];
    all_fits_dir(iArea).x_orituned_2Hz_pt16cpd = [];
    all_fits_dir(iArea).x_orituned_2Hz_pt04cpd = [];
    all_fits_dir(iArea).x_orituned_8Hz_pt16cpd = [];
    all_fits_dir(iArea).x_orituned_8Hz_pt04cpd = [];
    all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd = [];
    all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd = [];
    all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd = [];
    all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd = [];
    all_fits_dir(iArea).data_bestSFTF = [];
    all_fits_dir(iArea).data_bestSFTF_orituned = [];
    all_fits_dir(iArea).data_bestSFTF_2Hz_pt16cpd = [];
    all_fits_dir(iArea).data_bestSFTF_2Hz_pt04cpd = [];
    all_fits_dir(iArea).data_bestSFTF_8Hz_pt16cpd = [];
    all_fits_dir(iArea).data_bestSFTF_8Hz_pt04cpd = [];
    all_fits_dir(iArea).data_2Hz_pt16cpd = [];
    all_fits_dir(iArea).data_2Hz_pt04cpd = [];
    all_fits_dir(iArea).data_8Hz_pt16cpd = [];
    all_fits_dir(iArea).data_8Hz_pt04cpd = [];
    all_fits_dir(iArea).data_orituned_2Hz_pt16cpd = [];
    all_fits_dir(iArea).data_orituned_2Hz_pt04cpd = [];
    all_fits_dir(iArea).data_orituned_8Hz_pt16cpd = [];
    all_fits_dir(iArea).data_orituned_8Hz_pt04cpd = [];
    all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt16cpd = [];
    all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt04cpd = [];
    all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt16cpd = [];
    all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt04cpd = [];
    all_fits_dir(iArea).x_orituned_2Hz_bestSF = [];
    all_fits_dir(iArea).x_orituned_8Hz_bestSF = [];
    all_fits_dir(iArea).x_2Hz_bestSF = [];
    all_fits_dir(iArea).x_8Hz_bestSF = [];
    all_fits_dir(iArea).data_orituned_2Hz_bestSF = [];
    all_fits_dir(iArea).data_orituned_8Hz_bestSF= [];
    all_fits_dir(iArea).data_2Hz_bestSF = [];
    all_fits_dir(iArea).data_8Hz_bestSF= [];
    all_fits_dir(iArea).x_orituned_pt04cpd_bestTF = [];
    all_fits_dir(iArea).x_orituned_pt16cpd_bestTF= [];
    all_fits_dir(iArea).x_pt04cpd_bestTF= [];
    all_fits_dir(iArea).x_pt16cpd_bestTF= [];
    all_fits_dir(iArea).data_orituned_pt04cpd_bestTF= [];
    all_fits_dir(iArea).data_orituned_pt16cpd_bestTF= [];
    all_fits_dir(iArea).data_pt04cpd_bestTF= [];
    all_fits_dir(iArea).data_pt16cpd_bestTF= [];
    
    for iexp = 1:nexp;
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_prot = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        SFs = exp_list.SF_mat{iexp};
        TFs = exp_list.TF_mat{iexp};
        dirs = exp_list.dirs_mat{iexp};
        pair_run = exp_list.paired_mat{iexp};
        nframes = exp_list.nframes_mat{iexp};
        zoom = exp_list.zoom_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        all_fits_dir(iArea).expt(iexp).mouse = mouse;
        all_fits_dir(iArea).expt(iexp).date = date;
        all_fits_dir(iArea).expt(iexp).userun = userun;
        all_fits_dir(iArea).expt(iexp).count_prot = count_prot;
        all_fits_dir(iArea).expt(iexp).dirs = dirs;
        all_fits_dir(iArea).expt(iexp).zoom = zoom;
        all_fits_dir(iArea).expt(iexp).nSFTF = length(SFs).*length(TFs);
        all_fits_dir(iArea).expt(iexp).nframes = nframes;
        all_fits_dir(iArea).expt(iexp).pair_run = pair_run;

        fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_lbub);
        fn_fit  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);
        load(fn_fit); 
        fn_local  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local); 
        fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_resp);
        fn_reps= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
        all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid = Fit_struct(1).True.s_.sfsf_tftf_grid;
        nReps = sum(stim_reps(1,:));
        nSFTF = length(SFs).*length(TFs);
        for iSFTF = 1:nSFTF
            all_fits(iArea).expt(iexp).n(:,iSFTF) = [n_pix length(goodfit_ind(iSFTF).all) length(goodfit_ind(iSFTF).orituned) length(goodfit_ind(iSFTF).oriuntuned) length(goodfit_ind(iSFTF).dirtuned)];
        end

        for iCell = 1:size(lbub_fits,1)
            all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI = Fit_struct(iCell).True.s_.x(:,3);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).ori_angle = Fit_struct(iCell).True.s_.x(:,5);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI = Fit_struct(iCell).True.s_.x(:,8);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).dir_angle = Fit_struct(iCell).True.s_.x(:,10);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).pos = [i(iCell, :) j(iCell, :)];
            all_fits_dir(iArea).expt(iexp).bouton(iCell).best_SF = Fit_struct(iCell).True.s_.best_SF;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).best_TF = Fit_struct(iCell).True.s_.best_TF;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF = Fit_struct(iCell).True.s_.ind_bestSFTF;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF = Fit_struct(iCell).True.s_.data;
            if dirs == 16
                data = zeros(8,nSFTF);
                start = 1;
                for idir = 1:8
                	data(idir,:) = mean(Fit_struct(iCell).True.s_.data(start:start+1,:),1);
                    start= start+2;
                end
                Fit_struct(iCell).True.s_.data = data;
            end
            for iSFTF = 1:nSFTF
                if find(goodfit_ind(iSFTF).all == iCell)
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit(iSFTF) = 1;
                else
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit(iSFTF) = 0;
                end
                if find(goodfit_ind(iSFTF).orituned == iCell)
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit_orituned(iSFTF) = 1;
                else
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit_orituned(iSFTF) = 0;
                end
            end
            
            if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                all_fits_dir(iArea).expt(iexp).bouton(iCell).orituned = 1;
                if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).dirtuned == iCell)
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).dirtuned = 1;
                else
                    all_fits_dir(iArea).expt(iexp).bouton(iCell).dirtuned = 0;
                end
            elseif find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).oriuntuned == iCell)
                all_fits_dir(iArea).expt(iexp).bouton(iCell).oriuntuned = 1;
            end
            if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).all == iCell)
                all_fits_dir(iArea).x_bestSFTF = [all_fits_dir(iArea).x_bestSFTF; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                all_fits_dir(iArea).data_bestSFTF = [all_fits_dir(iArea).data_bestSFTF Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                    all_fits_dir(iArea).x_bestSFTF_orituned = [all_fits_dir(iArea).x_bestSFTF_orituned; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                    all_fits_dir(iArea).data_bestSFTF_orituned = [all_fits_dir(iArea).data_bestSFTF_orituned Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                end
                if Fit_struct(iCell).True.s_.best_SF == 0.16
                    if Fit_struct(iCell).True.s_.best_TF == 2
                        all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd = [all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                        all_fits_dir(iArea).data_bestSFTF_2Hz_pt16cpd = [all_fits_dir(iArea).data_bestSFTF_2Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd = [all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                            all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt16cpd = [all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        end
                    elseif Fit_struct(iCell).True.s_.best_TF == 8
                        all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd = [all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                        all_fits_dir(iArea).data_bestSFTF_8Hz_pt16cpd = [all_fits_dir(iArea).data_bestSFTF_8Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd = [all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                            all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt16cpd = [all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        end
                    end
                elseif Fit_struct(iCell).True.s_.best_SF == 0.04
                    if Fit_struct(iCell).True.s_.best_TF == 2
                        all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd = [all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                        all_fits_dir(iArea).data_bestSFTF_2Hz_pt04cpd = [all_fits_dir(iArea).data_bestSFTF_2Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd = [all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                            all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt04cpd = [all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        end
                    elseif Fit_struct(iCell).True.s_.best_TF == 8
                        all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd = [all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                        all_fits_dir(iArea).data_bestSFTF_8Hz_pt04cpd = [all_fits_dir(iArea).data_bestSFTF_8Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        if find(goodfit_ind(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd = [all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd; Fit_struct(iCell).True.s_.x(all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF,:)];
                            all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt04cpd = [all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,all_fits_dir(iArea).expt(iexp).bouton(iCell).ind_bestSFTF)];
                        end
                    end
                end
                if find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,:) == 2);
                    x = find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,:) == 2);
                    max_dF= max(Fit_struct(iCell).True.s_.data,[],1);
                    [bestSF_2Hz  bestSF_2Hz_ind] = max(max_dF(:,x),[],2);
                    if find(goodfit_ind(x(bestSF_2Hz_ind)).all == iCell)
                        all_fits_dir(iArea).x_2Hz_bestSF = [all_fits_dir(iArea).x_2Hz_bestSF; Fit_struct(iCell).True.s_.x(x(bestSF_2Hz_ind),:)];
                        all_fits_dir(iArea).data_2Hz_bestSF = [all_fits_dir(iArea).data_2Hz_bestSF Fit_struct(iCell).True.s_.data(:,x(bestSF_2Hz_ind))];
                        if find(goodfit_ind(x(bestSF_2Hz_ind)).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_2Hz_bestSF = [all_fits_dir(iArea).x_orituned_2Hz_pt16cpd; Fit_struct(iCell).True.s_.x(x(bestSF_2Hz_ind),:)];
                            all_fits_dir(iArea).data_orituned_2Hz_bestSF = [all_fits_dir(iArea).data_orituned_2Hz_bestSF Fit_struct(iCell).True.s_.data(:,x(bestSF_2Hz_ind))];
                        end
                    end
                end
                if find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,:) == 8);
                    x = find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,:) == 8);
                    max_dF= max(Fit_struct(iCell).True.s_.data,[],1);
                    [bestSF_8Hz  bestSF_8Hz_ind] = max(max_dF(:,x),[],2);
                    if find(goodfit_ind(x(bestSF_8Hz_ind)).all == iCell)
                        all_fits_dir(iArea).x_8Hz_bestSF = [all_fits_dir(iArea).x_8Hz_bestSF; Fit_struct(iCell).True.s_.x(x(bestSF_8Hz_ind),:)];
                        all_fits_dir(iArea).data_8Hz_bestSF = [all_fits_dir(iArea).data_8Hz_bestSF Fit_struct(iCell).True.s_.data(:,x(bestSF_8Hz_ind))];
                        if find(goodfit_ind(x(bestSF_8Hz_ind)).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_8Hz_bestSF = [all_fits_dir(iArea).x_orituned_8Hz_pt16cpd; Fit_struct(iCell).True.s_.x(x(bestSF_8Hz_ind),:)];
                            all_fits_dir(iArea).data_orituned_8Hz_bestSF = [all_fits_dir(iArea).data_orituned_8Hz_bestSF Fit_struct(iCell).True.s_.data(:,x(bestSF_8Hz_ind))];
                        end
                    end
                end
                if find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,:) == 0.04);
                    x = find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,:) == 0.04);
                    max_dF= max(Fit_struct(iCell).True.s_.data,[],1);
                    [bestTF_pt04cpd  bestTF_pt04cpd_ind] = max(max_dF(:,x),[],2);
                    if find(goodfit_ind(x(bestTF_pt04cpd_ind)).all == iCell)
                        all_fits_dir(iArea).x_pt04cpd_bestTF = [all_fits_dir(iArea).x_pt04cpd_bestTF; Fit_struct(iCell).True.s_.x(x(bestTF_pt04cpd_ind),:)];
                        all_fits_dir(iArea).data_pt04cpd_bestTF = [all_fits_dir(iArea).data_pt04cpd_bestTF Fit_struct(iCell).True.s_.data(:,x(bestTF_pt04cpd_ind))];
                        if find(goodfit_ind(x(bestTF_pt04cpd_ind)).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_pt04cpd_bestTF = [all_fits_dir(iArea).x_orituned_pt04cpd_bestTF; Fit_struct(iCell).True.s_.x(x(bestTF_pt04cpd_ind),:)];
                            all_fits_dir(iArea).data_orituned_pt04cpd_bestTF = [all_fits_dir(iArea).data_orituned_pt04cpd_bestTF Fit_struct(iCell).True.s_.data(:,x(bestTF_pt04cpd_ind))];
                        end
                    end
                end
                if find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,:) == 0.16);
                    x = find(Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,:) == 0.16);
                    max_dF= max(Fit_struct(iCell).True.s_.data,[],1);
                    [bestTF_pt16cpd  bestTF_pt16cpd_ind] = max(max_dF(:,x),[],2);
                    if find(goodfit_ind(x(bestTF_pt16cpd_ind)).all == iCell)
                        all_fits_dir(iArea).x_pt16cpd_bestTF = [all_fits_dir(iArea).x_pt16cpd_bestTF; Fit_struct(iCell).True.s_.x(x(bestTF_pt16cpd_ind),:)];
                        all_fits_dir(iArea).data_pt16cpd_bestTF = [all_fits_dir(iArea).data_pt16cpd_bestTF Fit_struct(iCell).True.s_.data(:,x(bestTF_pt16cpd_ind))];
                        if find(goodfit_ind(x(bestTF_pt16cpd_ind)).orituned == iCell)
                            all_fits_dir(iArea).x_orituned_pt16cpd_bestTF = [all_fits_dir(iArea).x_orituned_pt16cpd_bestTF; Fit_struct(iCell).True.s_.x(x(bestTF_pt16cpd_ind),:)];
                            all_fits_dir(iArea).data_orituned_pt16cpd_bestTF = [all_fits_dir(iArea).data_orituned_pt16cpd_bestTF Fit_struct(iCell).True.s_.data(:,x(bestTF_pt16cpd_ind))];
                        end
                    end
                end
                for iSFTF = 1:nSFTF
                    if find(goodfit_ind(iSFTF).all == iCell)
                        if Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,iSFTF) == 0.16
                            if Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 2
                                all_fits_dir(iArea).x_2Hz_pt16cpd = [all_fits_dir(iArea).x_2Hz_pt16cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_2Hz_pt16cpd = [all_fits_dir(iArea).data_2Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 8
                                all_fits_dir(iArea).x_data_pt16cpd = [all_fits_dir(iArea).x_8Hz_pt16cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_8Hz_pt16cpd = [all_fits_dir(iArea).data_8Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            end
                        elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,iSFTF) == 0.04
                            if Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 2
                                all_fits_dir(iArea).x_2Hz_pt04cpd = [all_fits_dir(iArea).x_2Hz_pt04cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_2Hz_pt04cpd = [all_fits_dir(iArea).data_2Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 8
                                all_fits_dir(iArea).x_8Hz_pt04cpd = [all_fits_dir(iArea).x_8Hz_pt04cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_8Hz_pt04cpd = [all_fits_dir(iArea).data_8Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            end
                        end
                    end
                    if find(goodfit_ind(iSFTF).orituned == iCell)
                        if Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,iSFTF) == 0.16
                            if Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 2
                                all_fits_dir(iArea).x_orituned_2Hz_pt16cpd = [all_fits_dir(iArea).x_orituned_2Hz_pt16cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_orituned_2Hz_pt16cpd = [all_fits_dir(iArea).data_orituned_2Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 8
                                all_fits_dir(iArea).x_orituned_8Hz_pt16cpd = [all_fits_dir(iArea).x_orituned_8Hz_pt16cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_orituned_8Hz_pt16cpd = [all_fits_dir(iArea).data_orituned_8Hz_pt16cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            end
                        elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(1,iSFTF) == 0.04
                            if Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 2
                                all_fits_dir(iArea).x_orituned_2Hz_pt04cpd = [all_fits_dir(iArea).x_orituned_2Hz_pt04cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_orituned_2Hz_pt04cpd = [all_fits_dir(iArea).data_orituned_2Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            elseif Fit_struct(iCell).True.s_.sfsf_tftf_grid(2,iSFTF) == 8
                                all_fits_dir(iArea).x_orituned_8Hz_pt04cpd = [all_fits_dir(iArea).x_orituned_8Hz_pt04cpd; Fit_struct(iCell).True.s_.x(iSFTF,:)];
                                all_fits_dir(iArea).data_orituned_8Hz_pt04cpd = [all_fits_dir(iArea).data_orituned_8Hz_pt04cpd Fit_struct(iCell).True.s_.data(:,iSFTF)];
                            end
                        end
                    end
                end
            end                    
        end                            
    end
end
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits_dir.mat');
    save(fn_out, 'all_fits_dir');
        
        
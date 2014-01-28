clear all
areas = ['PM'; 'LM'; 'AL'];
P = 2;
matrix = 'DIR8';
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
    all_boutons = [0 0];
    
    for iexp = 1:nexp;
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_prot = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        SFs = exp_list.SF_mat{iexp};
        TFs = exp_list.TF_mat{iexp};
        pair_run = exp_list.paired_mat{iexp};
        nframes = exp_list.nframes_mat{iexp};
        zoom = exp_list.zoom_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        all_fits_dir(iArea).expt(iexp).mouse = mouse;
        all_fits_dir(iArea).expt(iexp).date = date;
        all_fits_dir(iArea).expt(iexp).userun = userun;
        all_fits_dir(iArea).expt(iexp).count_prot = count_prot;
        all_fits_dir(iArea).expt(iexp).dirs = mouse;
        all_fits_dir(iArea).expt(iexp).zoom = zoom;
        all_fits_dir(iArea).expt(iexp).nSFTF = SFs.*TFs;
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
        
        nReps = sum(stim_reps(1,:));
        
        all_fits(iArea).expt(iexp).n = [size(lbub_fits,1) length(goodfit_ind)];
        
        all_boutons = all_boutons+all_fits(iArea).expt(iexp).n;
        
        for iCell = 1:size(lbub_fits,1)
            all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI = Fit_struct(iCell).True.s_.x(:,1);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).ori_angle = Fit_struct(iCell).True.s_.x(:,2);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI = Fit_struct(iCell).True.s_.x(:,3);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).dir_angle = Fit_struct(iCell).True.s_.x(:,4);
            all_fits_dir(iArea).expt(iexp).bouton(iCell).pos = [i(iCell, :) j(iCell, :)];
            all_fits_dir(iArea).expt(iexp).bouton(iCell).best_SF = Fit_struct(iCell).True.s_.best_SF;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).best_TF = Fit_struct(iCell).True.s_.best_TF;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF = Fit_struct(iCell).True.s_.data;
            all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF_bestSFTF = Fit_struct(iCell).True.s_.data(:,Fit_struct(iCell).True.s_.ind_bestSFTF);
            if find(goodfit_ind == iCell)>0
                all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit = 1;
            else
                all_fits(iArea).expt(iexp).bouton(iCell).goodfit = 0;
            end
        end
    end
    all_fits(iArea).n = all_boutons;
end
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits_dir.mat');
    save(fn_out, 'all_fits_dir');
        
        
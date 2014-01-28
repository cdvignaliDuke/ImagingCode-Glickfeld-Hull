clear all
areas = ['PM'; 'LM'; 'AL'; 'RL'; 'AM'];
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
for iArea = 1:5
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
        
        if run == 1
            if dirs ==1
                nCond = 25;
            elseif dirs ==2
                nCond = 50;
            end

            base = 'G:\users\lindsey\analysisLG\active mice';    
            outDir = fullfile(base, mouse,date);

            fn_reps= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
            load(fn_reps);
            
            fn_wheel= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_wheel.mat']);
            load(fn_wheel);
            
            fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
            load(fn_local)

            fn_stim = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
            load(fn_stim);

            fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
            load(fn_resp);

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
                ind_run = find(Wheel_mat(start:start-1+nRep,:)>0);
                ind_norun = find(Wheel_mat(start:start-1+nRep,:)==0);
                Ind_struct(iCond).run_trials = start-1 + ind_run;
                Ind_struct(iCond).norun_trials = start-1 + ind_norun;
                start = start+nRep;
            end
            Fit_struct_run = [];
            for count_shuf = 0:Nshuf
                fprintf('.')
                Im_mat_USE_run = zeros(size(resp_dFoverF,1), 25);
                Im_mat_USE_norun = zeros(size(resp_dFoverF,1), 25);
                dF_mat_run = zeros(size(resp_dFoverF,1), 25);
                dF_mat_norun = zeros(size(resp_dFoverF,1), 25);
                for iCond = 1:25        
                    ind_run = Ind_struct(iCond).run_trials;
                    ind_norun = Ind_struct(iCond).norun_trials;
                    if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                        ind_run_1 = ind_run(randsample(length(ind_run),length(ind_run),1));
                        ind_norun_1 = ind_norun(randsample(length(ind_norun),length(ind_norun),1));
                    else
                        ind_run_1 = ind_run;
                        ind_norun_1 = ind_norun;
                    end
                    Im_mat_USE_run(:,iCond) = mean(resp_dFoverF(:,ind_run_1),2);
                    Im_mat_USE_norun(:,iCond) = mean(resp_dFoverF(:,ind_norun_1),2);
                    dF_mat_run(:,iCond) = mean(resp_dF(:,ind_run_1),2);
                    dF_mat_norun(:,iCond) = mean(resp_dF(:,ind_norun_1),2);
                end

                start = 1;
                for iCell = 1:n_pix;
                    for run = 1:2
                        if run == 1
                            a = Im_mat_USE_run(iCell,:);
                            dF_mat = dF_mat_run;
                        elseif run == 2
                            a = Im_mat_USE_norun(iCell,:);
                            dF_mat = dF_mat_norun;
                        end
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
                                if run == 1
                                    eval(['Fit_struct_run(iCell).True_run.s_',' = s;']);
                                elseif run == 2
                                    eval(['Fit_struct_run(iCell).True_norun.s_',' = s;']);
                                end
                            else
                                if run == 1
                                    if ~isnan(Fit_struct_run(iCell).True_run.s_.x(1))
                                        SAVEALLDATA = 0;
                                        PLOTIT_FIT = 0;
                                        Fit_2Dellipse_LG
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_',' = s;']);
                                    else
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.x',' = nan(1,6);']);
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_50', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_50', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_10', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_10', '= nan;']);
                                    end
                                elseif run == 2
                                    if ~isnan(Fit_struct_run(iCell).True_norun.s_.x(1))
                                        SAVEALLDATA = 0;
                                        PLOTIT_FIT = 0;
                                        Fit_2Dellipse_LG
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_',' = s;']);
                                    else
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.x',' = nan(1,6);']);
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_50', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_50', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_10', '= nan;']);
                                        eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_10', '= nan;']);
                                    end
                                end 
                            end
                        else
                            if count_shuf == 0
                                if run == 1
                                    eval(['Fit_struct_run(iCell).True_run.s_.x',' = nan(1,6);']);
                                    eval(['Fit_struct_run(iCell).True_run.s_.SFhicut_50', '= nan(1,1);']);
                                    eval(['Fit_struct_run(iCell).True_run.s_.TFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).True_run.s_.SFhicut_10', '= nan;']);
                                    eval(['Fit_struct_run(iCell).True_run.s_.TFhicut_10', '= nan;']);
                                elseif run == 2
                                    eval(['Fit_struct_run(iCell).True_norun.s_.x',' = nan(1,6);']);
                                    eval(['Fit_struct_run(iCell).True_norun.s_.SFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).True_norun.s_.TFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).True_norun.s_.SFhicut_10', '= nan;']);
                                    eval(['Fit_struct_run(iCell).True_norun.s_.TFhicut_10', '= nan;']);
                                end
                            elseif count_shuf> 0
                                if run == 1
                                    eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.x',' = nan(1,6);']);
                                    eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_10', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_10', '= nan;']);
                                elseif run == 2
                                    eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.x',' = nan(1,6);']);
                                    eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_50', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_10', '= nan;']);
                                    eval(['Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_10', '= nan;']);
                                end
                            end
                        end
                    end
                end
            end

            fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct_run.mat']);   
            save(fn_out, 'Fit_struct_run')
        
            if Nshuf>1;
                fit_true_vec_norun = [];
                fit_true_vec_run = [];
                fit_shuf_vec_norun = [];
                fit_shuf_vec_run = [];

                for iCell = 1:n_pix
                    for run = 1:2
                        if run == 1
                            if ~isempty(Fit_struct_run(iCell).True_run)                
                                eval(['tmp = Fit_struct_run(iCell).True_run.s_.x;']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_run.s_.SFhicut_50];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_run.s_.TFhicut_50];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_run.s_.SFhicut_10];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_run.s_.TFhicut_10];']);
                                fit_true_vec_run(iCell,:) = tmp;
                            end
                        elseif run == 2
                            if ~isempty(Fit_struct_run(iCell).True_norun)                
                                eval(['tmp = Fit_struct_run(iCell).True_norun.s_.x;']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_norun.s_.SFhicut_50];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_norun.s_.TFhicut_50];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_norun.s_.SFhicut_10];']);
                                eval(['tmp = [tmp Fit_struct_run(iCell).True_norun.s_.TFhicut_10];']);
                                fit_true_vec_norun(iCell,:) = tmp;
                            end
                        end
                    end
                end

                for count_shuf = 1:Nshuf
                    for iCell = 1:n_pix
                        for run = 1:2
                            if run == 1
                                if length(Fit_struct_run(iCell).Shuf_run)>= count_shuf
                                    if ~isempty(Fit_struct_run(iCell).Shuf_run(count_shuf).s_)
                                        eval(['tmp = Fit_struct_run(iCell).Shuf_run(count_shuf).s_.x;']);
                                        eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_50];']);
                                        eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_50];']);
                                        eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_run(count_shuf).s_.SFhicut_10];']);
                                        eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_run(count_shuf).s_.TFhicut_10];']);
                                        %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                                        fit_shuf_vec_run(iCell,:,count_shuf) = tmp;
                                        end
                                    end
                                elseif run == 2
                                    if length(Fit_struct_run(iCell).Shuf_norun)>= count_shuf
                                        if ~isempty(Fit_struct_run(iCell).Shuf_norun(count_shuf).s_)
                                            eval(['tmp = Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.x;']);
                                            eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_50];']);
                                            eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_50];']);
                                            eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.SFhicut_10];']);
                                            eval(['tmp = [tmp Fit_struct_run(iCell).Shuf_norun(count_shuf).s_.TFhicut_10];']);
                                            %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                                            fit_shuf_vec_norun(iCell,:,count_shuf) = tmp;
                                        end 
                                    end
                                end
                            end
                        end
                    end

                    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
                    load(fn_out)

                    Npars = size(fit_shuf_vec_run,2);
                    lbub_fits_new = zeros(n_pix,Npars,5,3);
                    lbub_fits_new(:,:,:,1) = lbub_fits;
                    alpha_bound = .025;
                    for iCell = 1:n_pix
                        for run = 1:2
                            for count2 = 1:Npars
                                if run == 1
                                    tmp = squeeze(fit_shuf_vec_run(iCell,count2,:));
                                elseif run ==2
                                    tmp = squeeze(fit_shuf_vec_norun(iCell,count2,:));
                                end
                                [i,j] = sort(tmp);
                                ind_shuf_lb = ceil(Nshuf*alpha_bound);
                                ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
                                lbub_fits_new(iCell,count2,1,run+1) = i(ind_shuf_lb);
                                lbub_fits_new(iCell,count2,2,run+1) = i(ind_shuf_ub);
                                lbub_fits_new(iCell,count2,3,run+1) = mean(i); 
                                lbub_fits_new(iCell,count2,5,run+1) = std(i);
                            end
                        %now take means from truedata fit:
                            if run == 1
                                lbub_fits_new(iCell,:,4, run+1) = fit_true_vec_run(iCell,:);
                            elseif run == 2
                                lbub_fits_new(iCell,:,4, run+1) = fit_true_vec_norun(iCell,:);
                            end
                        end
                    end

                lbub_fits = lbub_fits_new;

                lbub_diff = squeeze(lbub_fits(:,:,2,:)-lbub_fits(:,:,1,:));
                goodfit_ind = [];
                goodfit_ind.all = [];
                goodfit_ind.run = [];
                goodfit_ind.norun = [];
                for run = 2:3
                    for iCell = 1:n_pix
                        if lbub_diff(iCell,4,run)<2 
                            if lbub_diff(iCell,5,run)<2
                                if run == 1
                                    goodfit_ind.all = [goodfit_ind.all iCell];
                                elseif run == 2
                                    goodfit_ind.run = [goodfit_ind.run iCell];
                                elseif run == 3
                                    goodfit_ind.norun = [goodfit_ind.norun iCell];
                                end
                            end
                        end
                    end
                end

                fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits_run_norun.mat']);   
                save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit_ind')
            end
        end
    end
end
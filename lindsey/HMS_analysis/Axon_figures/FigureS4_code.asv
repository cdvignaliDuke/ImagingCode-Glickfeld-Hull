fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 'S4';
areas = strvcat('PM', 'LM', 'AL');
area_order = [2;3;1];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

%get F of boutons
for iArea = 1:3
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    F_all = [];
    for iexp = 1:nexp;   
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_stim = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        load(fn_stim);
        
        F_off = mean(stim_off,3);
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
        	if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
                suby = pos(1)-1:pos(1)+1;
                subx = pos(2)-1:pos(2)+1;
                F_bouton = mean(mean(F_off(suby, subx),2),1);
                F_all = [F_all F_bouton];
            end
        end
    end
    Goodfits(iArea).F = F_all;
end

%get variability of responses
for iArea = 1:3
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    resp_var = [];
    for iexp = 1:nexp;   
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        
        fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_resp);
        fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);
        
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
        
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1;
                [peak_dF,dF_ind] = max(max(all_fits(iArea).expt(iexp).bouton(iCell).plotfit,[],1),[],2);
                begin = sum(stim_reps(:,1:dF_ind-1))+1;
                stop = sum(stim_reps(:,1:dF_ind));
                resp_var = [resp_var; var(resp_dFoverF(iCell, begin:stop))];
            end
        end
    end
    Goodfits(iArea).var = resp_var;
end
            
for iArea = 1:3
    Goodfits(iArea).std = sqrt(Goodfits(iArea).var);
    Goodfits(iArea).cv = Goodfits(iArea).std./Goodfits(iArea).dF;
end

edges_dF = [0 .2 .4 .8 1.6 3];
plot_dF = [.1 .3 .6 1.2 2.3];
F_avg= zeros(5,3);
F_sem= zeros(5,3);
sigma_avg= zeros(5,3);
sigma_sem= zeros(5,3);
var_avg = zeros(5,3);
var_sem= zeros(5,3);
cv_avg = zeros(5,3);
cv_sem= zeros(5,3);

figure;
for iArea = 1:3
    [n bin] = histc(Goodfits(iArea).dF, edges_dF);
    for ibin = 1:5
        ind = find(bin==ibin);
        F_avg(ibin,iArea) = mean(Goodfits(iArea).F(ind));
        F_sem(ibin,iArea) = std(Goodfits(iArea).F(ind),[],2)./sqrt(length(ind));
        sigma_avg(ibin,iArea) = mean(mean([Goodfits(iArea).sigma_SF(ind) Goodfits(iArea).sigma_TF(ind)],2),1);
        sigma_sem(ibin,iArea) = std(mean([Goodfits(iArea).sigma_SF(ind) Goodfits(iArea).sigma_TF(ind)],2),[],1)./sqrt(length(ind));
        var_avg(ibin,iArea) = mean(Goodfits(iArea).var(ind));
        var_sem(ibin,iArea) = std(Goodfits(iArea).var(ind),[],1)./sqrt(length(ind));
        cv_avg(ibin,iArea) = mean(Goodfits(iArea).cv(ind));
        cv_sem(ibin,iArea) = std(Goodfits(iArea).cv(ind),[],1)./sqrt(length(ind));
    end
    subplot(2,2,1)
    errorbar(log2(plot_dF), F_avg(:,iArea), F_sem(:,iArea), col(iArea,:));
    hold on
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    axis square
    ylabel('bouton F');
    ylim([0 100])
    subplot(2,2,2)
    errorbar(log2(plot_dF), sigma_avg(:,iArea), sigma_sem(:,iArea), col(iArea,:));
    hold on
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    ylabel('sigma');
    ylim([0 3])
    axis square
    subplot(2,2,3)
    errorbar(log2(plot_dF), var_avg(:,iArea), var_sem(:,iArea), col(iArea,:));
    hold on
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    ylabel('variance');
    axis square
    subplot(2,2,4)
    errorbar(log2(plot_dF), cv_avg(:,iArea), cv_sem(:,iArea), col(iArea,:));
    hold on
    xlim([-3.5 1.5]);
    xlabel('dF/F');
    ylabel('cv');
    axis square
end

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_dFbinning_F_sigma_cv.ps']);
        print(gcf, '-depsc2', fn_out);

%stats on F, sigma and cv
dF_all = [];
F_all = [];
cv_all = [];
sigma_all = [];
for iArea = 1:3
    dF_all = [dF_all; Goodfits(iArea).dF];
    F_all = [F_all; Goodfits(iArea).F'];
    cv_all = [cv_all; Goodfits(iArea).cv];
    sigma_all = [sigma_all; mean([Goodfits(iArea).sigma_TF Goodfits(iArea).sigma_SF],2)];
end

[n_dF bin_dF] = histc(dF_all, edges_dF);
low_dF_ind = find(bin_dF<3);
high_dF_ind = find(bin_dF>3);

[h_F, p_F] = ttest2(F_all(low_dF_ind,:), F_all(high_dF_ind,:));
[h_cv, p_cv] = ttest2(cv_all(low_dF_ind,:), cv_all(high_dF_ind,:));
[h_sigma, p_sigma] = ttest2(sigma_all(low_dF_ind,:), sigma_all(high_dF_ind,:));

figure;
F_all_avg= zeros(5,2);
sigma_all_avg= zeros(5,2);
var_all_avg = zeros(5,2);
cv_all_avg = zeros(5,2);

[n_1 bin_1] = histc(Goodfits(1).dF, edges_dF);
[n_2 bin_2] = histc(Goodfits(2).dF, edges_dF);
[n_3 bin_3] = histc(Goodfits(3).dF, edges_dF);
for ibin = 1:5
    ind_1 = find(bin_1==ibin);
    ind_2 = find(bin_2==ibin);
    ind_3 = find(bin_3==ibin);
    ind_all = [ind_1; ind_2; ind_3];
    F_all_avg(ibin,1) = mean([Goodfits(1).F(ind_1) Goodfits(2).F(ind_2) Goodfits(3).F(ind_3)]);
    F_all_avg(ibin,2) = std([Goodfits(1).F(ind_1) Goodfits(2).F(ind_2) Goodfits(3).F(ind_3)],[],2)./sqrt(length(ind_all));
    sigma_all_avg(ibin,1) = mean([mean([Goodfits(1).sigma_SF(ind_1) Goodfits(1).sigma_TF(ind_1)],2); mean([Goodfits(2).sigma_SF(ind_2) Goodfits(2).sigma_TF(ind_2)],2); mean([Goodfits(3).sigma_SF(ind_3) Goodfits(3).sigma_TF(ind_3)],2)],1);
    sigma_all_avg(ibin,2) = std([mean([Goodfits(1).sigma_SF(ind_1) Goodfits(1).sigma_TF(ind_1)],2); mean([Goodfits(2).sigma_SF(ind_2) Goodfits(2).sigma_TF(ind_2)],2); mean([Goodfits(3).sigma_SF(ind_3) Goodfits(3).sigma_TF(ind_3)],2)],[],1)./sqrt(length(ind_all));
    var_all_avg(ibin,1) = mean([Goodfits(1).var(ind_1); Goodfits(2).var(ind_2); Goodfits(3).var(ind_3)]);
    var_all_avg(ibin,2) = std([Goodfits(1).var(ind_1); Goodfits(2).var(ind_2); Goodfits(3).var(ind_3)],[],1)./sqrt(length(ind_all));
    cv_all_avg(ibin,1) = mean([Goodfits(1).cv(ind_1); Goodfits(2).cv(ind_2); Goodfits(3).cv(ind_3)]);
    cv_all_avg(ibin,2) = std([Goodfits(1).cv(ind_1); Goodfits(2).cv(ind_2); Goodfits(3).cv(ind_3)],[],1)./sqrt(length(ind_all));
end
subplot(2,2,1)
errorbar(log2(plot_dF), F_all_avg(:,1), F_all_avg(:,2), '-k');
hold on
xlim([-3.5 1.5]);
xlabel('dF/F');
axis square
ylabel('bouton F');
ylim([0 100])
subplot(2,2,2)
errorbar(log2(plot_dF), sigma_all_avg(:,1), sigma_all_avg(:,2), '-k');
hold on
xlim([-3.5 1.5]);
xlabel('dF/F');
ylabel('sigma');
ylim([0 3])
axis square
subplot(2,2,3)
errorbar(log2(plot_dF), var_all_avg(:,1), var_all_avg(:,2),'-k');
hold on
xlim([-3.5 1.5]);
ylim([0 .4])
xlabel('dF/F');
ylabel('variance');
axis square
subplot(2,2,4)
errorbar(log2(plot_dF), cv_all_avg(:,1), cv_all_avg(:,2), '-k');
hold on
xlim([-3.5 1.5]);
ylim([0 1])
xlabel('dF/F');
ylabel('cv');
axis square

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_dFbinning_F_sigma_cv.ps']);
        print(gcf, '-depsc2', fn_out);
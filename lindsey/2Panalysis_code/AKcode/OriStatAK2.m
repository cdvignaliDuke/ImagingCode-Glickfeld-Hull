function OriStatKO = OriStatAK2(data_table,norm_table,alpha)
% Modification of Kenichi Ohki's ori stats. Changes by AK (compatibility w/ 16 directions, ori by TF analysis) 

%%%%%%%%%%%%%% arg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tc: timecourses of multiple ROIs. (Nframes) x (Ncells).
%       output of get_tcoursesCR can be used.
% Noff: Number of off frames
% Non: Number of on frames
% nstim_per_run: number of stimuli (usually gratings) per run
% run_ind: run indexes which you want to include in the table.
%           0-based (first run is 0).
% rest_length: length (frames) of rest period before the beginnings of stimuli.
% alpha: threshold for significance

%%%%%%%%%%%%% statKO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           data
%
% data_table: [10x9 double] Nrun x Nconditions (8 directions + 1 baseline). 
% dir_ratio_change: [1x8 double] ratio change of singal for each direction
%                   from the baseline
% ori_ratio_change: [1x4 double] ratio change of singal for each orientation
% from the baseline
%
%           anova statistics
%
% p_value_resp: p-value for responsiveness
% p_value_sel: p-value for selectivity
% dir_sig: significance of best response > null response. 
% if significant =1, else =0 (P < alpha).
%
%           classical analysis
%
% best_dir: best responding direction. 1st - 8th.
% null_dir: opposite direction to best_dir. 1st - 8th.
% R_best_dir: response for best direction (ratio change).
% R_null_dir: response for null direction (ratio change).
% R_min_dir: min response for 8 directions (ratio change).
% DI: classical direction index (1-R(null_dir)/R(best_dir)).
% best_ori: best responding orientation. 1st - 4th.
% null_ori: orthogonal to best orientation. 1st - 4th.
% R_best_ori: max response for orientation (ratio change).
% R_null_ori: response for null orientation (ratio change).
% R_min_ori: min response for orientation (ratio change).
%
%           vector average
%
% dir_vector_angle: preferred direction obtained from vector averaging
% dir_vector_mag: magnitude of averaged vector for direction
% dir_vector_tune: direction tuning index obtained from vector averaging 
% ori_vector_angle: preferred orientation obtained from vector averaging
% ori_vector_mag: magnitude of averaged vector for orientation
% ori_vector_tune: orientation tuning index obtained from vector averaging 
%
%           miscellaneous
%
%           miscellaneous from anova
%
% misc.anova_table_resp: anova table for responsiveness 
% misc.anova_stats_resp: anova statistics for responsiveness
% misc.sig_epochs_resp: [Nx2 double] significantly different pairs (P < alpha)
%                   from multiple comparisons (Tukey's t-test)
% misc.anova_table_sel: anova table for selectivity
% misc.anova_stats_sel: anova statistics for selectivity
% misc.sig_epochs_sel: [Nx2 double] ignificantly different pairs (P < alpha)
%                   from multiple comparisons (Tukey's t-test)
%
%
%   Kenichi Ohki 09/20/04   modified 09/22/04
%

printflag_for_fitting =0;

Ncells = size(norm_table,3);

for i=1:Ncells
    
    OriStatKO(i).data_table = squeeze(data_table(:,:,i));
    OriStatKO(i).norm_table = squeeze(norm_table(:,:,i));
    
    OriStatKO(i).dir_ratio_change = nanmean(OriStatKO(i).norm_table)-1; % ratio change
    OriStatKO(i).ori_ratio_change = nansum(reshape (OriStatKO(i).dir_ratio_change, size(norm_table,2)/2, 2)')/2;

    % anova
    [OriStatKO(i).p_value_resp,...
        OriStatKO(i).misc.anova_table_resp,...
        OriStatKO(i).misc.anova_stats_resp,...
        OriStatKO(i).misc.sig_epochs_resp] = anova1KO (OriStatKO(i).data_table, alpha);

    [OriStatKO(i).p_value_sel,...
        OriStatKO(i).misc.anova_table_sel,...
        OriStatKO(i).misc.anova_stats_sel,...
        OriStatKO(i).misc.sig_epochs_sel] = anova1KO (OriStatKO(i).norm_table, alpha);
    
    % classical analysis
    [OriStatKO(i).best_dir,...
        OriStatKO(i).null_dir,...
        OriStatKO(i).R_best_dir,...
        OriStatKO(i).R_null_dir,...
        OriStatKO(i).R_min_dir,...
        OriStatKO(i).DI]...
            = dir_indexKO(OriStatKO(i).dir_ratio_change);
    [OriStatKO(i).best_ori,...
        OriStatKO(i).null_ori,...
        OriStatKO(i).R_best_ori,...
        OriStatKO(i).R_null_ori,...
        OriStatKO(i).R_min_ori,...
        OriStatKO(i).OI]...
            = dir_indexKO(OriStatKO(i).ori_ratio_change);
    
    % vector average analysis

    [OriStatKO(i).dir_vector_angle,...
        OriStatKO(i).dir_vector_mag,...
        OriStatKO(i).dir_vector_tune]...
        = vector_average(OriStatKO(i).dir_ratio_change);

    [OriStatKO(i).ori_vector_angle,...
        OriStatKO(i).ori_vector_mag,...
        OriStatKO(i).ori_vector_tune]...
        = vector_average(OriStatKO(i).ori_ratio_change);
   
    OriStatKO(i).ori_vector_angle = OriStatKO(i).ori_vector_angle /2;
    
        % swindale's tuning curve fitting
    if OriStatKO(i).p_value_sel < alpha
        [OriStatKO(i).best_dir_fit,...
                OriStatKO(i).null_dir_fit,...
                OriStatKO(i).R_best_dir_fit,...
                OriStatKO(i).R_null_dir_fit,...
                OriStatKO(i).DI_fit,...
                OriStatKO(i).DS,...
                OriStatKO(i).dir_tuning_width,...
                OriStatKO(i).dir_A1,...
                OriStatKO(i).dir_A2,...
                OriStatKO(i).misc.dir_k1,...
                OriStatKO(i).misc.dir_k2,... 
                OriStatKO(i).misc.dir_phi2,...
                OriStatKO(i).misc.dir_A,...
                OriStatKO(i).misc.dir_resnorm]...
                = estimate_dir_tuningKO (OriStatKO(i).dir_ratio_change,1,'spline',printflag_for_fitting);
        [OriStatKO(i).best_ori_fit,...
                OriStatKO(i).ori_tuning_width,...
                OriStatKO(i).ori_A1,...
                OriStatKO(i).misc.ori_k1,...
                OriStatKO(i).misc.ori_A,...
                OriStatKO(i).misc.ori_resnorm]...
                = estimate_ori_tuningKO (OriStatKO(i).ori_ratio_change, 1, 'spline',printflag_for_fitting);
    else
        OriStatKO(i).best_dir_fit =NaN;
        OriStatKO(i).null_dir_fit =NaN;
        OriStatKO(i).R_best_dir_fit =NaN;
        OriStatKO(i).R_null_dir_fit =NaN;
        OriStatKO(i).DI_fit =NaN;
        OriStatKO(i).DS =NaN;
        OriStatKO(i).dir_tuning_width =NaN;
        OriStatKO(i).dir_A1 =NaN;
        OriStatKO(i).dir_A2 =NaN;
        OriStatKO(i).misc.dir_k1 =NaN;
        OriStatKO(i).misc.dir_k2 =NaN; 
        OriStatKO(i).misc.dir_phi2 =NaN;
        OriStatKO(i).misc.dir_A =NaN;
        OriStatKO(i).misc.dir_resnorm =NaN;
        OriStatKO(i).best_ori_fit =NaN;
        OriStatKO(i).ori_tuning_width =NaN;
        OriStatKO(i).ori_A1 =NaN;
        OriStatKO(i).misc.ori_k1 =NaN;
        OriStatKO(i).misc.ori_A =NaN;
        OriStatKO(i).misc.ori_resnorm = NaN;
    end
  
    if find_pair_in_list(OriStatKO(i).misc.sig_epochs_sel, [OriStatKO(i).best_dir, OriStatKO(i).null_dir]);
        OriStatKO(i).dir_sig =1;
    else
        OriStatKO(i).dir_sig =0;
    end
    
end

OriStatKO = ori_CI_upton (OriStatKO,0.05);
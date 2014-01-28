function doOriStatKO(s,iter,total_iters,varargin)
% run Kenichi's ori stats.  New version by CR.
% saves an array of structures rather than a cell arry.

% reviesed June 30 2005 by KO
% p_values for orientation selectivity are included.

if s.Nstim_per__run  ~= 8
    return
end

rest_length =1;
alpha = 0.01;
printflag_for_interp=0;
printflag_for_fitting=0;
run_ind = [0:s.Nrep-1];

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
%           clay's dot product
%
% best_dir_dp: better responding direction orthogonal to best orientation,
%               obtained from Clay's dot product 
% null_dir_dp: less responding direction orthogonal to best orientation,
%               obtained from Clay's dot product 
% R_best_dir_dp: response to best_dir_dp obtained from Clay's dot product 
% R_null_dir_dp: response to null_dir_dp obtained from Clay's dot product
% DI_dp: Clay's direction index from dot product
%       (1-R_null_dir_dp / R_best_dir_dp)
%
%           kenichi's interpolation
%
% best_dir_interp: best direction from interpolation
% null_dir_interp: opposite direction to best_dir_interp
% R_best_dir_interp: response to best_dir_interp,
%                           obtained from interpolation
% R_null_dir_interp: response to null_dir_interp,
%                           obtained from interpolation
% DI_interp: direction index from interpolation
%               1-R(null_dir_interp)/R(best_dir_interp)          
%
%           swindale's gaussian fitting
%
% best_dir_fit: preferred direction obtained from curve fitting.
% null_dir_fit: oppposite direction to best_dir_fit
% R_best_dir_fit: response to best direction
% R_null_dir_fit: response to null direction
% DI_fit: direction index from curve fitting
%               (1-R_best_dir_fit/R_null_dir_fit).
% DS: Swindale's direction selectivity index.
% dir_tuning_width: direction tuning width from curve fitting.
% dir_A1: amplitude of the larger peak for direction (ratio change).
% dir_A2: amplitude of the smaller peak for direction (ratio change).
%
% best_ori_fit: preferred orientation obtained form curve fitting. 
% ori_tuning_width: orientation tuning width from curve fitting. 
% ori_A1: amplitude of the peak for orientation (ratio change).
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
%           miscellaneous from swindale's fitting
%
% misc.dir_k1: inverse (not exact) of direction tuning width for the larger peak. 
% misc.dir_k2: inverse (not exact) of direction tuning width for the smaller peak.
% misc.dir_phi2: direction of the second peak.
% misc.dir_A: all the fitting parameters
% misc.dir_resnorm: residual of curve fitting for direction.
% misc.ori_k1: inverse (not exact) of orientation tuning width.
% misc.ori_A: all the fitting parameters
% misc.ori_resnorm: residual of curve fitting for orientation.
%
%   Kenichi Ohki 09/20/04   modified 09/22/04
%

  [procoutfname, matoutfname]= getoutfilenameNewCR(s); % standard calling with one matfile and one procfile  
  load ([procoutfname,'_tcourse.mat'],'tc');   % THIS SHOULD HAVE ARRAY 'avg_img'
  % NOW THE REAL WORK BEGINS
  whos

nframes_per_run = (s.Noff+s.Non)*s.Nstim_per__run;
[Nframes, Ncells] = size(tc);

% selecting default parameters for epochs

epochs = select_epochs (s.Noff, s.Non, s.Nstim_per__run, rest_length);

%CR DELETED: should have structure array, not cellarray OriStatKO = cell(Ncells,1);
% changed all curly braces to parens.
for i=1:Ncells
    i
    % low cut filtering
    tc_lowcut = low_cutKO(tc(:,i), nframes_per_run*2, 'gaussian', 1);
    
    OriStatKO(i).data_table = average_epoch (tc_lowcut, epochs, nframes_per_run, run_ind);

    temp=sum(OriStatKO(i).data_table);
    OriStatKO(i).dir_ratio_change = temp(1:s.Nstim_per__run)/temp(s.Nstim_per__run+1)-1; % ratio change
    OriStatKO(i).ori_ratio_change = sum(reshape (OriStatKO(i).dir_ratio_change, s.Nstim_per__run/2, 2)')/2;

    % anova
    [OriStatKO(i).p_value_resp,...
        OriStatKO(i).misc.anova_table_resp,...
        OriStatKO(i).misc.anova_stats_resp,...
        OriStatKO(i).misc.sig_epochs_resp] = anova1KO (OriStatKO(i).data_table, alpha);

    [OriStatKO(i).p_value_sel,...
        OriStatKO(i).misc.anova_table_sel,...
        OriStatKO(i).misc.anova_stats_sel,...
        OriStatKO(i).misc.sig_epochs_sel] = anova1KO (OriStatKO(i).data_table(:,1:s.Nstim_per__run,:), alpha);

    temp = (OriStatKO(i).data_table(:,1:s.Nstim_per__run/2) + OriStatKO(i).data_table(:,s.Nstim_per__run/2+1:s.Nstim_per__run))/2;
    temp=[temp,OriStatKO(i).data_table(:,s.Nstim_per__run+1)];
    
    [OriStatKO(i).p_value_ori_resp,...
        OriStatKO(i).misc.anova_table_ori_resp,...
        OriStatKO(i).misc.anova_stats_ori_resp,...
        OriStatKO(i).misc.sig_epochs_ori_resp] = anova1KO (temp, alpha);
    [OriStatKO(i).p_value_ori_sel,...
        OriStatKO(i).misc.anova_table_ori_sel,...
        OriStatKO(i).misc.anova_stats_ori_sel,...
        OriStatKO(i).misc.sig_epochs_ori_sel]  = anova1KO (temp(:,1:s.Nstim_per__run/2), alpha);

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

    % clay's dot product analysis
    
    [OriStatKO(i).best_dir_dp,...
        OriStatKO(i).null_dir_dp,...
        OriStatKO(i).R_best_dir_dp,... 
        OriStatKO(i).R_null_dir_dp,...
        OriStatKO(i).DI_dp]...
        = clayDI (OriStatKO(i).dir_ratio_change, OriStatKO(i).ori_vector_angle);
    
    % kenichi's interpolate analysis
    [OriStatKO(i).best_dir_interp,...
        OriStatKO(i).null_dir_interp,...
        OriStatKO(i).R_best_dir_interp,...
        OriStatKO(i).R_null_dir_interp,...
        OriStatKO(i).DI_interp] = interp_dir_tuningKO (OriStatKO(i).dir_ratio_change,'fourier',printflag_for_interp);

   
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

close all; 

oristatfname = [procoutfname, '_OriStat.mat']

save (oristatfname, 'OriStatKO');
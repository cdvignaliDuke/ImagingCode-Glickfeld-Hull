function stats = tcStats(data, alpha)

%%%%%%%%%%%%%%%%% minimal statistics %%%%%%%%%%%%%%%%%%
% this gets only p_values and best response etc.
% 04/01/08 Kenichi Ohki

if nargin < 2 
    alpha = 0.05;
end

ncells=length(data);

for icell=1:ncells
    disp(icell)
    
    this = data{icell};
        
    [ntrials,nstimplus]=size(this);
    
    nstim = nstimplus - 1;

    temp=sum(this,1);
    
    stat.data_table=this;
    stat.dir_ratio_change = temp(1:nstim)/temp(nstim+1)-1; % ratio change

    if ntrials <=1
        stats(icell)=stat;
        continue;
    end

    
    % anova
    [stat.p_value_resp,...
        stat.misc.anova_table_resp,...
        stat.misc.anova_stats_resp,...
        stat.misc.sig_epochs_resp] = anova1mult (this, alpha);

    [stat.p_value_sel,...
        stat.misc.anova_table_sel,...
        stat.misc.anova_stats_sel,...
        stat.misc.sig_epochs_sel] = anova1mult (this(:,1:nstim,:), alpha);

    % best and worst response
    [stat.best_dir,...
        stat.min_dir,...
        stat.R_best_dir,...
        stat.R_min_dir,...
        stat.sel_index]...
            = find_best_stim(stat.dir_ratio_change);
        
    stats(icell)=stat;
end

return;
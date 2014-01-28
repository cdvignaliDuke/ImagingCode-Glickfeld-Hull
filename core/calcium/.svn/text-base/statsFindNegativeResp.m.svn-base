function negative_cells=statsFindNegativeResp (stats, alpha)

% find negatively responsive cells with a statistical significance 
% positvely responsive cells are defined when the following two conditions are satisfied:
% 1) p<alpha, ANOVA, across stims + baseline
% 2) at least one stim < baseline, p<alpha, post-hoc Tukey's t-test
% 2010.01.23    Kenichi Ohki

if nargin <2
    alpha=0.05;
end

Ncells=length(stats);
negative_cells=[];

for i=1:Ncells
    data=stats(i).data_table;
    [Nrep,Nstim]=size(data);
    baseline=Nstim;

    % anova
    [p_value_resp,...
        anova_table_resp,...
        anova_stats_resp,...
        sig_epochs_resp] = anova1mult (data, alpha);


  %  if p_value_resp >=alpha | isempty(sig_epochs_resp)
    if isempty(sig_epochs_resp)
        continue;
    else
        if ~isempty(find(sig_epochs_resp(:,1)==baseline))
            negative_cells=[negative_cells,i];
        end
    end
end
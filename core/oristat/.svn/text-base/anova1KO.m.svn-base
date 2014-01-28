function [p_value, anova_table, anova_stats, sig_epochs] = anova1KO (table, alpha);

% extenstion of anova1 in matlab.
% returns not only ANOVA results, but also multlcompare resutls.
%
% table: data table. (Number of runs) x(Number of epochs).
%           each cell contains an average value of timecourse
%           during the epoch.
% alpha: significance threshold for multiple comparisons.
%
% p_value: p_value of ANOVA for each cell. 
% anova_table: ANOVA tables for each cell. 
% anova_stats: ANOVA statistics for each cell. 
% sig_epochs: signigicantly different epoch pairs 
%               obtained by multiple comparisons. 
%               N x 2 matrix.
%               each row of matrix represents a pair of epochs,
%               which satisfies (i,1) > (i,2) (P < alpha)
%
%   Kenichi Ohki 09/16/04
%

[Nruns, Nepochs] = size(table);

[p_value, anova_table, anova_stats] =anova1(table,[],'off');
    
sig=multcompare(anova_stats, alpha,'off');
sig_epochs=[];
for j=1:size(sig,1)
    if (sig(j,3)*sig(j,5) >0)
        if (sig(j,3)>0)
            sig_epochs=[sig_epochs; [sig(j,1),sig(j,2)]];
        else
            sig_epochs=[sig_epochs; [sig(j,2),sig(j,1)]];
        end
    end
end



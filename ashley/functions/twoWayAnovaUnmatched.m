function anovaStruct = twoWayAnovaUnmatched(data,varNames)
% data should be a M xN cell array with M rows corresponding to grouping 1,
% and N columns corresponding to grouping 2
% Interaction is included
    [nrow,ncol] = size(data);
    n = cellfun(@length,data);
    grpIDcell = cellfun(@(x) ones(size(x)),data,'unif',0);
    rowGrpID = cell(size(grpIDcell));
    colGrpID = cell(size(grpIDcell));
    for irow = 1:nrow
        rowGrpID(irow,:) = cellfun(@(x) x.*irow,grpIDcell(irow,:),'unif',0);
    end
    for icol = 1:ncol
        colGrpID(:,icol) = cellfun(@(x) x.*icol,grpIDcell(:,icol),'unif',0);
    end
    data_t = cell2mat(data);
    rowID_t = cell2mat(rowGrpID);
    colID_t = cell2mat(colGrpID);
    [p,~,stats] = anovan(data_t(:)',{rowID_t(:)',colID_t(:)'},...
        'display','off','model','interaction');
%     [p,~,stats] = anovan(data_t(:)',{rowID_t(:)',colID_t(:)'},...
%         'display','off');
    posthoc = cell(2,1);
    for ip = 1:2
        if p(ip) < 0.05
            anovaTable = multcompare(stats,'dimension',ip);
            posthoc{ip} = anovaTable(:,[1:2,6]);
        end
    end

    anovaStruct.varTestName = cat(1,varNames,{'Int'});
    anovaStruct.p = p;
    anovaStruct.posthoc = posthoc;
end
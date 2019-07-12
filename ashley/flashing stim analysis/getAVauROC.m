function [avauroc, avauroc_test] = getAVauROC(responses1,responses2)
% responses should be trials x cells

ncells = size(responses1,2);
avauroc = nan(1,ncells);
avauroc_test = nan(1,ncells);
for icell = 1:ncells
    avauroc(icell) = roc_gh(responses1(:,icell),responses2(:,icell));
    avauroc_test(icell) = ranksum(responses1(:,icell),responses2(:,icell));
end

end

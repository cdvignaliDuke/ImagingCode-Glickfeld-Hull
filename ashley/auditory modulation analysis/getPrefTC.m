function prefTC = getPrefTC(stimSortedTC,pref_ind)
% stimSortedTC is a conditions by n stim cell array of timecourses for each
% cell 
% pref_ind is the preferred stimulus index

ncondition = size(stimSortedTC,2);
ncell = size(stimSortedTC{1,1},2);

prefTC = cell(1,ncondition);
for icondition = 1:ncondition
    for icell = 1:ncell
        p = pref_ind(icell);
        d = stimSortedTC{p,icondition}(:,icell);

        prefTC{icondition}(:,icell) = d;    
    end
end

end
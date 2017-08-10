function orthTC = getOrthTC(stimSortedTC,pref_ind)
% stimSortedTC is a conditions by n stim cell array of timecourses for each
% cell 
% pref_ind is the preferred stimulus index
nstim = max(pref_ind);
pref_shift = floor(nstim/2);
orth_ind = pref_ind+pref_shift;
orth_ind(orth_ind > nstim) = pref_ind(orth_ind > nstim) - pref_shift;

ncondition = size(stimSortedTC,2);
ncell = size(stimSortedTC{1,1},2);

orthTC = cell(1,ncondition);
for icondition = 1:ncondition
    for icell = 1:ncell
        p = orth_ind(icell);
        d = stimSortedTC{p,icondition}(:,icell);

        orthTC{icondition}(:,icell) = d;    
    end
end

end
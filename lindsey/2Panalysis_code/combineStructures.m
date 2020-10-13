function outS = combineStructures(expt1,expt2);
%merges fields from two structures into single structure with matched
%fields. discards fields that are not matched.
    
fNames1 = fieldnames(expt1);
fNames2 = fieldnames(expt2);
fNames_comb = intersect(fNames1,fNames2);
nF = length(fNames_comb);
nexp1 = size(expt1,2);
nexp2 = size(expt2,2);
for iF = 1:nF
	fN = cell2mat(fNames_comb(iF,:));
    for i = 1:nexp1
        outS(i).(fN) = expt1(i).(fN);
    end
    for i = 1:nexp2
         outS(i+nexp1).(fN) = expt2(i).(fN);
    end
end
    
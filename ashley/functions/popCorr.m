function pCorr = popCorr(dff,frRateHz,timebinS);
nc = size(dff,2);

nFr = size(dff,1);
nFrPerBin = round(timebinS*frRateHz);

checkBinSz = mod(nFr,nFrPerBin);

if checkBinSz > 0
   nFr = nFr-checkBinSz;
   dff = dff(1:nFr,:);
end

dff_bin = reshape(dff,nFrPerBin,nFr/nFrPerBin,nc);
dff_down = squeeze(sum(dff_bin,1));

pCorr = nan(1,nc);
for i = 1:nc
    j = logical(ones(1,nc));
    j(i) = false;

    fi = dff_down(:,i);
    fj = dff_down(:,j);
    pCorr(i) = corr(fi,sum(fj,2));
end
end
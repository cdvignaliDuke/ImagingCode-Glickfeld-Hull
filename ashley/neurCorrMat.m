function [nCorrMat,nFRMat] = neurCorrMat(dff,frRateHz,timebinS)
nc = size(dff,2);

nFr = size(dff,1);
nFrPerBin = round(timebinS*frRateHz);

checkBinSz = mod(nFr,nFrPerBin);

if checkBinSz > 0
   nFr = nFr-checkBinSz;
   dff = dff(1:nFr,:);
end

% dff_bin = reshape(dff,nFrPerBin,nFr/nFrPerBin,nc);
% dff_down = squeeze(sum(dff_bin,1));

% nCorrMat = corr(dff_down);
nCorrMat = corr(dff);

nFRMat_prod = nan(nc);
for i = 1:nc
    for j = 1:nc
        FR1 = mean(dff(:,i));
        FR2 = mean(dff(:,j));
        nFRMat_prod(i,j) = FR1.*FR2;
    end
end
nFRMat_prod_pos = nFRMat_prod;
nFRMat_prod_pos(nFRMat_prod_pos<0) = NaN;

nFRMat = sqrt(nFRMat_prod_pos);

end
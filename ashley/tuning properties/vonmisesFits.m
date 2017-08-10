function [cellFits,R_square] = vonmisesFits(data,dataBoot,tuningReliabilityThresh,directions)
% output orientation tuning matrix for each cell from data (n oris+1 x
% cells) and resampled data (n oris+1 x cells x number of bootstraps)
% cells that aren't reliably fit (90% confidence) are left as NaN.
nboot = size(dataBoot,3);
data4fit = cat(1,data,data(1,:))';
dataBoot4fit = permute(cat(1,dataBoot,dataBoot(1,:,:)),[2,1,3]);

if any(directions > 180)
    orisFit = [directions(1:4) 180];
    orisSmooth = orisFit(1):1:orisFit(end);
else
    orisFit = [directions 180];
    orisSmooth = orisFit(1):1:orisFit(end);    
end

[~, theta_dist, theta_90] = vonmisesReliableFit(data4fit, dataBoot4fit,orisFit,nboot);

nc = size(data4fit,1);
R_square = nan(1,nc);
cellFits = nan(length(orisSmooth),nc);
for icell = 1:nc
    if theta_90(icell) < tuningReliabilityThresh
        [b, k, R, u, ~,R_square(icell)] = miaovonmisesfit_ori(deg2rad(orisFit),data4fit(icell,:));
        cellFits(:,icell) = b+R.*exp(k.*(cos(2.*(deg2rad(orisSmooth)-u))-1));
    end
end

end


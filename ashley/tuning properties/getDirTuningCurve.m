function [tuningCurve,R_square] = getDirTuningCurve(theta,data)

theta_smooth = 0:1:360;
nCells = size(data,1);

y_fit = nan(length(theta_smooth),nCells);
R_square = nan(nCells,1);
for icell = 1:nCells
    [b, k, R1, R2,u1, u2,~,R_square(icell)] = ...
        miaovonmisesfit_dir(deg2rad(theta),data(icell,:));
    y_fit(:,icell) = b+R1.*exp(k.*(cos((deg2rad(theta_smooth)-u1))-1)) + ...
        b+R2.*exp(k.*(cos((deg2rad(theta_smooth)-u2))-1));
end
tuningCurve = y_fit;
end


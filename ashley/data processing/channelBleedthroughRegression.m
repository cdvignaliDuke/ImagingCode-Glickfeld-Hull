
function [redImg_corrected,s_star,si_regress,Bi_regress,Bi_star] = greenBleedIntoRed(doRegressOnly,redImg,data_r,data_g,mask)
    if size(data_r,3) > 20000
        data_r = data_r(:,:,1:20000);
        data_g = data_g(:,:,1:20000);
    end
%     redImg = mean(data_r,3);
    if nargin == 5
        red_mask_cell = mask;
    else
        red_mask_cell = maskFromMultiMaxDFFStack(redImg);
    end
    %%
    ri = cellfun(@(x) x(:),stackGetPixVals(data_r,red_mask_cell),'unif',0);
    gi = cellfun(@(x) x(:),stackGetPixVals(data_g,red_mask_cell),'unif',0);
    nc = length(ri);
    % get the slope and offset for the regression model for each roi, then find
    % mean slope
    Bi_regress = nan(nc,1);
    si_regress = nan(nc,1);
    for i = 1:nc
        lr = fitlm(gi{i},ri{i});
        Bi_regress(i) = lr.Coefficients.Estimate(1);
        si_regress(i) = lr.Coefficients.Estimate(2);
    end

    % compute the cost value for the regression model and get optimized s and B
    if doRegressOnly == 1
        s_star = nan;
        Bi_star = nan;
        redImg_corrected = redImg - (mean(data_g,3).*mean(si_regress));
    else
        [s_star,Bi_star] = regressCommonSlopeModel(gi,ri);
        redImg_corrected = redImg - (mean(data_g,3).*s_star);
    end

end
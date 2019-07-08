
function [redImg_corrected,s_star,s_regress] = greenBleedIntoRed(redImg,data_r,data_g)
    if size(data_r,3) > 20000
        data_r = data_r(:,:,1:20000);
        data_g = data_g(:,:,1:20000);
    end
%     redImg = mean(data_r,3);
    red_mask_cell = maskFromMultiMaxDFFStack(redImg);
    %%
    ri = cellfun(@(x) x(:),stackGetPixVals(data_r,red_mask_cell),'unif',0);
    gi = cellfun(@(x) x(:),stackGetPixVals(data_g,red_mask_cell),'unif',0);
    nc = length(ri);
    % get the slope and offset for the regression model for each roi, then find
    % mean slope
    Bi = nan(nc,1);
    si = nan(nc,1);
    for i = 1:nc
        lr = fitlm(gi{i},ri{i});
        Bi(i) = lr.Coefficients.Estimate(1);
        si(i) = lr.Coefficients.Estimate(2);
    end
    s_regress = mean(si);

    % compute the cost value for the regression model and get optimized s and B

    [s_star,Bi_star] = regressCommonSlopeModel(gi,ri);

    redImg_corrected = redImg - (mean(data_g,3).*s_star);
end
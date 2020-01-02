
function [img1_corrected,s_star,si_regress,Bi_regress,Bi_star] = ...
channelBleedRegression(doRegressOnly,maskImg,data_1,data_2,mask)
% calculate the fraction of bleedthough of data_2 into data_1
    if size(data_1,3) > 20000
        data_1 = data_1(:,:,1:20000);
        data_2 = data_2(:,:,1:20000);
    end
%     redImg = mean(data_r,3);
    if nargin == 5
        mask_cell = mask;
    else
        mask_cell = maskFromMultiMaxDFFStack(maskImg);
    end
    %%
    i1 = cellfun(@(x) x(:),stackGetPixVals(data_1,mask_cell),'unif',0);
    i2 = cellfun(@(x) x(:),stackGetPixVals(data_2,mask_cell),'unif',0);
    nc = length(i1);
    % get the slope and offset for the regression model for each roi, then find
    % mean slope
    Bi_regress = nan(nc,1);
    si_regress = nan(nc,1);
    for i = 1:nc
        lr = fitlm(i2{i},i1{i});
        Bi_regress(i) = lr.Coefficients.Estimate(1);
        si_regress(i) = lr.Coefficients.Estimate(2);
    end

    % compute the cost value for the regression model and get optimized s and B
    if doRegressOnly == 1
        s_star = nan;
        Bi_star = nan;
        img1_corrected = maskImg - (mean(data_2,3).*mean(si_regress));
    else
        [s_star,Bi_star] = regressCommonSlopeModel(i2,i1);
        img1_corrected = maskImg - (mean(data_2,3).*s_star);
    end

end
function img = getPixelCorrelationImage(data)
    [ypix,xpix,nfr] = size(data);
    corrImage = zeros(ypix,xpix);
    for iX = 1:xpix
        for iY = 1:ypix
            P0_TC = squeeze(data(iY,iX,:));

            PK_TC = nan(8,nfr);

            if iY > 1 & iY < ypix & iX > 1 & iX < xpix
                PK_TC(1:3,:) = squeeze(data(iY-1,iX-1:iX+1,:));
                PK_TC(4,:) = squeeze(data(iY,iX-1,:));
                PK_TC(5,:) = squeeze(data(iY,iX+1,:));
                PK_TC(6:8,:) = squeeze(data(iY+1,iX-1:iX+1,:));
            end
            PK_TC = PK_TC';

            P0_corr = zeros(1,8);
            for iP = 1:8
                P0_corr(iP) = corr(PK_TC(:,iP),P0_TC);
            end
            meanP0Corr = mean(P0_corr);
            if meanP0Corr == 0 | isnan(meanP0Corr)
                meanP0Corr = 0;
            end
            corrImage(iY,iX) = mean(P0_corr);
        end
    end
    img = corrImage;
end

function respImageFig = plotImageStackAsSubplots(imageStack,titles)
    nplots = size(imageStack,3);
    [subplotRows, subplotColumns] = optimizeSubplotDim(nplots);
    respImageFig = figure; colormap gray
    for iplot = 1:nplots
       subplot(subplotRows,subplotColumns,iplot)
       imagesc(imageStack(:,:,iplot))
       title(titles(iplot))
    end
end
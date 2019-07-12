function regImgStartFrame = compareRegImg_2color(rc,ds,SubNum,expDate,data_g,data_r)

nfr = size(data_g,3);
nMeanImages = 9;
[nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

randStarterFrames = sort(randsample(nfr-100,nMeanImages));
setFigParams4Print('landscape')
figure
suptitle([SubNum '-' expDate '-Green'])
for iimg = 1:nMeanImages
    subplot(nRows,nCols,iimg)
    imagesc(mean(data_g(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+99),3))
    title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
end
if ~exist(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds),'dir')
    mkdir(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds))
end
print(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    ['rand image samples_' SubNum '-' expDate]),'-dpdf','-fillpage')
savefig(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
    ['rand image samples_' SubNum '-' expDate]))
if nargin > 5
    figure
    suptitle([SubNum '-' expDate '-Red'])
    for iimg = 1:nMeanImages
        subplot(nRows,nCols,iimg)
        imagesc(mean(data_r(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+99),3))
        title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
    end
    print(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
        ['red_rand image samples_' SubNum '-' expDate]),'-dpdf','-fillpage')
    savefig(fullfile(rc.ashleyAnalysis,'FSAV Summaries',ds,...
        ['red_rand image samples_' SubNum '-' expDate]))
end

regImgStartFrame = input('Enter Registration Image Start Frame, ENTER INTO DS:');
end
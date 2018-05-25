clear all
close all
ds = 'adaptDREADDS_V1_SOM_control';
rc = behavConstsAV;
eval(ds)
slct_expt = 5:6;
%%
for iexp = slct_expt

mouse = expt(iexp).mouse;
subnum = mouse;
expDate = expt(iexp).date;
    
adaptFolder = expt(iexp).adapt{1};
fName = [adaptFolder '_000_000'];

data = loadsbx_choosepmt(1,mouse,expDate,adaptFolder,fName);
nfr = size(data,3); 

%%
fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

nMeanImages = 9;
[nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

randStarterFrames = sort(randsample(nfr-100,nMeanImages));
setFigParams4Print('landscape')
figure
suptitle([subnum '-' expDate])
for iimg = 1:nMeanImages
    subplot(nRows,nCols,iimg)
    imagesc(mean(data(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+100),3))
    title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
end
print(['Z:\Analysis\Expt Summaries\' ds '\rand image samples_' subnum '-' expDate],'-dpdf','-fillpage')
savefig(['Z:\Analysis\Expt Summaries\' ds '\rand image samples_' subnum '-' expDate])

clear data
end

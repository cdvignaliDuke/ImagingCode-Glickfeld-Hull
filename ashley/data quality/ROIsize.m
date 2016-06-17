datasetStr = '_AL';

if isempty(datasetStr)
    dataGroup = 'awFSAVdatasets';
else
    dataGroup = ['awFSAVdatasets' datasetStr];
end
eval(dataGroup)

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

% iexp = 1;
for iexp = 1:size(expt,2)

rc = behavConstsAV;
try
load(fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging',expt(iexp).date,expt(iexp).dirtuning,'mask&TCDir.mat'))
try
    maxDFoverF = readtiff(fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging',expt(iexp).date,expt(iexp).dirtuning,'maxDFoverF.tif'));
catch
    maxDFoverF = readtiff(fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging',expt(iexp).date,expt(iexp).dirtuning,'maxDF.tif'));
end

mask_cell_tbl = floor(tabulate(mask_cell(mask_cell~=0)));

ROIhist = figure;
hist(mask_cell_tbl(:,2),100)
xlabel('ROI size')
ylabel('n cells')
axis square

sz_edges = [0 80 max(mask_cell_tbl(:,2))+1];
[h sz_bin] = histc(mask_cell_tbl(:,2),sz_edges);
bins = unique(sz_bin);

cell_sz_ind = cell(length(bins),1);
for ibin = 1:length(bins)    
    cell_sz_ind{ibin} = mask_cell_tbl(find(sz_bin == bins(ibin)),1);
end

mask_dim = size(mask_cell);

mask_binned = zeros([mask_dim 3]);
for ibin = 1:length(bins)
   ind = ismember(mask_cell(:),cell_sz_ind{ibin});
   temp = mask_cell(:);
   temp(~ind) = 0;
   mask_binned(:,:,ibin) = reshape(temp,mask_dim);
end

ROIimg = figure;
for iplot = 1:length(bins)
    subplot(length(bins)+1,1,iplot)
    imagesc(mask_binned(:,:,iplot))
    title(['#pix < ' num2str(sz_edges(iplot+1))])
end
subplot(length(bins)+1,1,length(bins)+1)
imagesc(maxDFoverF)
title('max dF/F')
colormap gray

figure(ROIhist)
print([fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging',expt(iexp).date) 'ROIhistogram.pdf'], '-dpdf')
figure(ROIimg)
print([fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging',expt(iexp).date) 'ROIimg.pdf'], '-dpdf')
catch
    disp(['skipped ' expt(iexp).mouse '-' expt(iexp).date])
end

end
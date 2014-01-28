%extract time courses
stack = readtiff(fn);

ind = strfind(fn, num2str(userun))+length(num2str(userun));
name = fn(ind:end-4);

rois = zeros(npix*pts,n/pts,size(stack,3));
for iroi = 1:n/pts;
    for ipix = 1:36;
        rois(ipix,iroi,:) = stack(roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2),:);
    end
end

fov_tc = squeeze(mean(mean(stack,1),2));
fov_mean = squeeze(mean(fov_tc,1));
roi_avg = bsxfun(@minus, squeeze(mean(rois,1))', fov_tc)+fov_mean;
roi_F = mean(roi_avg,1);
roi_dF = bsxfun(@minus, roi_avg, roi_F);

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) name '_roi_TCs.mat']);
save(fn_out, 'roi_avg', 'roi_dF');

for iroi1 = 1:n/4;
    for iroi2 = 1:n/4;
        r_corr(iroi1, iroi2) = triu2vec(corrcoef(roi_dF(:,iroi1),roi_dF(:,iroi2)));
    end
end
figure;
imagesq(r_corr);
colormap hot;
colorbar;
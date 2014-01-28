rois = zeros(npix*pts,n/pts,size(stack,3));
for iroi = 1:n/pts;
    for ipix = 1:36;
        rois(ipix,iroi,:) = stack(roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2),:);
    end
end

fov_tc = squeeze(mean(mean(stack,1),2));
fov_mean = squeeze(mean(fov_tc,1));
roi_avg = bsxfun(@minus, squeeze(mean(rois,1))', fov_tc)+fov_mean;

b = 5;
nRoi = size(roi_avg,2);
r = zeros(size(stack,1), size(stack,2), nRoi);

for iroi = 1:nRoi;
    for iy = b+1:240-b
        fprintf('.');
        for ix = b+1:256-b
            sub = stack(iy-3:iy+3,ix-3:ix+3,:);
            sub_avg = mean(mean(sub,2),1);
            r(iy,ix,iroi)= triu2vec(corrcoef(roi_avg(:,iroi),sub_avg));
        end;
    end;
end;

figure;
ind = 1;
nsubs = ceil(sqrt(nRoi));
for iroi = 1:nRoi;
    subplot(nsubs,nsubs,ind); 
    imagesq(r(:,:,iroi));
    colormap('hot');
    ind =ind+1;
end

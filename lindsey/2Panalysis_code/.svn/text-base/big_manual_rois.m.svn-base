npix = 9;
roi = [];
nRoi = length(a);
points = [];
for iRoi = 1:nRoi;
    b = a{iRoi};
    points{iRoi} = bsxfun(@plus, b(:,2:3), [1 1]);
    clear b;
    roi_pix = zeros(length(points{iRoi})*npix,2);
    n = length(points{iRoi});
    for ipoint = 1:n;
        for ix = 1:3;
            for iy =1:3;
                b = points{iRoi};
                roi_pix(1+iy-1+(3*(ix-1))+((ipoint-1)*9),:) = [b(ipoint,2)-1+(iy-1),b(ipoint,1)-1+(ix-1)];
                clear b;
            end
        end
    end
    roi{iRoi} = roi_pix;
end

%display regions of interest
figure;
subs = ceil(sqrt(nRoi));
for iRoi = 1:nRoi
    ind = zeros(length(roi{iRoi}),1);
    roi_pix = roi{iRoi};
    for ipix = 1:length(roi{iRoi})
        ind(ipix,1) = sub2ind([240 256],roi_pix(ipix,1),roi_pix(ipix,2));
    end
    clear roi_pix
    fov = zeros(240,256);
    fov(ind)=1;
    subplot(subs,subs,iRoi);    
    imagesq(fov);
    colormap('hot');
end

fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_big_manualrois.mat']);
save(fn_out, 'roi')

fn = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn);

roi_avg = zeros(size(stack,3),nRoi);
for iRoi = 1:nRoi;
    b = roi{iRoi};            
    TC = zeros(length(b), size(stack,3));
    for ipix = 1:length(b);
        TC(ipix,:) = squeeze(stack(b(ipix,1),b(ipix,2),:));
    end
    TC_avg = squeeze(mean(TC,1));
    roi_avg(:,iRoi) = TC_avg';
end

roi_dF = bsxfun(@minus, roi_avg, mean(roi_avg,1));

fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_bigroi_TCs.mat']);
save(fn_out, 'roi_avg');

figure;
tcOffsetPlot(roi_dF);

for iroi1 = 1:nRoi;
    for iroi2 = 1:nRoi;
        r_corr(iroi1, iroi2) = triu2vec(corrcoef(roi_dF(iroi1,:),roi_dF(iroi2,:)));
    end
end
figure;
imagesq(r_corr);
colormap hot;
colorbar;
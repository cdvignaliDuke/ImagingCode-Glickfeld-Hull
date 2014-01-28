%make regions of interest
man_points = bsxfun(@plus, a(:,2:4), [1 1 0]);
points = man_points(:,1:2);
n = length(points);
roi_ic = ceil(man_points(1:4:end,3)/2);
npix = 9;
roi = zeros(length(points)*npix,2);
for ipoint = 1:n;
    for ix = 1:3;
        for iy =1:3;
            roi(1+iy-1+(3*(ix-1))+((ipoint-1)*9),:) = [points(ipoint,2)-1+(iy-1),points(ipoint,1)-1+(ix-1)];
        end
    end
end

%display regions of interest
figure;
nRoi = n/4;
ind = zeros(npix*4,n/4);
roi_stack = zeros(240,256,nRoi);
start =1;
for iroi = 1:nRoi
    for ipix = 1:36
        ind(ipix,iroi) = sub2ind([240 256],roi(ipix+((iroi-1)*36),1),roi(ipix+((iroi-1)*36),2));
    end
    fov = zeros(240,256);
    p = ind(:,iroi);
    fov(p)=1;
    if start == 17;
        figure;
        start = 1;
    end
    subplot(4,4,start);    
    imagesq(fov);
    colormap('hot');
    title([num2str(iroi) '  ic' num2str(roi_ic(iroi))]);
    start = start+1;
    roi_stack(:,:,iroi) = fov;
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_manualrois.mat']);
save(fn_out, 'points', 'roi', 'n', 'npix', 'roi_ic');

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) 'roi_stack.tif']);
writetiff(roi_stack, fn_out);
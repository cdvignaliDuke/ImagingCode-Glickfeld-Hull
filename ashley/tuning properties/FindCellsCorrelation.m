% find cells based on correlation in time and space
b = 5;
siz = size(data_reg);
corr_map = zeros(siz(1),siz(2));
for ix = b:siz(2)-b
    for iy = b:siz(1)-b
        TC = data_reg(iy,ix,:);
        surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
        R = corrcoef(TC,surround);
        corr_map(iy,ix) = R(1,2);
    end
end

figure; imagesq(corr_map); colormap(gray)

bwout = imCellEditInteractive(corr_map);
mask_cell = bwlabel(bwout);
function [x,y] = roiCentersXY(img)
    img(img > 0) = 1;
    mask = bwlabel(img);
    n = max(mask(:));
    
    y = nan(1,n);
    x = nan(1,n);
    for i = 1:n
        [r,c] = find(mask == i);
        y(i) = mean(r);
        x(i) = mean(c);
    end
end
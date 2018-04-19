function cImage = cropImage(img)
avgImg = mean(img,3);
%find coordinates to crop image based on the min and max x/y values of non-zero values
x_max = max(find(sum(avgImg,1)));
x_min = min(find(sum(avgImg,1)));
y_max = max(find(sum(avgImg,2)));
y_min = min(find(sum(avgImg,2)));

cImage = img([y_min:y_max], [x_min:x_max], :);

end
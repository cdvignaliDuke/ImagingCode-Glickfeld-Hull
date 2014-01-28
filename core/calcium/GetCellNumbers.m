function getCellNumbers (labelimg_fname, avgimg_fname)

load(labelimg_fname);
avg_img=imread(avgimg_fname);
avg_img=double(avg_img);
avg_img=avg_img./(max(avg_img(:)));

% get cell numbers by clicking
figure;
imshow(avg_img);
[xi,yi,p]=impixel;
Ncells=length(xi);
cellNumbers=zeros(Ncells,1);
for i=1:Ncells
    cellNumbers(i)=labelimg(yi(i),xi(i));
end

sort(cellNumbers)

figure

cells=zeros(size(labelimg,1),size(labelimg,2));
for i=1:length(cellNumbers)
    if cellNumbers(i)>0
        cells=cells+(labelimg==cellNumbers(i));
    end
end
im=imShade(avg_img, logical(cells));

imagesc(labelimg);
xy=getCellCoordinate(labelimg);
for i=1:length(cellNumbers)
    if cellNumbers(i)>0
        text(xy(1,cellNumbers(i))+5,xy(2,cellNumbers(i)),num2str(cellNumbers(i)));

    end
end



dirName = uigetdir( 'D:\Reid lab\','Directory for saving');
for i=[1:size(QQ,3)]
    if i>100, imwrite(QQ(:,:,i),[dirName '\z' num2str(i) '.tif'],'Compression','none');end;
    if (i<100)&&(i>9), imwrite(QQ(:,:,i),[dirName '\z0' num2str(i) '.tif'],'Compression','none');end;
    if i<10, imwrite(QQ(:,:,i),[dirName '\z00' num2str(i) '.tif'],'Compression','none');end;
end

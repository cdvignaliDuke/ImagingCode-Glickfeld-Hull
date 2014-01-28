function Index=getCellBasedIndexFromPixelMap(pixelmap,labelimg)
% get cell-based Index values from pixelmap
% pixelmap: any pixel map
% labelimg: labels for cells. Size should be the same as pixelmap
%
% 2010.02.15    Kenichi Ohki

Ncells=max(labelimg(:));
Index=zeros(Ncells,1);
for i=1:Ncells
    Index(i)=mean(pixelmap(find(labelimg==i)));
end
return;

function CO_value = getBlobIndex (blobfield, labelimg)

% assigns a number to each cell in a cell map, corresponding to its
% blobness (0=interblob, 1=blob)


%%  scale blobfield.
blobfield=double(blobfield);
blobfield=blobfield-min(blobfield(:));
blobfield=blobfield./max(blobfield(:));

figure(1)
imagesc(blobfield);
axis image;
axis off;

%%  extract blobfield data for each cell

no_cells=max(labelimg(:));
CO_value=zeros(1,no_cells);

for j=1:no_cells
    
    [r,c]=find(labelimg==j);
    ind=[r';c'];
    no_pixels=size(ind,2);

    for i=1:no_pixels
        a=ind(1,i);
        b=ind(2,i);
        c(i)=blobfield(a,b);
    end
    
    CO_value(j)=1-mean(c);
end

%map=ezCellMap(CO_value, labelimg);
%axis image;

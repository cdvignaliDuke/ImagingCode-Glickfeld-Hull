function npmasks = np_maskAK(binarymask,centroid,np_radius, space)

si=size(binarymask);
se = strel('disk',round(np_radius));
Ncells=size(centroid,2);

nhood = ones(space+1,space+1);
largebinary = imdilate(binarymask,nhood);

for i=1:Ncells
   %Create np_dia disk around the cell center
   mask = zeros(si(1),si(2));
   mask(round(centroid(2,i)),round(centroid(1,i))) = 1;
   mask = imdilate(mask,se);

   %subtract enlarged cells
   submask = mask .* largebinary;
   mask = mask - submask;

   npmasks(:,:,i)=mask;
end    
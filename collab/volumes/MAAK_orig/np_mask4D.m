function npmasks = np_mask4D(binarymask,cell_ind,shiftmask,centroid,z_skew,np_radius,space)

si=size(binarymask);
np_nhood = ball_nhood(np_radius,z_skew);
Ncells=size(centroid,2);

sp_nhood = ball_nhood(space,z_skew);
largebinary = imdilate(binarymask,sp_nhood);

for i=1:length(cell_ind)
%     %Create np_dia disk around the cell center
      mask = zeros(si(1),si(2),si(3));
      mask(round(centroid(2,cell_ind(i))),round(centroid(1,cell_ind(i))),round(centroid(3,cell_ind(i)))) = 1;
      mask = imdilate(mask,np_nhood);
% 
%         center=[round(centroid(1,cell_ind(i))),round(centroid(2,cell_ind(i))),round(centroid(3,cell_ind(i)))];
% 
%         xedge=center(1)+[0:(2*np_radius)]-np_radius;
%         yedge=center(2)+[0:(2*np_radius)]-np_radius;
%         zedge=center(3)+[0:(2*round(np_radius./z_skew))]-round(np_radius./z_skew);
% 
%         for j=xedge
%             j
%             for k=yedge
%                 for m=zedge
%                     if sqrt((j-center(1))^2+(k-center(2))^2+((m-center(3)).*z_skew)^2)<np_radius
%                         if min([j k m])>0
%                             mask(j,k,m)=1;
%                         end
%                     end
%                 end
%             end
%         end
        
    %subtract enlarged cells and shiftmask
    submask = mask .* largebinary;
    mask = mask - submask;
    submask = mask .* ~shiftmask;
    mask = mask - submask;
    
    npmasks(:,:,:,i)=mask;
end 

function nhood=ball_nhood(radius,z_skew)
edge=2*radius+1;
zedge=2*round(radius./z_skew)+1;
center=radius+1;
zcenter=round(radius./z_skew)+1;

nhood=zeros(edge,edge,zedge);
for i=1:edge
    for j=1:edge
        for k=1:zedge
            if sqrt((i-center)^2+(j-center)^2+((k-zcenter)*z_skew)^2)<radius
                nhood(i,j,k)=1;
            end
        end
    end
end
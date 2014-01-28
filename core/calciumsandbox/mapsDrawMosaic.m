function mapsDrawMosaic(maps,clim,siz)

[ny,nx,n1,n2]= size(maps);
maps = reshape(maps,ny,nx,n1*n2);

if nargin < 2
    clim = [min(maps(:)) max(maps(:))];
end

if nargin < 3
    if n2<=1
        siz = ceil(sqrt(n1))*[1 1];
    else
        siz = [n1,n2];
    end
end

%% code using mosaic (ok)

mosaic = zeros(ny*siz(1),nx*siz(2));

[height,width]=size(mosaic);

index = 1;
for id2 = 1:siz(2)
    for id1 = 1:siz(1) 
        if index > n1*n2;break;end;
        iy = (id1-1)*ny+1:(id1-1)*ny+ny;
        ix = (id2-1)*nx+1:(id2-1)*nx+nx;
        mosaic(iy,ix)=maps(:,:,index);
        index = index + 1;
    end
end

imshow(mosaic,clim);axis square;%h = colorbar;

if any(siz>1)
    xticks = nx/2:nx:width-nx/2;
    yticks = ny/2:ny:height-ny/2;
    set(gca,'visible','on');
    set(gca,'xtick',xticks,'xticklabel',1:siz(2));
    set(gca,'ytick',yticks,'yticklabel',1:siz(1));
end

box off;
set(gca,'plotboxas',get(gca,'dataas'));

return;

% %% code using subplots (looks like crap)
% figure;
% for itrial = 1:ntrials
%     for istim = 1:nstim 
%         ind = (itrial-1)*nstim+istim;
%         ax(ind)=subplot(ntrials,nstim,ind);
%         imagesc(maps(:,:,itrial,istim),clim);axis square;%h = colorbar;
%         
%         if itrial == 1
%             title(sprintf('S%i',istim));
%         end;
%         
%         if mod(istim,nstim)==1
%             ylabel(sprintf('T%i',itrial));
%         end
%     end
% end
% 
% set(ax,'xticklabel',[],'yticklabel',[]);
% 
% colormap bluered(1);

function im = drawCellMap(labelMask,values,lims,cms,sel,c,plotWhat)
%DRAWCELLMAP
% IM = DRAWCELLMAP(LABELMASK,VALUES,LIMS,CMS,SEL,C)

if nargin < 7 
    plotWhat = {};
end;

nCells = max(labelMask(:));

hold off;

if nargin < 6
    c = [0 0 0]+.1;
end
    
textColor = [1 1 1]*.5;

background = labelMask == 0;

if nargin > 4 & ~isempty(sel)
    background = maskSelect(labelMask,sel)==0;
    selectMask = labelMask .* maskSelect(labelMask,sel);
    unMask = maskSelect(labelMask,setdiff(1:nCells,sel));
else
    selectMask = labelMask;
    unMask = 0*labelMask;
end

if ~iscell(values);    values = {values};end;
if ~iscell(lims);    lims = {lims};end;
if ~iscell(cms); cms = {cms};end;

for ii = 1:length(values)
    % rescale to [-1, 1]
    if isempty(lims{ii});lims{ii}=[-1 1]*prctile(triu2vec(values{ii},1),99);end
    maps(:,:,ii) = imscale(mask2map(labelMask,values{ii}),lims{ii}); 
    inds(:,:,ii) = gray2ind(maps(:,:,ii)); % rescale to [1,64]
    ims(:,:,:,ii) = bsxfun(@times,ind2rgb(inds(:,:,ii),cms{ii}),labelMask>0);
end

im = sum(ims,4);

for iPlane = 1:3
    temp = im(:,:,iPlane);
    temp(find(background))=c(iPlane);
%     if ~isempty(find(unMask))
%         temp(find(unMask))= .8;
%     end
    im(:,:,iPlane) = temp;
end;
    
image(im);axis image;

%axis off;
set(gca,'xtick',[],'ytick',[]);

hold on;
[c,h]=imcontour(labelMask>0,1);
set(h,'color',.2*[1 1 1]);
hold off;

% make up colorbar matching first colormap
hh = colorbar('east');

pos = get(hh,'pos');
newpos = pos;
newpos(3) = pos(3)/2;
newpos(4) = pos(4)/4;
set(hh,'pos',newpos);   

% set(get(hh,'ylabel'),'string',num2str(lims{1}));


% image object of colorbar's colormap
set(hh,'ytick',linspace(0,1,9),'YAxisLocation','right','xcolor',[1 1 1]*.5,'ycolor',[1 1 1]*.5);
ticks = get(hh,'ytick');
set(hh,'ytick',ticks);

set(hh,'yticklabel',ticks * diff(lims{1}) + lims{1}(1) );

imH = findobj(get(hh,'children'),'tag','TMW_COLORBAR'); 
set(imH,'cdata',permute(cms{1},[1 3 2]));

if any(strcmp(plotWhat,'numbers'))
    rp = regionprops(selectMask);
    labList = setdiff(unique(selectMask(:)),0)';

    for iL=labList
        tH(iL) = text(rp(iL).Centroid(1), rp(iL).Centroid(2), ...
                      sprintf('%d',iL), ...
                      'HorizontalAlignment', 'center', ...
                      'VerticalAlignment', 'middle', ...
                      'Tag', 'CellNumbers', ...
                      'Color', textColor);
    end
end;

return;
    
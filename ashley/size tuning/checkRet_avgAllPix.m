ms = '764';
dt = '190423';
t = '1618';
imgFolder = '001';
imgName = '001_000_000';
nfr = 5000;
doShortSession = true;
%%
fn = fullfile('Z:\home\ashley\data',ms,'two-photon imaging',dt,imgFolder);
fnout = fullfile('Z:\home\ashley\Analysis',ms,'two-photon imaging',dt,imgFolder);
cd(fn);
d = squeeze(sbxread(imgName,0,nfr));
% d = squeeze(d(1,:,:,:));
tc = squeeze(mean(mean(d,1),2));
% d = loadsbx_choosepmt(1,ms,dt,imgFolder,imgName);
%%
mw = loadMworksFile(ms,dt,t);
 
%%
on = mw.nScansOn;
off = mw.nScansOff;
% ntrials = length(cell2mat_padded(mw.tGratingAzimuthDeg));
% ntcframes = size(d,3);
% nposstrials = ntcframes./(on+off);
% if ntrials > nposstrials
%     ntrials = floor(nposstrials);
% else
%     tc = tc(1:(ntrials*(on+off)));
% end
%%
taz = cell2mat_padded(mw.tGratingAzimuthDeg);
az = unique(taz);
naz = length(az);
tel = cell2mat_padded(mw.tGratingElevationDeg);
ntr = length(tel);
el = unique(tel);
nel = length(el);
npos = nel.*naz;

tpos = cell(1,ntr);
for iaz = 1:naz
    for iel = 1:nel
        ind = taz == az(iaz) & tel == el(iel);        
        tpos(ind) = {[num2str(az(iaz)) '/' num2str(el(iel))]};
    end
end

% pos_azel = unique(tpos);
pos_azel = cell(1,naz*nel);
for iaz = 1:naz
    if iaz == 1
        ipos = 1;
    end
    for iel = 1:nel
        pos_azel{ipos} = [num2str(az(iaz)) '/' num2str(el(iel))];
        ipos = ipos+1;
    end
end

%%
if doShortSession
    ntr = floor(nfr./(on+off));
    tpos = tpos(1:ntr);
end

%%
tc = tc(1:((on+off).*ntr));
tc_tr = reshape(tc,[on+off,ntr]);
f0 = mean(tc_tr(off/2:off,:),1);
dff = (tc_tr - f0)./f0;

dff_pos = nan(on+off,npos);
for i = 1:npos
    ind = strcmp(tpos,pos_azel{i});
    dff_pos(:,i) = mean(dff(:,ind),2);
end

%%
colors = nan(npos,3);
nbins = ceil(npos/8);
for i = 1:nbins
    if i == 1
        offset = 1;
    else
        offset = offset+7;
    end
    switch i
        case 1
            c = brewermap(13,'RdPu');
        case 2
            c = brewermap(13,'Blues');
        case 3
            c = brewermap(13,'BuPu');
        case 4
            c = brewermap(13,'Greys');
        case 5
            c = brewermap(13,'YlGn');
        case 6
            c = brewermap(13,'Reds');
        case 7
            c = brewermap(13,'Oranges');
    end
    colors(offset:offset+7,:) = c(3:end-3,:);
end
    
setFigParams4Print('portrait')
figure
for i = 1:npos
    hold on
    h = plot(dff_pos(:,i));
    h.Color = colors(i,:);
    h.LineWidth = 3;
end
L = legend(pos_azel,'location','northeastoutside');
L.FontSize = 8;
title(L,'Pos az/el')
if exist(fnout,'dir') == 0
    mkdir(fnout)
end
print(fullfile(fnout,'checkRet_meanPixTC'),'-dpdf','-fillpage')

dff_pos_mean = mean(dff_pos((on/2+off):(off+on),:),1);
dff_pos_sq = reshape(dff_pos_mean,[naz nel]);
pos_names_sq = reshape(pos_azel,[naz nel]);
figure; colormap(brewermap([],'YlOrRd'))
imagesc(flipud(dff_pos_sq))
figXAxis([],'az',[],1:naz,az)
figYAxis([],'el',[],1:nel,fliplr(el))
figAxForm
colorbar
print(fullfile(fnout,'checkRet_meanPixHM'),'-dpdf','-fillpage')
%%
[ypix,xpix,nfr] = size(d);
d = d(:,:,1:(ntr*(on+off)));
d_align = double(reshape(d,[ypix,xpix,on+off,ntr]));
f0 = mean(d_align(:,:,(off/2):off,:),3);
dff = (d_align - f0)./f0;

dff_img_pos = nan(ypix,xpix,npos);
for i = 1:npos
    ind = strcmp(tpos,pos_azel(i));
    dff_img_pos(:,:,i) = mean(mean(dff(:,:,(on/2+off):(off+on),ind),4),3);
end

figure
[nrows, ncols] = optimizeSubplotDim(npos);
for i = 1:npos
    subplot(nrows,ncols,i)
    imagesc(dff_img_pos(:,:,i))
    title(pos_azel(i))
end
print(fullfile(fnout,'checkRet_imgEaPos'),'-dpdf','-fillpage')

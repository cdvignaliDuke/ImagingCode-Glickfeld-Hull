% Plotting for paper figure 2 - size tuning at max contrast
% also does fig 4
% need to load from compareSizeTuneVisAreas_Ai9x, then can run this

%% superplot: example cells, mean curves, proportions, prefSize, SI
expInd = [];
areaInd = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExpRet(i))];
    switch cell2mat(expdata.area(i))
        case 'V1'
            areaInd = [areaInd repmat(1,1,nCellsExpRet(i))];
        case 'LM'
            areaInd = [areaInd repmat(2,1,nCellsExpRet(i))];
        case 'AL'
            areaInd = [areaInd repmat(3,1,nCellsExpRet(i))];
        case 'PM'
            areaInd = [areaInd repmat(4,1,nCellsExpRet(i))];
    end
end
figure(17);clf;
areas = ["V1","LM","AL","PM"];
x = []; yPS = []; ySI = []; ySS = []; x_matched=[]; ySI_matched=[];
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd(goodfit_ind_all),expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
            exCell = 3675;
        case 'LM'
            cutoff = 15; %lm cutoff at 15
            exCell = 2516;
        case 'AL'
            cutoff = 15; %alm cutoff at 15
            exCell = 2147; %
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            exCell = 4243;
    end
    ind_matched = intersect(ind,find(cellDists_all<15)); %15 deg matched size
    ind = intersect(ind,find(cellDists_all<cutoff));
    x = [x; i*ones(size(ind))];
    
    nExpi = length(expIndi);
    nCellsi = length(ind)
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    szs = 5*1.5.^(0:7); nSz=length(szs);
    cons = 0.1*2.^(0:3); nCon=length(cons);
    % this below is for raw data mean curves
    sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
    for iCell = 1:nCellsi
        %dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
        dum = sizeMean(:,nCon,iCell); % only at highest con
        normA = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm(:,nCon,iCell) = sizeMean(:,nCon,iCell)/normA; % normalize by this max for the individual cell
        sizeSEM_norm(:,nCon,iCell) = sizeSEM(:,nCon,iCell)/normA;
    end
    sizeMean_normall = mean(sizeMean_norm,3);
    normA = max(sizeMean_normall);
    sizeMean_normall = sizeMean_normall(:,nCon);%/norm;
    sizeSEM_normall = geomean(sizeSEM_norm,3);%/norm;
    % this below is for fits mean curves
%     sizeMean_norm = repmat(sizeFits(1,1).fitout1*0,size(sizeFits,1),1);
%     for iCell = 1:nCellsi
%         if sizeFits(iCell,nCon).Ftest
%             sizeMean_norm(iCell,:) = sizeFits(iCell,nCon).fitout2;
%         else
%             sizeMean_norm(iCell,:) = sizeFits(iCell,nCon).fitout1;
%         end
%     end
%     sizeMean_normall = mean(sizeMean_norm,1);
%     norm = max(sizeMean_normall);
%     sizeMean_normall = sizeMean_normall/norm;
%     %sizeSEM_normall = geomean(sizeSEM_norm,3)/norm;
%     sizeSEM_normall = std(sizeMean_norm,[],1);%/sqrt(nCellsi);

    cellFit = sizeFits_all(exCell,:);
    subplot(4,3,1+(i-1)*3) % plot example cell
    %errorbar(szs,sizeMean_all(:,nCon,exCell),sizeSEM_all(:,nCon,exCell),'.k')
    plot(sizeFits_all(exCell,nCon).szs0,sizeFits_all(exCell,nCon).data,'.k')
    hold on
    plot(szRng,sizeFits_all(exCell,nCon).fitout1,'b');
    plot(szRng,sizeFits_all(exCell,nCon).fitout2,'r');
%     if cellFit(nCon).Ftest
%         plot(szRng,sizeFits_all(exCell,nCon).fitout2,'k');
%     else
%         plot(szRng,sizeFits_all(exCell,nCon).fitout1,'k');
%     end
    %title(['Example cell from ' char(areas(i))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    maxAmp = max([sizeFits_all(exCell,:).maxResp1 sizeFits_all(exCell,:).maxResp2]);
    %ylim([-0.2*maxAmp 1.2*maxAmp+max(sizeSEM_all(:,nCon,exCell))])
    exLims = [0.16 0.08 0.3 0.4];
    switch i
        case 1
            yticks([-0.04:0.04:0.16])
        case 2
            yticks([-0.02:0.02:0.08])
        case 3
            yticks([-0.1:0.1:0.3])
        case 4
            yticks([-0.1:0.1:0.4])
    end
    ylim([-0.2*maxAmp exLims(i)])
    yLims=get(gca,'ylim');
    xLims=get(gca,'xlim');
    text(xLims(2),yLims(1),char(areas(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
    set(gca,'box','off','TickDir','out')
    if 0;%i==4
        hL = legend(num2str(cons'),'Location','best');
        set(hL,'FontSize',5);
    end
    
    % plot mean size-tune curve
    subplot(4,3,2+(i-1)*3)
    %errorbar(szs,sizeMean_normall(:,nCon),sizeSEM_normall(:,nCon),'k')
    plot(szs,sizeMean_normall,'k')
    %plot(szRng,sizeMean_normall,'k')
    hold on
    upper = sizeMean_normall + sizeSEM_normall;
    lower = sizeMean_normall - sizeSEM_normall;
    border = [upper' fliplr(lower')];
    x2 = [szs fliplr(szs)];
    %upper = sizeMean_normall + sizeSEM_normall;
    %lower = sizeMean_normall - sizeSEM_normall;
    %border = [upper fliplr(lower)];
    %x2 = [szRng fliplr(szRng)];
    fill(x2, border,'k','EdgeColor','none','facealpha',0.2)
    %title({sprintf('Mean Size Tune Curve:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Size (deg)')
    ylabel('dF/F (norm)')
    ylim([0 1.25])
    yticks([0:0.2:1])
    exLims = [1 1 1 1];
    
    ylim([0 exLims(i)])
    yLims=get(gca,'ylim');
    xLims=get(gca,'xlim');
    text(xLims(2),yLims(1),[char(areas(i))],'HorizontalAlignment','right','VerticalAlignment','bottom')
    set(gca,'box','off','TickDir','out')
    if 0%i==4;
        hL = legend(num2str(cons'),'Location','best');
        set(hL,'FontSize',6,'Orientation','horizontal');
    end
    
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    prefMean=zeros(1,nCon);prefSEM=prefMean;
    suppInd = reshape([sizeFits.suppInd],size(sizeFits));
    %suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
    suppMean=zeros(1,nCon);suppSEM=suppMean;
    yPS = [yPS;prefSize(:,nCon)];
    ySI = [ySI;suppInd(:,nCon)];
    ySS = [ySS;reshape([sizeFits(:,nCon).Ftest],size(sizeFits(:,nCon)))];
    
    x_matched = [x_matched; i*ones(size(ind_matched))];
    sizeFits_matched = sizeFits_all(ind_matched,:); %cell,con
    suppInd_matched = reshape([sizeFits_matched.suppInd],size(sizeFits_matched));
    ySI_matched = [ySI_matched;suppInd_matched(:,nCon)];
end

c_areas = categorical({'V1' 'LM' 'AL' 'PM'},{'V1' 'LM' 'AL' 'PM'});
[tbl,chi2stat,pval] = crosstab(x,ySS)
subplot(4,3,3)
phat = tbl(:,1)./sum(tbl,2);
perr = 1.96*sqrt(phat.*(1-phat)./sum(tbl,2));
b=bar(c_areas,[phat,1-phat],'stacked')
b.CData;
hold on
errorbar(1:4,phat,perr,'.k')
hold off
set(gca,'TickLength',[0.015 0.015])
%set(gca,'XTick',1:4,'XTickLabel',{['V1 (n=' num2str(sum(tbl(1,:))) ')'],['LM (n=' num2str(sum(tbl(2,:))) ')'],['AL (n=' num2str(sum(tbl(1,:))) ')'],['PM (n=' num2str(sum(tbl(1,:))) ')']})
set(gca,'box','off','TickDir','out')
ylabel('Frac. cells')
hold on

% Pref Size
[p,~,statPS] = anova1(yPS,x,'off');
[results, means] = multcompare(statPS,'CType','hsd','Display','off');
for i=1:4
    data_cell{i} = [];
    data_cell{i} = yPS(x==i);
end
y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i});
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(4,3,9);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 100 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('Pref size (deg)')

% SI
[p,~,statSI] = anova1(ySI,x,'off');
[results, means] = multcompare(statSI,'CType','hsd','Display','off');
for i=1:4
    data_cell{i} = [];
    ySI_temp = ySI(x==i);
    %ySI_temp = ySI_temp(ySI_temp>0);
    data_cell{i} = ySI_temp;
end
%stats only for SI>0
%ySI2 = [data_cell{1}; data_cell{2}; data_cell{3}; data_cell{4}]';
%x2 = [repmat(1,1,length(data_cell{1})) repmat(2,1,length(data_cell{2})) repmat(3,1,length(data_cell{3})) repmat(4,1,length(data_cell{4}))];
%
y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 2.5]); %,'BoundaryCorrection','reflection'
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(4,3,6);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 1.5 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
ax.XTick = [0:0.5:2.5];
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('SI')

% SI matched
[p,~,statSI_matched] = anova1(ySI_matched,x_matched,'off');
[results, means] = multcompare(statSI_matched,'CType','hsd','Display','off');
for i=1:4
    data_cell{i} = [];
    ySI_temp = ySI_matched(x_matched==i);
    %ySI_temp = ySI_temp(ySI_temp>0);
    data_cell{i} = ySI_temp;
end
%stats only for SI>0
%ySI2 = [data_cell{1}; data_cell{2}; data_cell{3}; data_cell{4}]';
%x2 = [repmat(1,1,length(data_cell{1})) repmat(2,1,length(data_cell{2})) repmat(3,1,length(data_cell{3})) repmat(4,1,length(data_cell{4}))];
%
y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 2.5]); %,'BoundaryCorrection','reflection'
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(4,3,12);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 1.5 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
ax.XTick = [0:0.5:2.5];
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('SI (matched)')

set(gcf,'Color','w')
set(gcf,'Position',[400 0 700 1000])

set(gcf, 'PaperPositionMode', 'auto');
filename = ['K:\ppts\new figs\superfig.eps']
%print(filename, '-depsc2')

%% one way anovas on all cells for size modulation
% for i=1:length(goodfit_ind_size_all)
%     iCell = goodfit_ind_size_all(i);
%     sizeTune = sizeTune_all{:,:,iCell};
%     sizeFit = sizeFits_all(iCell,:);
%     sizeMean = sizeMean_all(:,:,iCell);
%     sizeSEM = sizeSEM_all(:,:,iCell);
%     conStruct = conStruct_all(iCell);
%     
%     figure(1);clf;
%     subplot(2,2,1)
%     title(['cell ' num2str(iCell)])
%     hold on
%     for iCon=1:nCon
%         if ~sizeFit(iCon).Ftest
%             plot(szRng,sizeFit(iCon).fitout1)
%         else
%             plot(szRng,sizeFit(iCon).fitout2)
%         end
%     end
%     line([sizeFit(nCon).prefSize sizeFit(nCon).prefSize],[0 sizeFit(nCon).maxResp1],'LineStyle','--')
%     legend('0.1','0.2','0.4','0.8')
%     subplot(2,2,2)
%     plot(cons,conStruct.resp)
%     subplot(2,2,3)
%     plot(sizeFit(1).szs0,sizeFit(1).data,'k.')
%     hold on
%     errorbar(szs,sizeMean(:,1),sizeSEM(:,1),'.')
%     
%     p1 = anova1(sizeFit(1).data,sizeFit(1).szs0);
%     y = [sizeFit(1).data sizeFit(2).data sizeFit(3).data sizeFit(4).data];
%     g1 = [sizeFit(1).szs0 sizeFit(2).szs0 sizeFit(3).szs0 sizeFit(4).szs0];
%     g2 = [sizeFit(1).szs0*0+1 sizeFit(2).szs0*0+2 sizeFit(3).szs0*0+3 sizeFit(4).szs0*0+4];
%     p2 = anovan(y,{g1, g2},'model','interaction','varnames',{'size','con'},'display','off');
%     
%     figure(1)
%     subplot(2,2,4)
%     plot(1,1,'x')
%     hold on
%     text(1,10,num2str(p1))
%     text(1,6,num2str(p2))
%     ylim([0 12])
%     xlim([0 10])
%     pause
% end
%%
usecon = []; usecon2=[];
Rsqcutoff=0.2;
for i=1:length(goodfit_ind_size_all)
    iCell = goodfit_ind_size_all(i);
    sizeFit = sizeFits_all(iCell,:);
    p1 = anova1(sizeFit(1).data,sizeFit(1).szs0,'off');
    if p1<0.01
        usecon = [usecon i];
    end
    if sizeFit(1).Ftest
        if sizeFit(1).Rsq2>Rsqcutoff
            usecon2 = [usecon2 i];
        end
    else
        if sizeFit(1).Rsq1>Rsqcutoff
            usecon2 = [usecon2 i];
        end
    end
end
temp0 = areaInd(goodfit_ind_all(goodfit_ind_size_all(:)));
temp = areaInd(goodfit_ind_all(goodfit_ind_size_all(usecon)));
temp2 = areaInd(goodfit_ind_all(goodfit_ind_size_all(usecon2)));
temp3 = areaInd(goodfit_ind_all(goodfit_ind_size_all(intersect(usecon,usecon2))));
fprintf('No criteria (all good-size): V1=%d, LM=%d, AL=%d, PM=%d\n',sum(temp0==1),sum(temp0==2),sum(temp0==3),sum(temp0==4));
fprintf('ANOVA (P<0.01) criteria: V1=%d, LM=%d, AL=%d, PM=%d\n',sum(temp==1),sum(temp==2),sum(temp==3),sum(temp==4));
fprintf('Rsq>%0.1f criteria: V1=%d, LM=%d, AL=%d, PM=%d\n',Rsqcutoff,sum(temp2==1),sum(temp2==2),sum(temp2==3),sum(temp2==4));
fprintf('Criteria overlap: V1=%d, LM=%d, AL=%d, PM=%d\n',sum(temp3==1),sum(temp3==2),sum(temp3==3),sum(temp3==4));
%% fig 4
figure(20);clf;
legStrs = strings(1,length(areas));
cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
x = []; yPS = []; ySI = []; ySS = [];
dummyPeak = zeros(4,4);
for i = 1:length(areas)
    figure(20)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    % first find cells out of all in expInd - using goodfit_ind_all inds
    ind = intersect(find(ismember(expInd(goodfit_ind_all),expIndi)),goodfit_ind_size_all(usecon2));
    
    % cutoff by cellDist
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
            exCell = 284;
        case 'LM'
            cutoff = 15; %lm cutoff at 15
            exCell = 1290;
        case 'AL'
            cutoff = 15; %alm cutoff at 15
            exCell = 54;
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            exCell = 2289;
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    x = [x; i*ones(size(ind))];
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nCells_area(i) = nCellsi;
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    ism1_cell{i} = ism1; ism2_cell{i}=ism2;
    
    szs = 5*1.5.^(0:7); nSz=length(szs);
    cons = 0.1*2.^(0:3); nCon=length(cons);
    sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
    for iCell = 1:nCellsi
        %dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
        dum = sizeMean(:,nCon,iCell); % only at highest con
        normA = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/normA; % normalize by this max for the individual cell
        sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/normA;
    end
%     for iCell = 1:nCellsi
%         for iCon = 1:nCon
%             dum = sizeMean(:,iCon,iCell); % only at highest con
%             normA = max(dum(:)); % take max of all dF/F's including all cons
%             sizeMean_norm(:,iCon,iCell) = sizeMean(:,iCon,iCell)/normA; % normalize by this max for the individual cell
%             sizeSEM_norm(:,iCon,iCell) = sizeSEM(:,iCon,iCell)/normA;
%         end
%     end
    sizeMean_normall = mean(sizeMean_norm,3);
    normA = 1;%max(sizeMean_normall(:,nCon)); %not normalizing mean curves
    sizeMean_normall = sizeMean_normall/normA;
    sizeSEM_normall = geomean(sizeSEM_norm,3)/normA;
    
    %dummySI1(:,i) = 1-(sizeMean_normall(end,:) ./ max(sizeMean_normall));
    dummyPeak(i,:) = mean(max(sizeMean,[],1),3);
    
    subplot(5,2,1+(i-1)*2) % plot mean size-tune curve
    for iCon = 1:nCon
        errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon),'Color',[0 0 0]+(0.8-cons(iCon)))
        hold on
    end
    %title({sprintf('Mean Size Tune Curve:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Size (deg)')
    ylabel('dF/F (norm)')
    ylim([0 1])
    yticks([0:0.5:1.5])
    yLims=get(gca,'ylim');
    xLims=get(gca,'xlim');
    text(xLims(2),yLims(1),char(areas(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
    set(gca,'box','off','TickDir','out')
    if i==1;
        hL = legend('\color[rgb]{0.7,0.7,0.7} 0.1','\color[rgb]{0.6,0.6,0.6} 0.2','\color[rgb]{0.4,0.4,0.4} 0.4','\color[rgb]{0,0,0} 0.8','Location','ne');
        set(hL,'FontSize',6); %'Orientation','horizontal'
    end
    
    subplot(5,2,2+(i-1)*2) % plot proportions
    modelcounts = [sum(ism1); sum(ism2)]'/nCellsi;
    phat = modelcounts(:,1);
    perr = 1.96*sqrt(phat.*(1-phat)/nCellsi);
    b = bar(cons_c,modelcounts,'stacked');
    %b(1).CData = repmat([113 238 184]/255,4,1); % seafoam green
    %b(2).CData = repmat([34 139 34]/255,4,1); % forest green -wont work :(
    hold on
    errorbar(1:4,phat,perr,'.k')
    hold off
    %title(sprintf('Model Proportions:%s',areas(i)))
    xlabel('Contrast')
    ylabel('Frac. cells')
    yLims=get(gca,'ylim');
    xLims=get(gca,'xlim');
    text(xLims(2),yLims(2),char(areas(i)),'HorizontalAlignment','center','VerticalAlignment','top')
    set(gca,'box','off','TickDir','out')
    if i==1;
        hL = legend('SS','DOS','location','nw');
        set(hL,'FontSize',6);
    end
    modelcounts(4,2)
    
    subplot(5,2,10) % prefSize
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    prefMean=zeros(1,nCon);prefSEM=prefMean;
    suppInd = reshape([sizeFits.suppInd],size(sizeFits));
    %suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
    suppMean=zeros(1,nCon);suppSEM=suppMean;
    for iCon=1:nCon
        prefMean(iCon) = mean(prefSize(:,iCon));
        prefSEM(iCon) = std(prefSize(:,iCon))./sqrt(nCellsi);
        suppMean(iCon) = mean(suppInd(:,iCon));
        suppSEM(iCon) = std(suppInd(:,iCon))./sqrt(nCellsi);
    end
    errorbar(cons,prefMean,prefSEM);
    set(gca,'box','off','TickDir','out')
    hold on
    if i==4
        %hL = legend(legStrs,'location','best');
        %set(hL,'FontSize',6);
        %title('Pref Size by Area')
        xlabel('Contrast')
        ylabel('PrefSize')
        xlim([0 1])
        ylim([0 60])
    end %'location','southoutside','Orientation','horizontal' for bottom
    yPS = [yPS;prefSize(:,nCon)];
        
    subplot(5,2,9) % SI
    errorbar(cons,suppMean,suppSEM);
    set(gca,'box','off','TickDir','out')
    hold on
    legStrs(i)=sprintf('%s',areas(i)); %  (n=%d), n_{exp}=%d
    if i==4
        hL = legend(legStrs,'location','best');
        set(hL,'FontSize',6);
        %title('SI by area')
        xlabel('Contrast')
        ylabel('SI')
        xlim([0 1])
        ylim([0 1])
    end %'location','southoutside','Orientation','horizontal' for bottom
    ySI = [ySI;suppInd(:,nCon)];
    
    PS_cell{i} = prefSize;
    SI_cell{i} = suppInd;
    
    figure(1);if i==1;clf;end
    h = errorbar(cons,suppMean,suppSEM);
    set(gca,'box','off','TickDir','out')
    hold on
    plot(cons,dummySI1(:,i),'x','Color',h.Color)
    plot(cons,dummySI2(:,i),'^','Color',h.Color)
    plot(cons,dummySI3(:,i),'s','Color',h.Color)
    legStrs(i)=sprintf('%s',areas(i)); %  (n=%d), n_{exp}=%d
    if i==4
        hL = legend(legStrs,'location','best');
        set(hL,'FontSize',6);
        %title('SI by area')
        xlabel('Contrast')
        ylabel('SI')
        xlim([0 1])
        ylim([0 1])
    end %'location','southoutside','Orientation','horizontal' for bottom
end
set(gcf,'Color','w')
set(gcf,'Position',[400 0 400 1000])

set(gcf, 'PaperPositionMode', 'auto');

%%
figure(2);clf;
for i =1:4
    plot(cons,dummyPeak(i,:),'o-')
    hold on
end
title('Average peak response of cells, by contrast and area')
legend(areas)
xlabel('Contrast')
ylabel('Average peak response')

%%
ism2_grid = [ism2_cell{1};ism2_cell{2};ism2_cell{3};ism2_cell{4}];
ism2_locon = ism2_grid(:,1);
ism2_hicon = ism2_grid(:,nCon);

[tbl,chi2stat,pval] = crosstab(x,ism2_locon)
%%
for i = 1:4
    fprintf('Area:%s, n_cells=%d\nLow-con: m1=%d, m2=%d\nHi-con : m1=%d, m2=%d\n',areas(i),sum(x==i),sum(1-ism2_locon(x==i)),sum(ism2_locon(x==i)),sum(1-ism2_hicon(x==i)),sum(ism2_hicon(x==i)))
    [tbl,chi2stat,pval] = crosstab([ism2_locon(x==i)*0; ism2_hicon(x==i)*0+1],[ism2_locon(x==i); ism2_hicon(x==i)])
    pause
end

%%
SI_grid=vertcat(SI_cell{:});
% b = repmat(1:nCon,1,length(SI_grid));
% c = reshape(SI_grid,1,numel(SI_grid));
% z = reshape([x x x x]',1,numel(SI_grid));
% [p tbl stats terms] = anovan(c,{z b},'model','linear','varnames',["area","con"])

PS_grid=vertcat(PS_cell{:});
% b = repmat(1:nCon,1,length(PS_grid));
% c=reshape(PS_grid,1,numel(PS_grid));
% z=reshape([x x x x]',1,numel(PS_grid));
% anovan(c,{z b},'model','full','varnames',["area","con"])
%%
SI_locon = SI_grid(:,1);
[p,tbl,statSIlc]=kruskalwallis(SI_locon,x)
[results, means] = multcompare(statSIlc,'CType','hsd')
%%
PS_locon = PS_grid(:,1);
[p,tbl,statPSlc]=kruskalwallis(PS_locon,x)
[results, means] = multcompare(statPSlc,'CType','hsd')

%% stats on size x con for each area
xData = []; gCon = []; gAr = []; gSz=[];
for i = 1:length(goodfit_ind_size_all)
    if ~sum(usecon2==i)
        continue
    end
    iCell = goodfit_ind_size_all(i);
    sizeFit = sizeFits_all(iCell,:);
    ar = cell2mat(expdata.area(expInd(goodfit_ind_all(iCell))));
    switch ar
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case 'LM'
            cutoff = 15; %alm cutoff at 15
        case 'AL'
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    if cellDists_all(iCell)>=cutoff
        continue
    end
    for iCon = 1:nCon
        dum = sizeFit(iCon).data;
        xData = [xData dum];
        gSz = [gSz sizeFit(iCon).szs0];
        gCon = [gCon iCon+0*dum];
        gAr = [gAr find(strcmp(ar,areas))+0*dum];
    end
end
% 2-way anova by area
i=4; %area
[p tbl stats terms] = anovan(xData(gAr==i),{gCon(gAr==i) gSz(gAr==i)},'model','full','varnames',["con","size"])

%% stats on SI + PS vs con for each area
xSI = []; xPS=[]; gCon = []; gAr = [];
for i = 1:length(goodfit_ind_size_all)
    if ~sum(usecon2==i)
        continue
    end
    iCell = goodfit_ind_size_all(i);
    sizeFit = sizeFits_all(iCell,:);
    ar = cell2mat(expdata.area(expInd(goodfit_ind_all(iCell))));
    switch ar
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case 'LM'
            cutoff = 15; %alm cutoff at 15
        case 'AL'
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    if cellDists_all(iCell)>=cutoff
        continue
    end
    for iCon = 1:nCon
        xSI = [xSI sizeFit(iCon).suppInd];
        xPS = [xPS sizeFit(iCon).prefSize];
        gCon = [gCon iCon];
        gAr = [gAr find(strcmp(ar,areas))];
    end
end
%%
anovan(xSI,{gCon gAr},'model','full','varnames',["con","area"])
%%
anovan(xPS,{gCon gAr},'model','full','varnames',["con","area"])
%%
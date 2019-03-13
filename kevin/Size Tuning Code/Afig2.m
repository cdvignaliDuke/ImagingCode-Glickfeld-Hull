% Plotting for paper figure 2 - size tuning at max contrast
% also does fig 4
% need to load from compareSizeTuneVisAreas_Ai9x, then can run this

%% superplot: example cells, mean curves, proportions, prefSize, SI
figure(17);clf;
areas = ["V1","LM","AL","PM"];
x = []; yPS = []; ySI = []; ySS = [];
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
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
    ind = intersect(ind,find(cellDists_all<cutoff));
    x = [x; i*ones(size(ind))];
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    szs = 5*1.5.^(0:7); nSz=length(szs);
    cons = 0.1*2.^(0:3); nCon=length(cons);
    sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
    for iCell = 1:nCellsi
        dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
        sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
    end
    sizeMean_normall = mean(sizeMean_norm,3);
    norm = max(sizeMean_normall(:,nCon));
    sizeMean_normall = sizeMean_normall/norm;
    sizeSEM_normall = geomean(sizeSEM_norm,3)/norm;
    
    cellFit = sizeFits_all(exCell,:);
    subplot(4,3,1+(i-1)*3) % plot example cell
    errorbar(szs,sizeMean_all(:,nCon,exCell),sizeSEM_all(:,nCon,exCell),'.k')
    hold on
    if cellFit(nCon).Ftest
        plot(szRng,sizeFits_all(exCell,nCon).fitout2,'k');
    else
        plot(szRng,sizeFits_all(exCell,nCon).fitout1,'k');
    end
    %title(['Example cell from ' char(areas(i))])
    xlabel('Size (deg)')
    ylabel('dF/F')
    maxAmp = max([sizeFits_all(exCell,:).maxResp1 sizeFits_all(exCell,:).maxResp2]);
    ylim([-0.2*maxAmp 1.2*maxAmp+max(sizeSEM_all(:,nCon,exCell))])
    yLims=get(gca,'ylim');
    xLims=get(gca,'xlim');
    text(xLims(2),yLims(1),char(areas(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
    set(gca,'box','off','TickDir','out')
    if 0;%i==4
        hL = legend(num2str(cons'),'Location','best');
        set(hL,'FontSize',5);
    end
    
    subplot(4,3,2+(i-1)*3) % plot mean size-tune curve
    %errorbar(szs,sizeMean_normall(:,nCon),sizeSEM_normall(:,nCon),'k')
    plot(szs,sizeMean_normall(:,nCon),'k')
    hold on
    upper = sizeMean_normall(:,nCon) + sizeSEM_normall(:,nCon);
    lower = sizeMean_normall(:,nCon) - sizeSEM_normall(:,nCon);
    border = [upper' fliplr(lower')];
    x2 = [szs fliplr(szs)];
    fill(x2, border,'k','EdgeColor','none','facealpha',0.2)
    %title({sprintf('Mean Size Tune Curve:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Size (deg)')
    ylabel('dF/F (norm)')
    ylim([0 1.3])
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
    suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
    suppMean=zeros(1,nCon);suppSEM=suppMean;
    yPS = [yPS;prefSize(:,nCon)];
    ySI = [ySI;suppInd(:,nCon)];
    ySS = [ySS;reshape([sizeFits(:,nCon).Ftest],size(sizeFits(:,nCon)))];
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

[p,~,statPS] = anova1(yPS,x,'off');
[results, means] = multcompare(statPS,'CType','hsd','Display','off');
figure(17)
subplot(4,3,6) % prefSize
%boxplot(yPS,x,'Labels',cellstr(areas))
%hold on
errorbar(c_areas,means(:,1),means(:,2),'.k');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
%xlabel('Area')
ylabel('PrefSize (deg)')
xlim([0.5 4.5])
ylim([0 50])
% hold on
% plot([3 4],[39 39],'-k', 'LineWidth',1)
% plot([3.3 3.5 3.7],[41 41 41],'*k')
% plot([1 2],[39 39],'-k', 'LineWidth',1)
% plot([1.3 1.5 1.7],[41 41 41],'*k')
% plot([1 3],[44 44],'-k', 'LineWidth',1)
% plot([1.8 2.0 2.2],[46 46 46],'*k')
% plot([2 4],[49 49],'-k', 'LineWidth',1)
% plot([2.8 3.0 3.2],[51 51 51],'*k')
% plot([1 4],[54 54],'-k', 'LineWidth',1)
% plot([2.3 2.5 2.7],[56 56 56],'*k')
% ylim([0 57])

[p,~,statSI] = anova1(ySI,x,'off');
[results, means] = multcompare(statSI,'CType','hsd','Display','off');
figure(17)
subplot(4,3,9) % SI
%boxplot(ySI,x,'Labels',cellstr(areas))
%hold on
errorbar(1:4,means(:,1),means(:,2),'.k');
set(gca,'box','off','TickDir','out')
set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
% legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
% hL = legend(legStrs,'location','best');
% set(hL,'FontSize',6);
%xlabel('Area')
ylabel('SI')
xlim([0.5 4.5])
ylim([0 1])
% hold on
% plot([3 4],[1 1],'-k', 'LineWidth',1)
% plot([3.3 3.5 3.7],[1.05 1.05 1.05],'*k')
% plot([1 2],[1 1],'-k', 'LineWidth',1)
% plot([1.3 1.5 1.7],[1.05 1.05 1.05],'*k')
% plot([1 3],[1.15 1.15],'-k', 'LineWidth',1)
% plot([1.8 2.0 2.2],[1.2 1.2 1.2],'*k')
% plot([2 4],[1.3 1.3],'-k', 'LineWidth',1)
% plot([2.8 3.0 3.2],[1.35 1.35 1.35],'*k')
% plot([1 4],[1.45 1.45],'-k', 'LineWidth',1)
% plot([2.3 2.5 2.7],[1.5 1.5 1.5],'*k')
% ylim([0 1.6])

set(gcf,'Color','w')
set(gcf,'Position',[400 0 700 1000])

set(gcf, 'PaperPositionMode', 'auto');
filename = ['K:\ppts\new figs\superfig.eps']
%print(filename, '-depsc2')

%% fig 4
figure(19);clf;
legStrs = strings(1,length(areas));
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
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
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    szs = 5*1.5.^(0:7); nSz=length(szs);
    cons = 0.1*2.^(0:3); nCon=length(cons);
    sizeMean_norm = sizeMean*0; sizeSEM_norm = sizeSEM*0;
    for iCell = 1:nCellsi
        dum = sizeMean(:,:,iCell); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm(:,:,iCell) = sizeMean(:,:,iCell)/norm; % normalize by this max for the individual cell
        sizeSEM_norm(:,:,iCell) = sizeSEM(:,:,iCell)/norm;
    end
    sizeMean_normall = mean(sizeMean_norm,3);
    norm = max(sizeMean_normall(:,nCon));
    sizeMean_normall = sizeMean_normall/norm;
    sizeSEM_normall = geomean(sizeSEM_norm,3)/norm;
    
    
    subplot(5,2,1+(i-1)*2) % plot mean size-tune curve
    for iCon = 1:nCon
        errorbar(szs,sizeMean_normall(:,iCon),sizeSEM_normall(:,iCon),'Color',[0 0 0]+(0.8-cons(iCon)))
        hold on
    end
    %title({sprintf('Mean Size Tune Curve:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    xlabel('Size (deg)')
    ylabel('dF/F (norm)')
    ylim([0 1.3])
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
    bar(cons_c,modelcounts,'stacked')
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
    
    subplot(5,2,9) % prefSize
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    prefMean=zeros(1,nCon);prefSEM=prefMean;
    suppInd = reshape([sizeFits.suppInd],size(sizeFits));
    suppInd(suppInd<0)=0;suppInd(suppInd>1)=1;
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
        
    subplot(5,2,10) % SI
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
end
set(gcf,'Color','w')
set(gcf,'Position',[400 0 400 1000])

set(gcf, 'PaperPositionMode', 'auto');
% Plotting for paper fig 3 - contrast sensitivity
% need to load up data from compareSizetuneVisAreas_Ai9x.m
% and run contrast response fitting (ends ~line 236~)
%% figure 3
% plot of average contrast response curve in each area, and c50 histogram
% both for @prefSize and @20 deg

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

%close all

legStrs = strings(1,length(areas));
legStrs5=legStrs;legStrs8=legStrs;legStrs10=legStrs;legStrs11=legStrs;
legStrsPC=legStrs;
con_cell = cell(1,4); con_cell20 = con_cell;
x = []; yPS = []; ySI = []; ySS = [];
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd(goodfit_ind_all),expIndi)),goodfit_ind_size_all(:));
    
    % cutoff by cellDist
    % try looking with different cutoffs
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
    %     nExp_area(i) = nExpi;
    %     nCells_area(i) = nCellsi;
    %     sizeMean = sizeMean_all(:,:,ind); % (size,con,cell)
    %     sizeSEM = sizeSEM_all(:,:,ind);
    %     sizeFits = sizeFits_all(ind,:); %cell,con
    
    cons_c = categorical({'0.1' '0.2' '0.4' '0.8'});
    conStruct = conStruct_all(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
    
    
    figure(17);if i==1;clf;end
    %figure 5: average contrast response in each area
    %subplot(2,2,1)
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
    cut = find([conStruct.Rsq]>0.9);
    %cut = find([conStruct.Rsq]);
    legStrs5(i)=sprintf('%s (n=%d)',areas(i),length(cut));
    %legStrs5(i)=sprintf('%s',areas(i));
    conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
    conResp_norm = conResp./conResp(:,nCon);
    con_cell{i} = conResp_norm;
    conMean = mean(conResp_norm,1);
    conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
    %figure(5);if i==1;clf;end
    ax = gca;
    %subplot(2,2,i)
    %for iCell = 1:nCellsi
    %    p1 = plot(cons,conResp_norm(iCell,:),'r-');
    %    p1.Color(4) = 0.1;
    %    hold on
    %end
    hold on
    ax.ColorOrderIndex = i;
    errorbar(cons,conMean,conSEM,'.')
    %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %title('Mean contrast response @prefSize')
    set(gca,'box','off','TickDir','out')
    xlabel('Contrast')
    ylabel('dF/F (norm) @ pref size')
    xlim([0 1])
    ylim([0 1.2])
    if i==4;legend(legStrs5,'location','se');end%'southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
    % fit
    conResp_norm = conResp_norm';
    cRi = conResp_norm(:);
    cons_exp = repmat(cons,1,length(cut));
    lb = [0 0 0.1 1];
    ub = [Inf Inf 0.8 Inf];
    SStot = sum((cRi-mean(cRi)).^2);
    R2best = -Inf;
    x0 = [cRi(1) mean(cRi) 0.2 3]; %BL Rmax C50 n
    [cF, res] = lsqcurvefit(conModelH,x0,cons_exp',cRi,lb,ub,opts);
    R2 = 1-res/SStot;
    
    fitout = conModelH(cF,conRng);
    R50 = fitout(1)+(fitout(end)-fitout(1))/2;
    i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
    C50 = conRng(i50);
    
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(conRng,fitout,'-','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(C50,R50,'x','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    
    figure(18);if i==1;clf;end
    %if sum(choosefig==8) %figure 8: average contrast response at 20 deg in each area
    subplot(2,2,3)
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
    cut = find([conStruct.Rsq20]>0.9);
    legStrs8(i)=sprintf('%s (n=%d)',areas(i),length(cut));
    %legStrs8(i)=sprintf('%s',areas(i));
    conResp = reshape([conStruct(cut).resp20],nCon,length(cut))';
    conResp_norm = conResp./conResp(:,nCon);
    con_cell20{i} = conResp_norm;
    conMean = mean(conResp_norm,1);
    conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
    %figure(8);if i==1;clf;end
    ax = gca;
    %subplot(1,2,i)
    %for iCell = 1:nCellsi
    %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
    %   p1.Color(4) = 0.1;
    %   hold on
    %end
    hold on
    ax.ColorOrderIndex = i;
    errorbar(cons,conMean,conSEM,'.')
    %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %title('Mean contrast response @20deg')
    set(gca,'box','off','TickDir','out')
    xlabel('Contrast')
    ylabel('dF/F (norm) @ 20 deg')
    xlim([0 1])
    ylim([0 1.2])
    if i==4;legend(legStrs8,'location','se');end%'southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
    
    % fit
    conResp_norm = conResp_norm';
    cRi = conResp_norm(:);
    cons_exp = repmat(cons,1,length(cut));
    lb = [0 0 0.1 1];
    ub = [Inf Inf 0.8 Inf];
    SStot = sum((cRi-mean(cRi)).^2);
    R2best = -Inf;
    x0 = [cRi(1) mean(cRi) 0.2 3]; %BL Rmax C50 n
    [cF, res] = lsqcurvefit(conModelH,x0,cons_exp',cRi,lb,ub,opts);
    R2 = 1-res/SStot;
    
    fitout = conModelH(cF,conRng);
    R50 = fitout(1)+(fitout(end)-fitout(1))/2;
    i50 = find(abs(fitout - R50) == min(abs(fitout - R50)),1);
    C50 = conRng(i50);
    
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(conRng,fitout,'-','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(C50,R50,'x','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    
    if i==1
        figure(19);clf
        exCell = 3676;
        cellFit = sizeFits_all(exCell,:);
        % plot example cell
%         for iCon = 1:nCon
%             errorbar(szs,sizeMean_all(:,iCon,exCell),sizeSEM_all(:,iCon,exCell),'.')
%             hold on
%         end
%         ax = gca;
%         ax.ColorOrderIndex = 1;
        maxR=0;max20=0;
        iPref = find(szRng == sizeFits_all(exCell,nCon).prefSize);
        i20 = find(round(szRng)==20,1);
        subplot(1,2,1)
        hold on
        for iCon = 1:nCon
            if cellFit(nCon).Ftest
                plot(szRng,sizeFits_all(exCell,iCon).fitout2,'Color',[0 0 0]+(0.8-cons(iCon)));
                maxR = max([maxR max(sizeFits_all(exCell,iCon).fitout2)]);
                max20 = max([max20 sizeFits_all(exCell,iCon).fitout2(i20)]);
            else
                plot(szRng,sizeFits_all(exCell,iCon).fitout1,'Color',[0 0 0]+(0.8-cons(iCon)));
                maxR = max([maxR max(sizeFits_all(exCell,iCon).fitout1)]);
                max20 = max([max20 sizeFits_all(exCell,iCon).fitout1(i20)]);
            end
        end
        for iCon = 1:nCon
            if cellFit(nCon).Ftest
                plot(szRng(iPref),sizeFits_all(exCell,iCon).fitout2(iPref),'ok');
                plot(szRng(i20),sizeFits_all(exCell,iCon).fitout2(i20),'sk');
            else
                plot(szRng(iPref),sizeFits_all(exCell,iCon).fitout1(iPref),'ok');
                plot(szRng(i20),sizeFits_all(exCell,iCon).fitout1(i20),'sk');
            end
        end
        %title(['Example cell from ' char(areas(i))])
        xlabel('Size (deg)')
        ylabel('dF/F')
        xticks([0:25:100])
        maxAmp = max([sizeFits_all(exCell,:).maxResp1 sizeFits_all(exCell,:).maxResp2]);
        ylim([-0.2*maxAmp 1.2*maxAmp+max(sizeSEM_all(:,nCon,exCell))])
        yLims=get(gca,'ylim');
        xLims=get(gca,'xlim');
        text(xLims(2),yLims(1),char(areas(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
        set(gca,'box','off','TickDir','out')
        
        line(repmat(sizeFits_all(exCell,nCon).prefSize,2), [0 maxR],'Color','r');
        line(repmat(szRng(i20),2), [0 max20],'Color','b');
        hL = legend(num2str(cons'),'Location','best');
        set(hL,'FontSize',8);
        
        axis square
        
        subplot(1,2,2)
        conStruct_ex = conStruct_all(exCell);
        hold on
        line(repmat(conStruct_ex.C50r,2),[0 conModelNB(conStruct_ex.fit,conStruct_ex.C50r)],'Color','r')
        line(repmat(conStruct_ex.C50r20,2),[0 conModelNB(conStruct_ex.fit20,conStruct_ex.C50r20)],'Color','b')
        plot(cons,conStruct_ex.resp,'ok')
        plot(conRng, conModelNB(conStruct_ex.fit,conRng),'r')
        plot(cons,conStruct_ex.resp20,'sk')
        plot(conRng, conModelNB(conStruct_ex.fit20,conRng),'b')
        
        xlim([0 1])
        xlabel('Contrast')
        xticks([0:.25:1])
        ylabel('dF/F')
        
        yLims=get(gca,'ylim');
        xLims=get(gca,'xlim');
        text(xLims(2),yLims(1),char(areas(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
        set(gca,'box','off','TickDir','out')
        legend(["Pref Size","20 deg"])
        axis square
        
    end
end

%% stats on C50, BLr
areas = ["V1","LM","AL","PM"];
Rsqcutoff = 0.9;
x =[]; yC50=x;
x20 =[]; yC5020=x;
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
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    sizeFits = sizeFits_all(ind,:); %cell,con
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    conStruct = conStruct_all(ind);
    C50r = [conStruct.C50r]';
    C50r20 = [conStruct.C50r20]';
    
    %prefCut = find((prefSize(:,nCon)>10).*(prefSize(:,nCon)<30));
    
    Rsq = [conStruct.Rsq];
    cut = find(Rsq>Rsqcutoff);
    %cut = intersect(cut,prefCut);
    nCut(i) = length(cut);
    Rsq20 = [conStruct.Rsq20];
    cut20 = find(Rsq20>Rsqcutoff);
    %cut20 = intersect(cut20,prefCut);
    nCut20(i) = length(cut20);
    
    x = [x i+0*cut];
    yC50 = [yC50; C50r(cut)];
    x20 = [x20 i+0*cut20];
    yC5020 = [yC5020; C50r20(cut20)];
    
end

c_areas = categorical({'V1' 'LM' 'AL' 'PM'},{'V1' 'LM' 'AL' 'PM'});
fprintf(['C50 @pS, nCut=' num2str(nCut)])
[p,~,statC50] = anova1(yC50,x,'off')
[p,tbl,statC50]=kruskalwallis(yC50,x)
[results, means] = multcompare(statC50,'CType','hsd','Display','off')
figure(17)
% subplot(2,2,2) % C50 @prefSize
% %boxplot(yC50,x,'Labels',cellstr(areas))
% %hold on
% errorbar(c_areas,means(:,1),means(:,2),'.k');
% set(gca,'box','off','TickDir','out')
% set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
% xlabel('Area')
% ylabel('C_{50}')
% xlim([0.5 4.5])
% ylim([0 1])
%5 sideways violin
for i=1:4
    data_cell{i} = [];
    data_cell{i} = yC50(x==i);
end
y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i},'BoundaryCorrection','reflection','Support',[0-eps 1+eps]); %,'BoundaryCorrection','reflection'
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(1,1,1);%2,2,2);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 1 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('C50 @prefSize')

fprintf(['C50 @20d, nCut=' num2str(nCut20)])
[p,~,statC5020] = anova1(yC5020,x20,'off')
[results, means] = multcompare(statC5020,'CType','hsd','Display','off')
figure(18)
% subplot(2,2,4) % C50 @20deg
% %boxplot(yC5020,x20,'Labels',cellstr(areas))
% %hold on
% errorbar(c_areas,means(:,1),means(:,2),'.k');
% set(gca,'box','off','TickDir','out')
% set(gca,'XTick',1:4,'XTickLabel',areas,'TickLength',[0.015 0.015])
% xlabel('Area')
% ylabel('C_{50}')
% xlim([0.5 4.5])
% ylim([0 1])
%5 horizontal violin
for i=1:4
    data_cell{i} = [];
    data_cell{i} = yC5020(x20==i);
end
y_mean = [mean(data_cell{1}) mean(data_cell{2}) mean(data_cell{3}) mean(data_cell{4})];
y_std = [std(data_cell{1}) std(data_cell{2}) std(data_cell{3}) std(data_cell{4})];
for i = 1:length(areas)
    [fi xi] = ksdensity(data_cell{i});
    fnorm(:,i) = fi/max(fi)*0.3;
    xinorm(:,i) = xi;
end
ax = subplot(2,2,4);
colors = get(ax,'ColorOrder');
for i=1:length(areas)
    hold on
    h5(i)=fill([xinorm(:,i);flipud(xinorm(:,i))],[fnorm(:,i)+(5-i);flipud((5-i)-fnorm(:,i))],[1 1 1],'EdgeColor','k');
    p(1)=plot([y_mean(i) y_mean(i)],[interp1(xinorm(:,i),fnorm(:,i)+(5-i),y_mean(i)), interp1(flipud(xinorm(:,i)),flipud((5-i)-fnorm(:,i)),y_mean(i)) ],'k','LineWidth',2);
    h5(i).FaceColor = colors(i,:);
end
axis([0 1 0.5 length(areas)+0.5]);
legend off
ax.YTick = [1:4];
ax.YTickLabel = fliplr(areas);
set(gca,'box','off','TickDir','out')
ylabel('Area')
xlabel('C50 @20deg')

set(gcf,'Color','w')
set(gcf,'Position',[400 0 800 600])

%% two way anova for contrast response (con x area)
% generate data
conStruct_new = conStruct_all(goodfit_ind_size_all);
xCon = []; gCon = []; gAr = []; xC50=[]; gAr2=[];
xCon20 = []; gCon20 = []; gAr20 = []; xC5020=[]; gAr202=[];
for iCell=1:length(conStruct_new)
    ar = cell2mat(expdata.area(expInd(goodfit_ind_all(goodfit_ind_size_all(iCell)))));
    switch ar
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
        case {'LM','AL'}
            cutoff = 15; %alm cutoff at 15
        case 'PM'
            cutoff = 20; %pm cutoff at 20
    end
    if cellDists_all(goodfit_ind_size_all(iCell))>=cutoff
        continue
    end
    if conStruct_new(iCell).Rsq>0.9
        dum = conStruct_new(iCell).resp;
        xCon = [xCon dum./max(dum)];
        gCon = [gCon cons];
        gAr = [gAr find(strcmp(ar,areas))+0*dum];
        xC50 = [xC50 conStruct_new(iCell).C50r];
        gAr2 = [gAr2 find(strcmp(ar,areas))];
    end
    if conStruct_new(iCell).Rsq20>0.9
        dum20 = conStruct_new(iCell).resp20;
        xCon20 = [xCon20 dum20./max(dum20)];
        gCon20 = [gCon20 cons];
        gAr20 = [gAr20 find(strcmp(ar,areas))+0*dum];
        xC5020 = [xC5020 conStruct_new(iCell).C50r20];
        gAr202 = [gAr202 find(strcmp(ar,areas))];
    end
end
%%
% 2 way anova with grouping by con and area
[p tbl stats terms] = anovan(xCon,{gCon gAr},'model','full','varnames',["con","area"])
% [p tbl stats terms] = anovan(xSz,{gCon gAr},'model','linear','varnames',["con","area"]) % for no interaction term
[p tbl stats terms] = anovan(xCon20,{gCon20 gAr20},'model','full','varnames',["con","area"])
%% tests for homogeneity
vartestn(xCon',gCon','TestType','BrownForsythe') %LeveneAbsolute %BrownForsythe
vartestn(xCon',gAr','TestType','BrownForsythe')
vartestn(xCon20',gCon20','TestType','BrownForsythe')
vartestn(xCon20',gAr20','TestType','BrownForsythe')
vartestn(xCon',gCon','TestType','LeveneAbsolute')
vartestn(xCon',gAr','TestType','LeveneAbsolute')
vartestn(xCon20',gCon20','TestType','LeveneAbsolute')
vartestn(xCon20',gAr20','TestType','LeveneAbsolute')
%% C50 Kruskal wallis

[p,tbl,stats]=kruskalwallis(xC50,gAr2)
%[p,~,stats] = anova1(yC5020,xC50,'off')
[results, means] = multcompare(stats,'CType','hsd','Display','off')
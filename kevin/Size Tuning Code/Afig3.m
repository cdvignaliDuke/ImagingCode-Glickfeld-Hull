% Plotting for paper fig 3 - contrast sensitivity
% need to load up data from compareSizetuneVisAreas_Ai9x.m
% and run contrast response fitting (ends ~line 236~)
%% figure 3
% plot of average contrast response curve in each area, and c50 histogram
% both for @prefSize and @20 deg
figure(18);clf;

fprintf('Examine cells from each area:\n')
areas = ["V1","LM","AL","PM"];

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

%close all

legStrs = strings(1,length(areas));
legStrs5=legStrs;legStrs8=legStrs;legStrs10=legStrs;legStrs11=legStrs;
legStrsPC=legStrs;
x = []; yPS = []; ySI = []; ySS = [];
for i = 1:length(areas)
    fprintf(['Area #' num2str(i) ' : ' char(areas(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
    
    % cutoff by cellDist
    % try looking with different cutoffs
    switch areas(i)
        case 'V1'
            cutoff = 10; %v1 cutoff at 10
            excells = [631 2128];
            % case {'LM','AL'}
            %     cutoff = 15; %alm cutoff at 15
        case 'LM'
            cutoff = 10; %alm cutoff at 15
            excells = [1861 1863];
        case 'AL'
            cutoff = 10; %alm cutoff at 15
            excells = [2395 1777];
        case 'PM'
            cutoff = 10; %pm cutoff at 20
            excells = [1952 2292];
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
    
    %if sum(choosefig==5) %figure 5: average contrast response in each area
    subplot(2,2,1)
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
    cut = find([conStruct.Rsq]>0.9);
    legStrs5(i)=sprintf('%s (n=%d)',areas(i),length(cut));
    conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
    conResp_norm = conResp./conResp(:,nCon);
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
    errorbar(cons,conMean,conSEM)
    %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %title('Mean contrast response @prefSize')
    set(gca,'box','off','TickDir','out')
    xlabel('Contrast')
    ylabel('norm. dF/F @ pref size')
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
    plot(conRng,fitout,':','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(C50,R50,'x','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    
    
    %if sum(choosefig==8) %figure 5: average contrast response at 20 deg in each area
    subplot(2,2,3)
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
    cut = find([conStruct.Rsq20]>0.9);
    legStrs8(i)=sprintf('%s (n=%d)',areas(i),length(cut));
    conResp = reshape([conStruct(cut).resp20],nCon,length(cut))';
    conResp_norm = conResp./conResp(:,nCon);
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
    errorbar(cons,conMean,conSEM)
    %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    %title('Mean contrast response @20deg')
    set(gca,'box','off','TickDir','out')
    xlabel('Contrast')
    ylabel('norm. dF/F @20deg')
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
    plot(conRng,fitout,':','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot(C50,R50,'x','HandleVisibility','off')
    ax = gca;
    ax.ColorOrderIndex = i;
    plot([C50 C50],[0 R50],'--','HandleVisibility','off')
    
end

%% stats on C50, BLr
areas = ["V1","LM","AL","PM"];
Rsqcutoff = 0.9;
x =[]; yC50=x; yBLr=x;
x20 =[]; yC5020=x; yBLr20=x;
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
    C50f = 0*ind;
    %baseline = C50f;
    %Rmax = C50f;
    BLr = C50f; C50f20=C50f; BLr20=C50f;
    for iCell = 1:nCellsi
        C50f(iCell) = conStruct(iCell).fit(3);
        %baseline(iCell) = conStruct(iCell).fit(1);
        %Rmax(iCell) = conStruct(iCell).fit(2);
        %BLr = baseline(cut)./(baseline(cut)+Rmax(cut));
        BLr(iCell) = conStruct(iCell).fit(1)./(conStruct(iCell).fit(1)+conStruct(iCell).fit(2));
        C50f20(iCell) = conStruct(iCell).fit20(3);
        BLr20(iCell) = conStruct(iCell).fit20(1)./(conStruct(iCell).fit20(1)+conStruct(iCell).fit20(2));
    end
    C50r = [conStruct.C50r]';
    C50r20 = [conStruct.C50r20]';
    
    prefCut = find((prefSize(:,nCon)>10).*(prefSize(:,nCon)<30));
    
    Rsq = [conStruct.Rsq];
    cut = find(Rsq>Rsqcutoff);
    cut = intersect(cut,prefCut);
    nCut(i) = length(cut);
    Rsq20 = [conStruct.Rsq20];
    cut20 = find(Rsq20>Rsqcutoff);
    cut20 = intersect(cut20,prefCut);
    nCut20(i) = length(cut20);
    
    x = [x; i+0*cut];
    yC50 = [yC50; C50r(cut)];
    yBLr = [yBLr; BLr(cut)];
    x20 = [x20; i+0*cut20];
    yC5020 = [yC5020; C50r20(cut20)];
    yBLr20 = [yBLr20; BLr20(cut20)];
    
end

c_areas = categorical({'V1' 'LM' 'AL' 'PM'},{'V1' 'LM' 'AL' 'PM'});
fprintf(['C50 @pS, nCut=' num2str(nCut)])
[p,~,statC50] = anova1(yC50,x,'off')
[results, means] = multcompare(statC50,'CType','hsd','Display','off')
figure(18)
subplot(2,2,2) % C50 @prefSize
boxplot(yC50,x,'Labels',cellstr(areas))
hold on
errorbar(c_areas,means(:,1),means(:,2),'.k');
set(gca,'box','off','TickDir','out')
xlabel('Area')
ylabel('C_{50}')
%xlim([0 1])
ylim([0 1])

fprintf(['C50 @20d, nCut=' num2str(nCut20)])
[p,~,statC5020] = anova1(yC5020,x20,'off')
[results, means] = multcompare(statC5020,'CType','hsd','Display','off')
figure(18)
subplot(2,2,4) % C50 @20deg
boxplot(yC5020,x20,'Labels',cellstr(areas))
hold on
errorbar(c_areas,means(:,1),means(:,2),'.k');
set(gca,'box','off','TickDir','out')
xlabel('Area')
ylabel('C_{50}')
%xlim([0 1])
ylim([0 1])

set(gcf,'Color','w')
set(gcf,'Position',[400 0 800 600])
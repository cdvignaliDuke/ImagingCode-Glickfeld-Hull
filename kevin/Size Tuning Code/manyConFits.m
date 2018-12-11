% test script to fit contrast response curves at many sizes
% plan:
% collect mean fit parameters including C50 and BLr at all sizes 5-35 (can adjust)

%% begin with sizes 5-35
szTakes = 5:35;

s4 = zeros(1,4);
s = zeros(1);
conStructI = struct('resp',s4,'resp20',s4,'fit',s4,'C50r',s,'Rsq',s,'x0',s4);
conStructI(nCellsTot) = conStructI;
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
conRng = 0.001:0.001:1;
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'

C50f_all = zeros(length(szTakes),2); %2 for V1+PM
C50r_all = C50f_all; BLr_all = C50f_all;
for iSz = 1:length(szTakes)
    sz = szTakes(iSz);
    
    % con resp for all cells
    fprintf('\nExtracting at size: %d...', sz)
    conStructX = conStructI;
    
    parfor iCell=1:nCellsTot
        if ~sum(iCell==goodfit_ind_size_all)
            if sum(iCell==[1 nCellsTot]) % check first and last to reset zeros to blank
                conStructX(iCell).resp = [];conStructX(iCell).fit = [];conStructX(iCell).C50r = [];conStructX(iCell).Rsq = [];conStructX(iCell).x0 = [];
            end
            continue % do not fit unless goodfit_size
        end
        
        indSz = find(min(abs(szRng-sz))==abs(szRng-sz),1);
        for iCon = 1:nCon
            if sizeFits_all(iCell,iCon).Ftest
                conStructX(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout2(indSz);
            else
                conStructX(iCell).resp(iCon) = sizeFits_all(iCell,iCon).fitout1(indSz);
            end
        end
        
        % fit at given size
        cRi = conStructX(iCell).resp;
        lb = [0 0 0.1 1];
        ub = [Inf Inf 0.8 Inf];
        SStot = sum((cRi-mean(cRi)).^2);
        R2best = -Inf;
        for i=1%1:4
            x0 = [cRi(1) max(cRi) 0.1+0.1*i 3]; %BL Rmax C50 n
            [cF, res] = lsqcurvefit(conModelH,x0,cons,cRi,lb,ub,opts);
            R2 = 1-res/SStot;
            if R2>R2best
                R2best = R2;
                cFbest = cF;
                x0best = x0;
            end
        end
        cF = cFbest;
        R2 = R2best;
        conStructX(iCell).fit = cF;
        conStructX(iCell).Rsq = R2;
        conStructX(iCell).x0 = x0best;
    end
    fprintf('done, taking means')
    
    conRng = 0.001:0.001:1;
    opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
    for i = 1:length(areas)
        fprintf(['Area #' num2str(i) ' : ' char(areas(i))])
        % select exps matching area
        expIndi = find(cellfun(@(x) strcmp(x,areas(i)), expdata.area, 'UniformOutput', 1));
        % find cells with correct exp inds, take only good fit cells
        ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_size_all);
        
        % cutoff by cellDist
        % try looking with different cutoffs
        switch areas(i)
            case 'V1'
                cutoff = 10; %v1 cutoff at 10
                excells = [515 529]; %m1:515, 1153; m2:374
            case 'PM'
                cutoff = 20; %pm cutoff at 20
                excells = [1029 652]; %m1: 590; m2:577, 591
        end
        ind = intersect(ind,find(cellDists_all<cutoff));
        
        nExpi = length(expIndi);
        nCellsi = length(ind);
        nExp_area(i) = nExpi;
        nCells_area(i) = nCellsi;
        
        conStruct = conStructX(ind);
        
        legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
        
        % plot mean contrast response curves and then make fits
        cut = find([conStruct.Rsq]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(5);if i==1;clf;end
        ax = gca;
        %subplot(1,2,i)
        for iCell = 1:length(cut)
            if i==1
                p1 = plot(cons,conResp_norm(iCell,:),'r-');
                p1.Color(4) = 0.05;
            else
                p1 = plot(cons,conResp_norm(iCell,:),'b-');
                p1.Color(4) = 0.05;
            end
            hold on
        end
        hold on
        if i==1
            errorbar(cons,conMean,conSEM,'k')
        else
            errorbar(cons,conMean,conSEM,'k--')
        end
        ax.ColorOrderIndex = i;
        %title({sprintf('Mean con resp @%d deg, Area:%s',sz,areas(i));['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})
        title(sprintf('Mean con resp @%d deg',sz))
        xlabel('Contrast')
        ylabel('norm. dF/F @ pref size')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
        
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
        
        C50f_all(iSz,i) = cF(3);
        C50r_all(iSz,i) = C50;
        BLr_all(iSz,i) = cF(1)/(cF(1)+cF(2));
        nCells_all(iSz,i) = length(cut);
    end
    pause(0.01)
end

figure(4);clf;
subplot(2,2,1)
plot(szTakes,C50f_all,'-')
title('C50f')
xlabel('Matched size (deg)')
subplot(2,2,3)
plot(szTakes,C50r_all,'-')
title('C50r')
xlabel('Matched size (deg)')
subplot(2,2,2)
plot(szTakes,BLr_all,'-')
title('BLr')
xlabel('Matched size (deg)')
subplot(2,2,4)
plot(szTakes,nCells_all,'-')
title('nCells')
xlabel('Matched size (deg)')
legend(["V1","PM"])

%%
fprintf('Examine cells from each area:\n')
areas = ["V1","PM"];

nExp_area = zeros(size(areas));
nCells_area = nExp_area;

%close all
choosefig = [3 5 8 9 10 11];
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd;
% 5=conresp; 6=ex.cells; 7=medianfits; 8=conresp matched size @20deg,
% 9=prefSize but PS within 10-30, 10=conresp but PS within 10-30
% 11= conResp @20 and PS within 10-30deg
legStrs = strings(1,length(areas)); legStrs2=legStrs; legStrsPC=legStrs;
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
            excells = [515 529]; %m1:515, 1153; m2:374
        case 'PM'
            cutoff = 20; %pm cutoff at 20
            excells = [1029 652]; %m1: 590; m2:577, 591
    end
    ind = intersect(ind,find(cellDists_all<cutoff));
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_area(i) = nExpi;
    nCells_area(i) = nCellsi;
    sizeTune = sizeTune_all{:,:,ind}; % (size,con,cell)
    sizeMean = sizeMean_all(:,:,ind);
    sizeSEM = sizeSEM_all(:,:,ind);
    sizeFits = sizeFits_all(ind,:); %cell,con
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    ism1 = reshape(~[sizeFits.Ftest],size(sizeFits));
    ism2 = reshape([sizeFits.Ftest],size(sizeFits));
    
    prefSize = reshape([sizeFits.prefSize],size(sizeFits));
    
    cons_c = categorical({'0.05' '0.1' '0.2' '0.4' '0.8' '1'});
    conStruct = conStructX(ind);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nCellsi,nExpi);
    
    if sum(choosefig==1)
        % C50 of mean fit in each area across various sizes
        figure(1);if i==1;clf;end
        histogram(szTakes, c50_bySize)
        
    end
    
    if sum(choosefig==5) %figure 5: average contrast response in each area
        
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
        fitout50rect = abs(fitout - R50);
        i50 = find(fitout50rect == min(fitout50rect),1);
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
    
    if sum(choosefig==6) %figure 6: example cells from each area, with fits
        figure(6);if i==1;clf;end
        subplot(2,2,2*(i-1)+1)
        dum = sizeMean_all(:,:,excells(1)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(1))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(1))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Non-suppressed cell in %s',areas(i)))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        if sum(i==[1 3]);ylabel('dF/F (norm)');end
        ylim([0 1.2])
        subplot(2,2,2*(i-1)+2)
        dum = sizeMean_all(:,:,excells(2)); % take all sizeMean values for cell
        %dum = sizeMean(:,nCon,iCell); % only at highest con
        norm = max(dum(:)); % take max of all dF/F's including all cons
        sizeMean_norm = sizeMean_all(:,:,excells(2))/norm; % normalize by this max for the individual cell
        sizeSEM_norm = sizeSEM_all(:,:,excells(2))/norm;
        for iCon = 1:nCon
            errorbar(szs,sizeMean_norm(:,iCon),sizeSEM_norm(:,iCon))
            hold on
        end
        title(sprintf('Suppressed cell in %s',areas(i)))
        if sum(i==[3 4]);xlabel('Size (deg)');end
        %ylabel('dF/F (norm)')
        ylim([0 1.2])
        if i==2;legend(num2str(cons'));end
    end
    
    if sum(choosefig==7) % contrast-resp fits
        % look at different data:
        % histogram of iC50 (best guess) vs C50fit and C50 rec
        % cross plots of these, C50rec-C50fit
        % Rsq histogram, choose cutoff, look at ratio by area
        % C50rec across areas (boxplot + mean)
        C50f = 0*ind;
        for iCell=1:length(ind)
            C50f(iCell) = conStruct(iCell).fit(3);
        end
        C50r = [conStruct.C50r];
        Rsq = [conStruct.Rsq];
        cut = find(Rsq>0.9);
        figure(7);if i==1;clf;end
        subplot(1,2,i)
        plot(C50f,C50r,'.')
        xlabel('C50f')
        ylabel('C50r')
        title({sprintf('Area:%s',areas(i));['(n=' num2str(length(cut)) ', n_{exp}=' num2str(nExpi) ')']})
    end
    
    if sum(choosefig==8) %figure 5: average contrast response at 20 deg in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(8);if i==1;clf;end
        ax = gca;
        ax.ColorOrderIndex = i;
        %subplot(1,2,i)
        %for iCell = 1:nCellsi
        %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %   p1.Color(4) = 0.1;
        %   hold on
        %end
        hold on
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response @20deg')
        xlabel('Contrast')
        ylabel('norm. dF/F @20deg')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
        
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
    
    prefCut = find((prefSize(:,6)>10).*(prefSize(:,6)<30));
    nPC = length(prefCut);
    legStrsPC(i)=sprintf('%s (n=%d, n_{exp}=%d)',areas(i),nPC,nExpi);
    if sum(choosefig==9)
        figure(9);if i==1;clf;end %figure 9 = prefSize vs con, cut 10-30
        
        prefMean=zeros(1,nCon);prefSEM=prefMean;
        for iCon=1:nCon
            prefMean(iCon) = mean(prefSize(prefCut,iCon));
            prefSEM(iCon) = std(prefSize(prefCut,iCon))./sqrt(nPC);
        end
        errorbar(cons,prefMean,prefSEM);
        hold on
        title('Mean Preferred Size by Area, only 10-30 deg prefSize')
        xlabel('Contrast')
        ylabel('PrefSize')
        xlim([0 1])
        ylim([0 60])
        if i==2;legend(legStrsPC,'location','southoutside','Orientation','horizontal');end %'location','eastoutside' %'location','southoutside','Orientation','horizontal' for bottom
    end
    
    if sum(choosefig==11) %figure 5: average contrast response at 20 deg in each area
        conRng = 0.001:0.001:1;
        opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
        cut = find([conStruct.Rsq]>0.9);
        cut = intersect(cut,prefCut);
        legStrs2(i)=sprintf('%s (n=%d)',areas(i),length(cut));
        conResp = reshape([conStruct(cut).resp],nCon,length(cut))';
        conResp_norm = conResp./conResp(:,nCon);
        conMean = mean(conResp_norm,1);
        conSEM = std(conResp_norm,[],1)./sqrt(length(cut));
        figure(11);if i==1;clf;end
        ax = gca;
        ax.ColorOrderIndex = i;
        %subplot(1,2,i)
        %for iCell = 1:nCellsi
        %   p1 = plot(cons,conResp_norm(iCell,:),'r-');
        %   p1.Color(4) = 0.1;
        %   hold on
        %end
        hold on
        errorbar(cons,conMean,conSEM)
        %title({sprintf('Contrast response - Area:%s',areas(i));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
        title('Mean contrast response @20deg, only 10-30 deg prefSize')
        xlabel('Contrast')
        ylabel('norm. dF/F @20deg')
        xlim([0 1])
        ylim([0 1.2])
        if i==2;legend(legStrs2,'location','southoutside','Orientation','horizontal');end %'location','southoutside','Orientation','horizontal' for bottom
        
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
        fitout50rect = abs(fitout - R50);
        i50 = find(fitout50rect == min(fitout50rect),1);
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
end

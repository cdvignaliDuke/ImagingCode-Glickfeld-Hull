SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;

x = [];
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));

TFminmax0 =  log2([1 15]);
SFminmax0 =  log2([.02 .32]);
        
TFminmax2 = [TFminmax0(1)-4 TFminmax0(2)+4];
SFminmax2 = [SFminmax0(1)-4 SFminmax0(2)+4];
% h1 = figure;

area_order = [2;3;1];
expt = [8; 4; 6];

h1 = figure;
h2 = figure;


for iArea = 3
    area = areas(area_order(iArea),:);
    iexp = expt(area_order(iArea));
    n = all_fits(area_order(iArea)).expt(iexp).n(1);
    ind = [];
    for iCell = 1:n
        if all_fits(area_order(iArea)).expt(iexp).bouton(iCell).goodfit == 1
            ind = [ind iCell];
        end
    end
    for iCell = 1:length(ind)
        fprintf([num2str(iCell) ' ']);
        xfit = [all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).dF_fit all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).sigma_SF all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).sigma_TF log2(all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).SF_fit) log2(all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).TF_fit) all_fits(area_order(iArea)).expt(iexp).bouton(ind(iCell)).xi_fit exp(-1*(1/8))];
        figure(h1);
        h0 = ezplot(@(x,y) Gauss2D_ellipseMA_forplotting(x,y,xfit),TFminmax2,SFminmax2);
        title(area);
        xlim(TFminmax0);
        ylim(SFminmax0);
        axis square;
        hold on
        hold on
        XData0 = get(h0,'XData');
        YData0 = get(h0,'YData');
        figure(h2);
        if length(XData0)>2
            if length(XData0)>10 & iscell(XData0)==0
                XData = XData0; 
                YData = YData0;
            else
                length(XData0)
                XData = [cell2mat(XData0(1))]; % cell2mat(XData0(2))]; 
                YData = [cell2mat(YData0(1))]; % cell2mat(YData0(2))]; 
            end

            h3 = patch(XData([1:2:end 1]),YData([1:2:end 1]),'k');
            set(h3,'FaceColor',[0 0 0],'FaceAlpha',.01,'EdgeColor',[1 1 1]*.4);
            axis([min(TFminmax0) max(TFminmax0) min(SFminmax0) max(SFminmax0)]);
            hold on
            axis square
            title(area);
            xlim(TFminmax0);
            ylim(SFminmax0);
        end   
    end
end


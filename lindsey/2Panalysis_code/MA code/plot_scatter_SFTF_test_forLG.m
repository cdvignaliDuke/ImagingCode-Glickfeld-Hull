%now plot SF vs TF, AL and PM on same axis, V1 separately, plot contour
%plots

%use same criteria as for speed tuning: 
ind_AL = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==1 & Mask_Sigma & Mask_rej2 & Mask_noextremeSFTF);
ind_PM = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==2 & Mask_Sigma & Mask_rej2 & Mask_noextremeSFTF);
ind_V1 = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==3 & Mask_Sigma & Mask_rej2 & Mask_noextremeSFTF);

ind_AL = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==1  & Mask_rej2 & Mask_noextremeSFTF);
ind_PM = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==2  & Mask_rej2 & Mask_noextremeSFTF);
ind_V1 = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==3  & Mask_rej2 & Mask_noextremeSFTF);

SF_vec01 = [.01 .02 .04 .08 .16 .32]';
TF_vec01 = [.25 .5 1 2 4 8 15 24]';
[sfsf,tftf]=meshgrid(SF_vec01,TF_vec01);
grid2.sfsf = sfsf;
grid2.tftf = tftf;

x = [];
x(:,1) = log2(grid2.sfsf(:));
x(:,2) = log2(grid2.tftf(:));

SFminmax = [log2(.01) log2(.64)];
TFminmax = [log2(.25) log2(32)];

figure
OLD = 1;
if OLD == 1
    h1 = figure;
    h2 = figure;
    ind_ALB = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==1  & Mask_rej2 );
    ind_PMB = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==2  & Mask_rej2 );
    ind_V1B = find(Big_ParVec_USErun(:,25)>0  & Big_ParVec_USErun(:,30)==0 & Big_ParVec_USErun(:,27)==3  & Mask_rej2 );

    NAL = length(ind_ALB);
    NV1 = length(ind_V1B);
    NPM = length(ind_PMB);
    Nmax = max([NAL,NV1,NPM])
    
    for count2 = 1:3; %1:3
        if count2 == 1
            ind_USE2 = ind_V1B;
            strtitle = 'V1';
            C = [0 0 0];
            N1 = NV1;
        elseif count2 == 2
            ind_USE2 = ind_ALB;
            strtitle = 'AL';
            C = [1 0 1];
            N1 = NAL;
        elseif count2 == 3
            ind_USE2 = ind_PMB;
            strtitle = 'PM';
            C = [0 1 0];
            N1 = NPM;
        end
        if ~isempty(who('h1'))
            close(h1);
        end
        
        if ~isempty(who('h2'))
            close(h2);
        end
        h1 = figure;
        h2 = figure;
        h1b = subplot(1,1,1);
        
        TFminmax0 =  log2([.5 24]);
        SFminmax0 =  log2([.02 .32]);
        
        TFminmax2 = [TFminmax0(1)-4 TFminmax0(2)+4];
        SFminmax2 = [SFminmax0(1)-4 SFminmax0(2)+4];
     
           TFminmax3 = [TFminmax0(1)-.25 TFminmax0(2)+.25];
        SFminmax3 = [SFminmax0(1)-.25 SFminmax0(2)+.25];
     
        for count1 = 1:length(ind_USE2)
%            xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) .606];
          %  xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) .1];
          %  xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) .606];
          %  xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) exp(-1*(1/4))];
            xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) exp(-1*(1/8))];
%            xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) .5];
%            xfit(4) = -4;
%            xfit(5) = 2;
            %ezcontour(@(x,y) Gauss2D_ellipseMA_forplotting(x,y,xfit));
            figure(h1);
            h0 = ezplot(@(x,y) Gauss2D_ellipseMA_forplotting(x,y,xfit),TFminmax2,SFminmax2);            
            hold on
            hold on            
            XData0 = get(h0,'XData');
            YData0 = get(h0,'YData');
%            hi = 1
%            pause(100)
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
            set(h3,'FaceColor',C,'FaceAlpha',.06*Nmax./N1,'EdgeColor',[1 1 1]*.4);
%            set(h3,'FaceColor',C,'FaceAlpha',1/N1,'EdgeColor',[1 1 1]*.4);
            axis([min(TFminmax0) max(TFminmax0) min(SFminmax0) max(SFminmax0)]);
            hold on
            end
           
        end
           line([min(TFminmax0) max(TFminmax0)],[min(SFminmax0) min(SFminmax0)],'Color','k')
            line([min(TFminmax0) max(TFminmax0)],[max(SFminmax0) max(SFminmax0)],'Color','k')
            line([min(TFminmax0) min(TFminmax0)],[min(SFminmax0) max(SFminmax0)],'Color','k')
            line([max(TFminmax0) max(TFminmax0)],[min(SFminmax0) max(SFminmax0)],'Color','k')
            
            axis([min(TFminmax3) max(TFminmax3) min(SFminmax3) max(SFminmax3)]);
          
        title(strtitle)

      figure(h1)
   file_print = [PWD_PRINT,'I2_ellipses_',strtitle,'_',str_run_PLOT4];
    print('-dtiff',[file_print,'.tif']); print('-depsc',[file_print,'.ps']);
      figure(h2)
   file_print = [PWD_PRINT,'I3_ellipses_',strtitle,'_',str_run_PLOT4];
    print('-dtiff',[file_print,'.tif']); print('-depsc',[file_print,'.ps']);

    end
 end


%now make a plot of the mean parameters
    figure
h000 =     subplot(1,1,1)
for count2 = 1:3
    if count2 == 1
        ind_USE2 = ind_V1;
        strtitle = 'V1';
        C = [0 0 0];
    elseif count2 == 2
        ind_USE2 = ind_AL;
        strtitle = 'AL';
        C = [0 0 .5];
    elseif count2 == 3
        ind_USE2 = ind_PM;
        strtitle = 'PM';
        C = [0 .5 0];
    end
    xfitvec = [];

    HalfHeight = .5;
    %HalfHeight = .606;
    for count1 = 1:length(ind_USE2)
        xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) HalfHeight];
%        xfit = [Big_ParVec_USErun(ind_USE2(count1),44:49) 0.05];
        xfit(4) = -4;
        xfit(5) = 2;
        xfitvec = [xfitvec; xfit];
    end
    
  %  ind = find(xfitvec(:,2)>1 & xfitvec(:,3)>1);
    ind = find(xfitvec(:,2)>0);
  %  xfit = median(xfitvec(ind,:),1);    
    xfit = mean(xfitvec(ind,:),1);    
 %   xfit = median(xfitvec(:,:),1);    
    %ezcontour(@(x,y) Gauss2D_ellipseMA_forplotting(x,y,xfit));
    h = ezplot(@(x,y) Gauss2D_ellipseMA_forplotting(x,y,xfit),TFminmax,SFminmax);
    set(h,'Color',C,'LineWidth',2)
    hold on
    axis([-1 5 -7 -1]);
    axis equal
end
    legend('V1','AL','PM');
       axis([-1 5 -7 -1]);
  xlabel('TF (octaves)');
            ylabel('SF (octaves)');
    XT = get(h000,'XTickLabel');
set(h000,'XTickLabel',num2str(str2num(deblank(XT)) - 2));
    YT = get(h000,'YTickLabel');
set(h000,'YTickLabel',num2str(str2num(deblank(YT)) - (-4)));
    title(['AvgContour at half-height (bandpass cells only): ' ,str_run_PLOT4]);
    file_print = [PWD_PRINT,'I_ellipses_','_',str_run_PLOT4];
    print('-dtiff',[file_print,'.tif']); print('-depsc',[file_print,'.ps']);
homeDir = 'Z:\All_Staff\home\miaomiao\labmeeting\151027';
xlsfile = fullfile(homeDir, 'retinotopic-caculation.xlsx');
[Gdata, Gtext, Graw] = xlsread(xlsfile,8);
ID = cell2mat(Graw(2:end,1));
ID_unique = unique(ID);
areas = {'V1','AL','LM','PM','RL','AM'};
Colors(1,:) = [0 0 0];
Colors(2:6,:)= lines(5);
X_val = cell2mat(Graw(2:end,3:2:13));
Y_val = cell2mat(Graw(2:end,4:2:14));
center = NaN(length(ID_unique),6,2);
long_axis = NaN(length(ID_unique),6);
short_axis = NaN(length(ID_unique),6);
region_area =  NaN(length(ID_unique),6); % in mm^2
distance = NaN(length(ID_unique),6,6);
ang = NaN(length(ID_unique),6);% the gangle of the ellipse
phi = NaN(length(ID_unique),6);
center_zero = NaN(length(ID_unique),6,2); % recenter via making V1 centered on [0 0]
ang_zero = NaN(length(ID_unique),6);
figure
for i = 1:length(ID_unique)
    idx = ID ==ID_unique(i);
    for i_area = 1:6
        %         scatter(X_val(idx,i_area),Y_val(idx,i_area),'MarkerEdgeColor',Colors(i_area,:))
        if sum(isnan(X_val(idx,i_area)))~=8 % deal with no such area identified
           
        
            ellipse_t = fit_ellipse(X_val(idx,i_area),Y_val(idx,i_area),1,Colors(i_area,:));%Colors(i_area,:)
            hold on
            center(i,i_area,1) = ellipse_t.X0_in;
            center(i,i_area,2) = ellipse_t.Y0_in;
            
            
            long_axis(i,i_area) = ellipse_t.long_axis;
            short_axis(i,i_area) = ellipse_t.short_axis;
            phi(i,i_area) = ellipse_t.phi;
            ang(i,i_area) = ellipse_t.ang;
            region_area(i,i_area) = ellipse_t.long_axis.*ellipse_t.short_axis.*pi./(4.*1000000); % mm^2
        end
        
    end
    
    for ii_area = 1:6
        
        for iii_area = 1:6
            if sum(isnan(X_val(idx,ii_area)))~=8 && sum(isnan(X_val(idx,iii_area)))~=8 
            % caculate distance
            distance(i,ii_area,iii_area) = pdist([squeeze(center(i,ii_area,:))';squeeze(center(i,iii_area,:))'],'euclidean');
            end
        end
        
        center_zero(i,ii_area,1) = squeeze(center(i,ii_area,1)) - squeeze(center(i,1,1));
        center_zero(i,ii_area,2) = squeeze(center(i,ii_area,2)) - squeeze(center(i,1,2));
        ang_zero(i,ii_area) = ang(i,ii_area) - ang(i,1);
        
        end 
    end
    
    

ylabel('um')
xlabel('um')


%% plot zeroed version
tau = 464; %1026; 
 figure
for i = 1%:length(ID_unique)
    for i_area = 1:6
        if ~isnan(short_axis(i,i_area))
        ellipse(short_axis(i,i_area)./2,long_axis(i,i_area)./2,ang(i,i_area),squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),'k')
        hold on
        
        end
        if i_area ==2
        ellipse(tau,tau,0,squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),'b:')
        end 
    end
end
%  xlim([-2000 2000])
% ylim([-2000 2500])
ylabel('um')
xlabel('um')
axis square
%% plot non zeroed version and store the images XY for each region,and a circular for each region for tau of spread:579.65 um
figure
DataXY = NaN(length(ID_unique),6,2,2001);
LEDXY = NaN(length(ID_unique),6,2,2001);
for i = 1:length(ID_unique)
    for i_area = 1:6
        if ~isnan(short_axis(i,i_area))
        h= ellipse(short_axis(i,i_area)./2,long_axis(i,i_area)./2,ang(i,i_area),squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),Colors(i_area,:),2000);
        DataXY(i,i_area,1,1:2001) = round(h.XData); % each value is a um
        DataXY(i,i_area,2,1:2001) = round(h.YData);
        
        h= ellipse(tau,tau,0,squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),'c',2000);
        LEDXY(i,i_area,1,1:2001) = round(h.XData);
        LEDXY(i,i_area,2,1:2001) = round(h.YData);
        
        hold on
        
        
        end
    end
end
ylabel('um')
xlabel('um')
%% adjust for big tau
% figure
% DataXY = NaN(length(ID_unique),6,2,10000);
% LEDXY = NaN(length(ID_unique),6,2,10000);
% for i = 1:length(ID_unique)
%     for i_area = 1:6
%         if ~isnan(short_axis(i,i_area))
%         h= ellipse(short_axis(i,i_area)./2,long_axis(i,i_area)./2,ang(i,i_area),squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),Colors(i_area,:),10000);
%         DataXY(i,i_area,1,1:10001) = round(h.XData); % each value is a um
%         DataXY(i,i_area,2,1:10001) = round(h.YData);
%         
%         h= ellipse(tau,tau,0,squeeze(center(i,i_area,1)),squeeze(center(i,i_area,2)),'c',10000);
%         LEDXY(i,i_area,1,1:10001) = round(h.XData);
%         LEDXY(i,i_area,2,1:10001) = round(h.YData);
%         
%         hold on
%         
%         
%         end
%     end
% end
% ylabel('um')
% xlabel('um')
% %% make it into pixel values in um, data sieze(5000*5000), needs to shift the x axis to avoid negative index
% value = 7000; 
% Data_pixel = NaN(length(ID_unique),6,value,value);
% LED_pixel = NaN(length(ID_unique),6,value,value);
% for  i = 1:length(ID_unique)
%     for i_area = 1:6
%         if ~isnan(short_axis(i,i_area))
%         a=[];       
%         a = sub2ind([value,value], squeeze(DataXY(i,i_area,2,1:10001))+700,squeeze(DataXY(i,i_area,1,1:10001))+700);
%         matrix = zeros(value,value);
%         mm = reshape(matrix, value*value,1);
%         mm(a) = 1;
%         BW = flipud(reshape(mm, value,value));
%         BW = im2bw(BW,0.5);                   %# binarize image
%         BW = imdilate(BW,strel('square',3)); %# dilation
%         BW = imfill(BW,'holes');             %# fill inside silhouette
%         BW = imerode(BW,strel('square',3));  %# erode
%         BW = bwperim(BW,8);                  %# get perimeter
%         BW2 = imfill(BW,'holes');
%         Data_pixel(i,i_area,1:value,1:value) = BW2;
%         
%         a=[];
%         a = sub2ind([value,value], squeeze(LEDXY(i,i_area,2,1:10001))+700,squeeze(LEDXY(i,i_area,1,1:10001))+700);
%         matrix = zeros(value,value);
%         mm = reshape(matrix, value*value,1);
%         mm(a) = 1;
%         BW = flipud(reshape(mm, value,value));
%         BW = im2bw(BW,0.5);                   %# binarize image
%         BW = imdilate(BW,strel('square',3)); %# dilation
%         BW = imfill(BW,'holes');             %# fill inside silhouette
%         BW = imerode(BW,strel('square',3));  %# erode
%         BW = bwperim(BW,8);                  %# get perimeter
%         BW2 = imfill(BW,'holes');
%         LED_pixel(i,i_area,1:value,1:value) = BW2;
%         end 
%     end
% end
%for smaller tau
Data_pixel = NaN(length(ID_unique),6,5000,5000);
LED_pixel = NaN(length(ID_unique),6,5000,5000);
for  i = 1:length(ID_unique)
    for i_area = 1:6
        if ~isnan(short_axis(i,i_area))
        a=[];       
        a = sub2ind([5000,5000], squeeze(DataXY(i,i_area,2,1:2001)),squeeze(DataXY(i,i_area,1,1:2001))+500);
        matrix = zeros(5000,5000);
        mm = reshape(matrix, 5000*5000,1);
        mm(a) = 1;
        BW = flipud(reshape(mm, 5000,5000));
        BW = im2bw(BW,0.5);                   %# binarize image
        BW = imdilate(BW,strel('square',3)); %# dilation
        BW = imfill(BW,'holes');             %# fill inside silhouette
        BW = imerode(BW,strel('square',3));  %# erode
        BW = bwperim(BW,8);                  %# get perimeter
        BW2 = imfill(BW,'holes');
        Data_pixel(i,i_area,1:5000,1:5000) = BW2;
        
        a=[];
        a = sub2ind([5000,5000], squeeze(LEDXY(i,i_area,2,1:2001)),squeeze(LEDXY(i,i_area,1,1:2001))+500);
        matrix = zeros(5000,5000);
        mm = reshape(matrix, 5000*5000,1);
        mm(a) = 1;
        BW = flipud(reshape(mm, 5000,5000));
        BW = im2bw(BW,0.5);                   %# binarize image
        BW = imdilate(BW,strel('square',3)); %# dilation
        BW = imfill(BW,'holes');             %# fill inside silhouette
        BW = imerode(BW,strel('square',3));  %# erode
        BW = bwperim(BW,8);                  %# get perimeter
        BW2 = imfill(BW,'holes');
        LED_pixel(i,i_area,1:5000,1:5000) = BW2;
        end 
    end
end
%% get the ovelayed pixel between areas
P_cover = NaN(length(ID_unique),6,6);
for i= 1:length(ID_unique)
    for i_led = 1:6
        if ~isnan(short_axis(i,i_led))
        for i_area = 1:6
            if ~isnan(short_axis(i,i_area))
                LED_temp =  squeeze(LED_pixel(i,i_led,:,:)); 
                Area_temp = squeeze(Data_pixel(i,i_area,:,:)); 
                cover = [];
                cover = LED_temp&Area_temp;
                T_area = sum(Area_temp(:));
                T_cover = sum(cover(:));
                P_cover(i,i_led,i_area) = T_cover./T_area; 
            end 
        end 
        end
    end 
end 
%% plot percentage coverage
m_pcover = squeeze(nanmean(P_cover,1)).*100;
std_pcover = squeeze(nanstd(P_cover,[],1)).*100;
n_pcover = squeeze(sum(~isnan(P_cover),1));
sem_pcover = std_pcover./sqrt(n_pcover);
figure
subplot(3,2,1)
imagesc(m_pcover,[0 100])
colormap(brewermap([],'Greys'))
colorbar
set(gca,'TickDir','out','XTick',1:1:6,'XTickLabel',{'V1','AL','LM','PM','RL','AM'},'YTick',1:1:6,'YTickLabel',{'V1','AL','LM','PM','RL','AM'})
ylabel('LED center area')
xlabel('Areas')
for i=1:4
subplot(3,2,i+2)
h=errorbar(1:1:6,m_pcover(i,:),sem_pcover(i,:),'k');
h.Marker = 'o';
h.LineStyle = 'none';
set(gca,'TickDir','out','XTick',1:1:6,'XTickLabel',{'V1','AL','LM','PM','RL','AM'},'YTick',0:50:100)
ylabel('% of area coverage')
xlabel('Areas')
ylim([-10 110])
xlim([0 7])
title(['LED center area:' areas{i} ])

end 

% plot the mean of percent coverage
%% bin the imagesc?
[c,idx]=histc(m_pcover(:),[0 0.1 2 50 98 101]);
bin_pcover = reshape(idx,6,6);

figure
imagesc(bin_pcover,[1 5])
colormap(brewermap(5,'Greys'))

colorbar
set(gca,'TickDir','out','XTick',1:1:6,'XTickLabel',{'V1','AL','LM','PM','RL','AM'},'YTick',1:1:6,'YTickLabel',{'V1','AL','LM','PM','RL','AM'})

ylabel('LED center area')
xlabel('Areas')
axis square

%% plot the distance across areas, sort from the shortest to longest
% deal with NaN values for std
sizedata = 60;
dist_mean = squeeze(nanmean(distance,1))./1000;
dist_std = squeeze(nanstd(distance,[],1))./1000;


figure
errorbar(1,dist_mean(6,4),dist_std(6,4)./sqrt(size(short_axis,1)),'k')
hold on
scatter(1,dist_mean(6,4),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)



errorbar(2,dist_mean(3,2),dist_std(3,2)./sqrt(size(short_axis,1)),'k')
hold on
scatter(2,dist_mean(3,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)


errorbar(3,dist_mean(5,2),dist_std(5,2)./sqrt(sum(~isnan(short_axis(:,5)))),'k')
hold on
scatter(3,dist_mean(5,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)


errorbar(4,dist_mean(3,1),dist_std(3,1)./sqrt(size(short_axis,1)),'k')
hold on
scatter(4,dist_mean(3,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)

errorbar(5,dist_mean(2,1),dist_std(2,1)./sqrt(size(short_axis,1)),'k')
hold on
scatter(5,dist_mean(2,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)


errorbar(6,dist_mean(4,1),dist_std(4,1)./sqrt(size(short_axis,1)),'k')
hold on
scatter(6,dist_mean(4,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)

errorbar(7,dist_mean(4,2),dist_std(4,2)./sqrt(size(short_axis,1)),'k')
hold on
scatter(7,dist_mean(4,2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)

errorbar(8,dist_mean(4,3),dist_std(4,3)./sqrt(size(short_axis,1)),'k')
hold on
scatter(8,dist_mean(4,3),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',sizedata)
hline(0.57965,'k:')
ylim([0 2.8])
xlim([0.5 8.5])
set(gca,'XTick',1:1:8,'YTick',0:0.4:2.8,'XTicklabel',{'AM-PM','AL-LM','AL-RL','LM-V1','AL-V1','PM-V1','PM-AL','PM-LM'},'TickDir','out')
ylabel('Distance between areas(mm)')

%% plot the areas across regions
figure
h=errorbar (1:1:4,mean(region_area,1),std(region_area,[],1)./sqrt(size(region_area,1)),'k');
h.LineStyle = 'none';
hold on
scatter(1:1:4,mean(region_area,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1])
set(gca,'XTick',1:4,'XTicklabel',{'V1','AL','LM','PM'},'TickDir','out')
ylim([0 0.7])
ylabel('area (mm^2)')

%% plot overlayed all ellipse
figure
for i = 1:length(ID_unique)
    for i_area = 1:6
        if ~isnan(short_axis(i,i_area))
        ellipse(short_axis(i,i_area)./2000,long_axis(i,i_area)./2000,0,i_area-1,0,Colors(i_area,:))
        end
        hold on
        text(i_area-1.08, 0,areas{i_area},'Color',Colors(i_area,:),'FontSize',8)
    end
end
%  xlim([-2000 2000])
 xlim([-0.8 5.6])
set(gca,'YTick',-0.8:0.4:0.8,'XTick',-0.8:0.4:5.6,'TickDir','out')
ylabel('mm')
xlabel('mm')

%% plot the schematic for showing the oval fitting

images=readtiff('Z:\All_Staff\home\miaomiao\labmeeting\170907\170203_i558_el[10]az[10 40]-1.tif');
scale_pixel =0.0655; 
 figure
 subplot(1,2,1)
 imagesc(images)
 colormap gray
 
 subplot(1,2,2)

    for i_area = 1:5
        if ~isnan(short_axis(9,i_area))
        ellipse((short_axis(9,i_area)./2)*scale_pixel,(long_axis(9,i_area)./2)*scale_pixel,ang(9,i_area),squeeze(center(9,i_area,1))*scale_pixel,squeeze(center(9,i_area,2))*scale_pixel,'k')
        hold on
        scatter(squeeze(center(9,i_area,1))*scale_pixel,squeeze(center(9,i_area,2))*scale_pixel,'.','MarkerEdgeColor',[0 0 0])
        
        end
        
    end
    
 xlim([1 251])
 ylim([1 250])
ylabel('um')
xlabel('um')

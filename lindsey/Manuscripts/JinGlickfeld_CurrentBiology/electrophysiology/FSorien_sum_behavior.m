% only for behavior paper 
areas={'V1','LM','AL','PM'};
colors = [0 0 0; lines(5)];
edges = -0.2:0.020:0.35;
color_con = [0 0 0;0.2 0.6 0.4; 1 0 0];
Ix_window = edges>0 & edges <= 0.13;
base_bin = find(edges==0);
% change everything to mean response
for i=1:length(areas)
    Gdata=[];
    Gtext=[];
    Graw = [];
    [Gdata, Gtext, Graw] = xlsread('Z:\All_Staff\home\miaomiao\Analysis\electrode\ducumentation\FSOrien.xlsx',i);
    all_cells(i,1)=0;
  
    
    % get response and adaptation index for responsive cells 
    for i_con = 1:3
    All.response{i,i_con} = []; % 
    All.response_sem{i,i_con} = [];  
    All.CirV{i,i_con} = []; % circular variance
    All.pref_ori{i,i_con} = []; 
    
    % for cells responsive to 0 deg
    Select.response{i,i_con} = []; % 
    Select.response_sem{i,i_con} = [];  
    
    Select.trace{i,i_con} = []; 
    Select.area{i,i_con} = []; 
    Select.areaN{i,i_con} = []; 
   
    end 
    Select.latency{i,1} = []; 
    All.latency{i,1} = []; % store the latency
    All.SNR{i,1} = []; % for all first
    All.ori_idx{i,1} = []; % cells that are reliable fit in the control
    All.ori_tuned{i,1} = []; % cells significantly modulated by orientations p values
    All.ori_tuned_idx{i,1} = []; % index
    First.pref{i,1} = [];
    First.OSI{i,1} = [];
    All.corr_idx{i,1} = []; % cells positive corr btw first and adpt 
    All.corr_idx{i,2} = []; % corr btw first and test
    All.FA_adpt_idx{i,1} = []; % store all adaptation index values   
    All.FT_adpt_idx{i,1} = []; 
    
    for i_exp = 1:size(Gdata,1)
        
        load(cell2mat(Gtext(i_exp+1,4)))
        all_cells(i,1) = all_cells(i,1)+  length(FS.responseIX); 
        %for responsive cells and reliable fit in the first condition only
        % select cell either response to collapsed first, or either of the
        % orientation, use Bonferroni correction
        p_thresh = 0.05/7; 
        idx = [];
        idx = FS.stats_areaFP<p_thresh|(nansum(FS.stats_oriPN<p_thresh,2)>0)';
              
   
        All.response{i,1} = [All.response{i,1}; FS.response.first.mean(idx,:)]; % get all the response
        All.response_sem{i,1} = [All.response_sem{i,1}; FS.response.first.sem(idx,:)];
        All.response{i,2} = [All.response{i,2}; FS.response.adpt.mean(idx,:)]; % get all the response
        All.response_sem{i,2} = [All.response_sem{i,2}; FS.response.adpt.sem(idx,:)];
        All.response{i,3} = [All.response{i,3}; FS.response.target.mean(idx,:)]; % get all the response
        All.response_sem{i,3} = [All.response_sem{i,3}; FS.response.target.sem(idx,:)];
        
        
        
      
       
        % get index for cells significantly modulated by orien, using
        % kruskalwallis test
        id_template = [];
        id_template = find(idx==1);
        temp_p = []; 
        temp_snr = []; 
        for i_template = 1:sum(idx)             
            
            num_trial = cellfun(@length,FS.response.first.raw(id_template(i_template),:));
            ori_idx =[];
            [~,ori_idx]=max(FS.response.first.mean(id_template(i_template),:)); 
            temp_snr(i_template,1) = FS.response.first.mean(id_template(i_template),ori_idx)./(FS.response.first.sem(id_template(i_template),ori_idx)*sqrt(num_trial(ori_idx)));
            
            data = [];
            for i_ori=1:6                     
                data = [data; FS.response.first.raw{id_template(i_template),i_ori}];    
            end
            group = [ones(num_trial(1),1); ones(num_trial(2),1)+1;ones(num_trial(3),1)+2; ones(num_trial(4),1)+3 ; ones(num_trial(5),1)+4; ones(num_trial(6),1)+5];
            
            
            %[temp_p(i_template,1),~,~] = kruskalwallis(data,group);
            [temp_p(i_template,1),~,~] = anova1(data,group);
            close all
            
        end
        % get the SNR
        All.SNR{i,1} = [All.SNR{i,1}; temp_snr]; 
        All.ori_tuned{i,1}=[All.ori_tuned{i,1}; temp_p]; 
        All.ori_tuned_idx{i,1} = [All.ori_tuned_idx{i,1}; temp_p<0.05]; 
        All.latency{i,1} = [All.latency{i,1} FS.Latency(idx)]; 
        % calculate vector strenth
        temp = [];
        temp = FS.response.first.mean(idx,:); 
        temp = temp-repmat(min(temp,[],2),1,6); 
        
        Vector = [];
        Vector = dot(temp,repmat(exp((2*pi*sqrt(-1)*FS.BOrien)./180),size(temp,1),1),2)./(sum(temp,2));
        All.CirV{i,1} = [All.CirV{i,1}; 1-abs(Vector)]; 
        temp_ori = [];
        temp_ori = rad2deg(angle(Vector));
        temp_ori (temp_ori<0) = 360+temp_ori(temp_ori<0); 
        temp_ori = temp_ori/2; 
        All.pref_ori{i,1} = [ All.pref_ori{i,1};temp_ori]; 
        
        temp = [];
        temp = FS.response.adpt.mean(idx,:); 
        temp = temp-repmat(min(temp,[],2),1,6); 
        
        Vector = [];
        Vector = dot(temp,repmat(exp((2*pi*sqrt(-1)*FS.BOrien)./180),size(temp,1),1),2)./(sum(temp,2));
        All.CirV{i,2} = [All.CirV{i,2}; 1-abs(Vector)]; 
        temp_ori = [];
        temp_ori = rad2deg(angle(Vector));
        temp_ori (temp_ori<0) = 360+temp_ori(temp_ori<0); 
        temp_ori = temp_ori/2; 
        All.pref_ori{i,2} = [ All.pref_ori{i,2};temp_ori]; 
        
        temp = [];
        temp = FS.response.target.mean(idx,:); 
        temp = temp-repmat(min(temp,[],2),1,6); 
        
        Vector = [];
        Vector = dot(temp,repmat(exp((2*pi*sqrt(-1)*FS.BOrien)./180),size(temp,1),1),2)./(sum(temp,2));
        All.CirV{i,3} = [All.CirV{i,3}; 1-abs(Vector)]; 
        temp_ori = [];
        temp_ori = rad2deg(angle(Vector));
        temp_ori (temp_ori<0) = 360+temp_ori(temp_ori<0); 
        temp_ori = temp_ori/2; 
        All.pref_ori{i,3} = [ All.pref_ori{i,3};temp_ori]; 
        
        % for cells response to 0 deg orientation
            
        idx=[];
        FS.stats_oriN(isnan(FS.stats_oriN(:,1)),1) = 0;
        idx =logical(FS.stats_oriN(:,1));
        % get the response trace
        Select.trace{i,1} = [Select.trace{i,1}; squeeze(FS.Traces.first.raw_mean(idx,1,:))]; 
        Select.trace{i,2} = [Select.trace{i,2}; squeeze(FS.Traces.adpt.raw_mean(idx,1,:))]; 
        Select.trace{i,3} = [Select.trace{i,3}; squeeze(FS.Traces.target.raw_mean(idx,1,:))]; 
        
        Select.latency{i,1} = [Select.latency{i,1} FS.Latency(idx)];
        
        Select.response{i,1} = [Select.response{i,1}; FS.response.first.mean(idx,1)]; %
        Select.response_sem{i,1} = [Select.response_sem{i,1}; FS.response.first.sem(idx,1)];
        Select.response{i,2} = [Select.response{i,2}; FS.response.adpt.mean(idx,1)]; %
        Select.response_sem{i,2} = [Select.response_sem{i,2}; FS.response.adpt.sem(idx,1)];
        Select.response{i,3} = [Select.response{i,3}; FS.response.target.mean(idx,1)]; %
        Select.response_sem{i,3} = [Select.response_sem{i,3}; FS.response.target.sem(idx,1)];
        % get the mean response within 0-130 ms window
        tmp=[];
        tmp = squeeze(FS.Traces.first.raw_mean(idx,1,:)); % get all the trace first
        Select.area{i,1} =[Select.area{i,1}; mean(tmp(:,Ix_window),2)-tmp(:,base_bin)]; 
        tmp=[];
        tmp = squeeze(FS.Traces.adpt.raw_mean(idx,1,:)); % get all the trace first
        Select.area{i,2} =[Select.area{i,2}; mean(tmp(:,Ix_window),2)-tmp(:,base_bin)]; 
        tmp=[];
        tmp = squeeze(FS.Traces.target.raw_mean(idx,1,:)); % get all the trace first
        Select.area{i,3} =[Select.area{i,3}; mean(tmp(:,Ix_window),2)-tmp(:,base_bin)]; 
        % get the new area response, peak+-2 bins (100 ms)
        Select.areaN{i,1} = [Select.areaN{i,1}; FS.response.first.area_mean(idx,1)]; 
        Select.areaN{i,2} = [Select.areaN{i,2}; FS.response.adpt.area_mean(idx,1)]; 
        % get the OSI measurea
        idx = [];
        idx = FS.stats_Ffit&FS.fit.first.R_square>=0.5;  % reliable fit cells
        
        First.pref{i,1} = [First.pref{i,1}  FS.fit.first.prefer(idx)];
        
        First.OSI{i,1} = [First.OSI{i,1} FS.fit.first.OSI(find(idx==1))];
        
        
        
     
    
    end
end 
%%  chisquare test for number of responsive cells
Num_resp = [cellfun(@length,All.SNR) all_cells-cellfun(@length,All.SNR)]; 


 x=Num_resp;
 
 e = sum(x,2)*sum(x)/sum(x(:)); % expected
 X2 = (x-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(x)-[1 1]); % degree of freedom
 p = 1-chi2cdf(X2,df);
 
%% for multiple comparison
p_multi = zeros(4,4); 
for i=1:4
    for ii=1:4
        if ii>i
            [~,p_multi(i,ii)] = chi2cont(Num_resp([i ii],:)');
          
        end 
    end
end


%% percentage of responsive cells
figure
p_responsive = cellfun(@length,All.SNR)./all_cells; 
p_tuned = cellfun(@sum,All.ori_tuned)./cellfun(@length,All.SNR);
% chi square test across areas
subplot(2,2,1)
for i=1:4
    bar(i,p_responsive(i),'FaceColor',colors(i,:),'EdgeColor',colors(i,:))
    hold on
end
ylim([0 1])
xlim ([0.5 4.5])
ylabel('% of visual driven cells')
set(gca,'TickDir','Out','XTick',1:4, 'XTickLabel',areas)
axis square

subplot(2,2,2)
% plot SNR
for i_area = 1:4
    [f,x]= ecdf(All.SNR{i_area,1});
    plot(x,f,'color',colors(i_area,:))
    hold on
end
 xlim([0 3.5])
xlabel('SNR: Mean/std')
ylabel ('Fraction of cells')
set(gca,'TickDir','out')
title('For preferred ori')
axis square


subplot(2,2,3)

for i=1:4
    bar(i,p_tuned(i),'FaceColor',colors(i,:),'EdgeColor',colors(i,:))
    hold on
end
ylim([0 1])
xlim ([0.5 4.5])
ylabel('% of ori tuned cells')
set(gca,'TickDir','Out','XTick',1:4, 'XTickLabel',areas)
axis square
title('using one way anova')

subplot(2,2,4)
for i_area = 1:4
    [f,x]= ecdf(All.CirV{i_area,1});
    plot(x,f,'color',colors(i_area,:))
    hold on
end
xlabel('Circular variance')
ylabel ('Fraction of cells')
set(gca,'TickDir','out')
axis square
%% stats of one way anova

num_cell = cellfun(@length, All.CirV(:,1));

data = [];
for i_area=1:4
    data = [data;  All.SNR{i_area,1}];
end
group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3];


%[temp_p(i_template,1),~,~] = kruskalwallis(data,group);
[temp_p,tbl,stats] = anova1(data,group);
p = multcompare(stats);
%% plot average respons trace
binwidth = FS.binwidth; 
edge_idx = edges<=0.3&edges>=-0.1;
figure
for i=1:4
    for j=1:2
        tmp_mean = mean(Select.trace{i,j},1);
        tmp_sem = std(Select.trace{i,j},[],1)./sqrt(size(Select.trace{i,j},1));
        
        shadedErrorBar_ch(edges(edge_idx)+(j-1)*0.5,(tmp_mean(edge_idx)./binwidth)-15*(i-1),tmp_sem(edge_idx)./binwidth,{'color',colors(i,:),'linewidth',1},0);
        hold on
        
    end
    
end
vline([0 0.1 0+0.5 0.1+0.5],'k:')

%% distribution of orientation preference
for i_area =1:4
histogram(First.pref{i_area,1})
hold on
end 
figure
for i_area =1:4
histogram(All.pref_ori{i_area,1})
hold on
end 
%% plot 0deg responsive cells
%response_FR = cellfun(@(x) x./FS.binwidth,  Select.response, 'UniformOutput', 0); 

response_FR = cellfun(@(x) x./FS.binwidth,  Select.areaN, 'UniformOutput', 0); 
n_cell = cellfun(@length, response_FR(:,1)); 
figure
% response to first stimulus
subplot(1,2,1)
plotSpread(response_FR(:,1), 'distributionColors', colors(1:4,:));
% ylim([0 140])
hold on
%errorbar(1:4, cellfun(@mean, response_FR(:,1)), cellfun(@std,response_FR(:,1))./sqrt(n_cell),'LineStyle','none','Marker','o')
% plot median and 25th and 75th quantile
tmp=[];
tmp = response_FR(:,1); 
tmp_median = cellfun(@median,tmp); 
tmp_quantile = cell2mat(cellfun(@(X) quantile(X,[0.25 0.75]), tmp,'UniformOutput',false)); 
for i_area =1:4
    errorbar(i_area,tmp_median(i_area),tmp_quantile(i_area,1),tmp_quantile(i_area,2),'Marker','o','color',colors(i_area,:))
    hold on
end 
 ylim([0 60])
set(gca,'XTick',1:4,'XTickLabel', areas)
axis square
ylabel('Response (Hz)')
title('Stim_1_s_t')
subplot(1,2,2)
plotSpread(response_FR(:,2), 'distributionColors', colors(1:4,:));

hold on
%errorbar(1:4, cellfun(@mean, response_FR(:,2)), cellfun(@std,response_FR(:,2))./sqrt(n_cell),'LineStyle','none','Marker','o')
% plot median and 25th and 75th quantile
tmp=[];
tmp = response_FR(:,2); 
tmp_median = cellfun(@median,tmp); 
tmp_quantile = cell2mat(cellfun(@(X) quantile(X,[0.25 0.75]), tmp,'UniformOutput',false)); 
for i_area =1:4
    errorbar(i_area,tmp_median(i_area),tmp_quantile(i_area,1),tmp_quantile(i_area,2),'Marker','o','color',colors(i_area,:))
    hold on
end 
 ylim([0 60])
set(gca,'XTick',1:4,'XTickLabel', areas)
axis square
ylabel('Response (Hz)')
title('Stim_4_-_5_t_h')


%% anova test across areas

num_cell = cellfun(@length,response_FR(:,1));

data = [];
for i=1:4
    
    %data = [data;  Select.latency{i}'];
         data = [data;  response_FR{i,2}];
    %data = [data;All.FT_adpt_idx{i,1}- All.FA_adpt_idx{i,1}];
end
group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3; ones(num_cell(5),1)+4; ones(num_cell(6),1)+5 ];

[p1,tbl,stats] = kruskalwallis(data,group);
b= multcompare(stats,'CType','dunn-sidak');
%%
[p,h] = ranksum(All.FT_adpt_idx{5,1},0)
%% measure the difference between adaptation index for cells that are adapted
New_adapt = {}; 
for i = 1:5
    New_adapt{i,1} = All.FT_adpt_idx{i,1}(All.FA_adpt_idx{i,1}<0) - All.FA_adpt_idx{i,1}(All.FA_adpt_idx{i,1}<0); 
    
end 
%% stats

num_cell = cellfun(@length,New_adapt(1:5,1));

data = [];
for i=1:5
    
    %data = [data; Adpt_idx_new{1,i}];
         data = [data; New_adapt{i,1}];
    %data = [data;All.FT_adpt_idx{i,1}- All.FA_adpt_idx{i,1}];
end
group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ; ones(num_cell(5),1)+4];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3; ones(num_cell(5),1)+4; ones(num_cell(6),1)+5 ];

[p1,tbl,stats] = kruskalwallis(data,group);
b= multcompare(stats,'CType','dunn-sidak');

%% plot how adaptation changes the circular variance and preferred orientation 

figure
for i_area = 1:5
    subplot(2,3,i_area)
    Temp = [];
    Temp = All.pref_ori{i_area,2}-All.pref_ori{i_area,1};
    Temp(Temp>90) =180-Temp(Temp>90);
    Temp(Temp<-90) =Temp(Temp<-90)+180;
    FA_ori{i_area,1} = Temp;

    scatter(All.pref_ori{i_area,1},Temp, 20,'MarkerEdgeColor',color_con(2,:))
    hold on
    Temp = [];
    Temp = All.pref_ori{i_area,3}-All.pref_ori{i_area,1};
    Temp(Temp>90) =180-Temp(Temp>90);
    Temp(Temp<-90) =Temp(Temp<-90)+180;
    FT_ori{i_area,1} = Temp;
    scatter(All.pref_ori{i_area,1},Temp, 20,'MarkerEdgeColor',color_con(3,:))
    xlabel('Pref ori-control')
    ylabel('Ori diff')
end 


%% plot the mean regression lines??

%% stats
[p,h]=signrank(All.Aintercept{i_area,1}(All.FA_corr{i_area,1}>0.5),0); 
[p,h] = ranksum(All.Tintercept{i_area,1}(All.FT_corr{i_area,1}>0.5),All.Aintercept{i_area,1}(All.FA_corr{i_area,1}>0.5));

[p,h]=signrank(All.FT_adpt_idx{5,1}, 0); 

%% get the FS dataset for behavior
areas={'V1','LM','AL','PM'};
colors = [0 0 0; lines(3)];
edges = -0.2:0.020:0.35;

for i=1:length(areas)
    Gdata=[];
    Gtext=[];
    Graw = [];
    [Gdata, Gtext, Graw] = xlsread('Z:\All_Staff\home\miaomiao\Analysis\electrode\ducumentation\FS.xlsx',i);
    ii_cell = 1;
    All_cell{i,1} = 0;
    layers{i,1} = [];
    response{i,1} = []; % normalized to collapsed first
    response{i,2} = []; % 500 ms off
    resp_Hz{i,1} = []; % unnormalized ones
    resp_Hz{i,2} = [];
    spon{i,1} = []; % already in Hz
    
    Date{i,1} = [];
    ID{i,1} = [];
    genotype{i,1} = [];
    full{i,1} = [];
    naive{i,1} = []; 
    
    First{i,1} = []; % first response
    
    for i_exp = 1:size(Gdata,1)
        
        load(cell2mat(Gtext(i_exp+1,2)))
        All_cell{i,1} = All_cell{i,1}+length(FS.responseIX); 
        
%         load(cell2mat(Gtext(i_exp+1,5)))
        idx=[];
        idx = FS.BS_IX; % for cells response to first pulse
           
        % get the mouse information
        Date{i,1} = [Date{i,1} repmat(Gdata(i_exp,1),1,length(idx))]; % each cell has a date tag
        ID{i,1} = [ID{i,1} repmat(Gtext(i_exp+1,6),1,length(idx))];
        genotype{i,1}= [genotype{i,1} repmat(Gtext(i_exp+1,9),1,length(idx))];
        full{i,1} = [full{i,1} repmat(Gdata(i_exp,7),1,length(idx))]; % 1 is full screen condition
        naive{i,1} = [naive{i,1} repmat(Gdata(i_exp,3),1,length(idx))]; % 1 is naive
        
        First{i,1} = [First{i,1}; FS.response.first.mean(idx)./FS.binwidth];
        % 250/500ms adapted 4-5 average
        resp_Hz{i,1} = [resp_Hz{i,1}; squeeze(mean(FS.response.mean(idx,1,4:5),3))./FS.binwidth]; 
        resp_Hz{i,2} = [resp_Hz{i,2}; squeeze(mean(FS.response.mean(idx,2,4:5),3))./FS.binwidth];    
        
        
        
    end
    
    
end

%% plot the response
n_cell = cellfun(@length, First); 
figure
% response to first stimulus
subplot(1,2,1)
plotSpread(First, 'distributionColors', colors(1:4,:));
 ylim([0 150])
 ylabel('Response (Hz)')
hold on
errorbar(1:4, cellfun(@mean, First), cellfun(@std,First)./sqrt(n_cell),'LineStyle','none','Marker','o')
set(gca,'XTick',1:4,'XTickLabel', areas)
title('Stim_1_s_t')
axis square
subplot(1,2,2)
plotSpread(resp_Hz(:,1), 'distributionColors', colors(1:4,:));
 ylim([0 150])
hold on
errorbar(1:4, cellfun(@mean, resp_Hz(:,1)), cellfun(@std,resp_Hz(:,1))./sqrt(n_cell),'LineStyle','none','Marker','o')
set(gca,'XTick',1:4,'XTickLabel', areas)
axis square
ylabel('Response (Hz)')
title('Stim_4_-_5_t_h')


%%
num_cell = cellfun(@length,resp_Hz(:,1));

data = [];
for i=1:4
    
    %data = [data; Adpt_idx_new{1,i}];
         data = [data; resp_Hz{i,1}];
    %data = [data;All.FT_adpt_idx{i,1}- All.FA_adpt_idx{i,1}];
end
group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3 ];
%group = [ones(num_cell(1),1); ones(num_cell(2),1)+1;ones(num_cell(3),1)+2; ones(num_cell(4),1)+3; ones(num_cell(5),1)+4; ones(num_cell(6),1)+5 ];

[p1,tbl,stats] = kruskalwallis(data,group);
b= multcompare(stats,'CType','dunn-sidak');


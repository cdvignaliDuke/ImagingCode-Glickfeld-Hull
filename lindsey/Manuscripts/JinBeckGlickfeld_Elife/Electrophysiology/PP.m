
%% new pair pulse analysis
function  [PP]=PP(input,Trigger,spikes,LFP)
%%
global Fs stimOnTimeS Offs
Fs=double(Trigger.SamplingFreq);

[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI(Trigger,input);
%[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI_nolabjack(Trigger,input);
% use the LFP identify the layers of each template
baseline=0.2;
binwidth = 0.020;

 [Data,psth,Data_full,psth_full]=align_isi(spikes,baseline,binwidth,Time_start_new);

%% now separate by different off conditions of the psth_full 
StimOffs = double(cell2mat(input.tStimOffTimeMs));
Offs = unique(StimOffs);
edges_full = -baseline:binwidth:5;
edges = -baseline:binwidth:0.35;
base = find(edges==0); % use this to subtract the peak response
basebins = round(baseline./binwidth);
Ix_window = edges>=0 & edges <= 0.12; % for the timewindow 120ms after visual start
on_edges =  edges(Ix_window);

for i_template = 1:length(psth.area)
    
    First = [];
    First = cell2mat(cellfun(@(x) x(1,:),psth.raw{i_template,1},'UniformOutput',false));
    PP.Traces.raw_mean(i_template,6,:) = mean(First,1);
    PP.Traces.raw_sem(i_template,6,:) = std(First,[],1)./sqrt(size(First,1));
    
    % find the peak of the First, apply it to all the pulses
    [~, ~,PP.responseIX(i_template,1)]=findpeakmiao(mean(First,1),Ix_window,on_edges,basebins);
    
    PP.response.raw{i_template,6} = First(:,PP.responseIX(i_template,1))-mean(First(:,base)); % already subtracted by bin 0
    PP.response.mean(i_template,6) = mean(First(:,PP.responseIX(i_template,1))-mean(First(:,base)));
    PP.response.sem(i_template,6) = std(First(:,PP.responseIX(i_template,1))-mean(First(:,base)))./sqrt(size(First,1));
    
    for i_off = 1:length(Offs)
    % try randomly select subset of trials
%     temp = psth_full.raw{i_template,1}(StimOffs==Offs(i_off),:);
%     y = datasample(temp,28,1);
    
    PP.full.Trace_mean(i_template,i_off,:) = mean( psth_full.raw{i_template,1}(StimOffs==Offs(i_off),:),1);
    PP.full.Trace_sem(i_template,i_off,:) = std( psth_full.raw{i_template,1}(StimOffs==Offs(i_off),:),[],1)./sqrt(sum(StimOffs==Offs(i_off)));
    
%     PP.full.Trace_mean(i_template,i_off,:) = mean( y,1);
%     PP.full.Trace_sem(i_template,i_off,:) = std(y,[],1)./sqrt(28);
    

    Second = [];
    Second = cell2mat(cellfun(@(x) x(2,:),psth.raw{i_template,1},'UniformOutput',false));
    Second = Second(StimOffs==Offs(i_off),:);
    
    PP.Traces.raw_mean(i_template,i_off,:) = mean(Second,1);
    PP.Traces.raw_sem(i_template,i_off,:) = std(Second,[],1)./sqrt(sum(StimOffs==Offs(i_off)));
    
  
    
    PP.response.raw{i_template,i_off} = Second(:,PP.responseIX(i_template,1))-mean(Second(:,base)); % already subtracted by bin 0
    
    PP.response.normalize.raw{i_template,i_off} = PP.response.raw{i_template,i_off}./PP.response.mean(i_template,6);
    PP.response.normalize.mean(i_template,i_off) = mean( PP.response.normalize.raw{i_template,i_off});
    PP.response.normalize.sem(i_template,i_off) = std(PP.response.normalize.raw{i_template,i_off})./sqrt(sum(StimOffs==Offs(i_off)));
    
    
    PP.response.mean(i_template,i_off) = mean(Second(:,PP.responseIX(i_template,1))-mean(Second(:,base)));
    PP.response.sem(i_template,i_off) = std(Second(:,PP.responseIX(i_template,1))-mean(Second(:,base)))./sqrt(sum(StimOffs==Offs(i_off)));
    
    
    end
    
   
    
end
% only include the cell that has significant response to the first pulse

PP.cellIX = []; % cells that responde to the first pulse for every conditions pulse one last pulse
PP.cellanova = [];% cells that responde to collapsed first pulse and also no difference between conditions using annova
PP.anovastrict = []; % cells that responde to collapsed first pulse and also no difference between conditions using annova, and must has response to the second pulse after 4s
% PP.cellttest = []; % cells that responde to the collapsed first pulse and also no different between conditions using ttest
% PP.cellttestloose = []; 
h=[];
p=[];


i=1;
j=1;
jj=1;
% n=1;
% ii=1;
for i_template = 1:length(psth.area)
    
    H = ttest2(psth.First{i_template},psth.baseline{i_template},'Tail','right'); % the mean of X is greater than the mean of Y; 
    
    response = NaN(max([sum(StimOffs==Offs(1)),sum(StimOffs==Offs(2)),sum(StimOffs==Offs(3)),sum(StimOffs==Offs(4)),sum(StimOffs==Offs(5))]),length(Offs));
    
    for  i_off = 1:length(Offs)
    [h(i_off),p(i_off)]=ttest2(psth.First{i_template}(StimOffs==Offs(i_off)),psth.baseline{i_template}(StimOffs==Offs(i_off)),'Tail','right'); % the mean of X is greater than the mean of Y;
    response(1:sum(StimOffs==Offs(i_off)),i_off) = psth.First{i_template}(StimOffs==Offs(i_off));
   
    end
%     % ttest between each other groups
%     k=1;
%     for i_off = 1:(length(Offs)-1)
%         for ii_off = (i_off+1):length(Offs)
%           
%             [a(k),b(k)]=ttest2(psth.First{i_template}(StimOffs==Offs(i_off)),psth.First{i_template}(StimOffs==Offs(ii_off)),'Tail','right'); 
%             [A(k),B(k)]=ttest2(psth.First{i_template}(StimOffs==Offs(i_off)),psth.First{i_template}(StimOffs==Offs(ii_off)),'Tail','left');% the mean of X and the mean of Y is not equal;
%             [C(k),~]=ttest2(psth.First{i_template}(StimOffs==Offs(i_off)),psth.First{i_template}(StimOffs==Offs(ii_off)));
%             k=k+1;
%         end 
%     end 
    
    % for the target after 4s also should show significant
    % response
    for i_off = [5]
    [h(i_off+1),p(i_off+3)]=ttest2(psth.target{i_template}(StimOffs==Offs(i_off)),psth.target_baseline{i_template}(StimOffs==Offs(i_off)),'Tail','right');
    end
    if sum(h)==6
       
        PP.cellIX(i) = i_template;
        i=i+1;
    end
    
    [P] = anova1(response);
    close all
    
    % need to response to the first pulse, no significant reponse to all
    % first pulses and also has significant response to 4s target
    if H==1 && P >0.05 
       PP.cellanova(j)= i_template;
       j=j+1;
    end 
    if H==1 && P >0.05 && h(6)==1 
       PP.anovastrict(jj)= i_template;
       jj=jj+1;
    end 
    
%     if H==1 && sum(a)==0 && sum(A)==0
%        PP.cellttest(n) = i_template;
%        n=n+1;
%         
%     end
%     
%     if H==1 && sum(C)==0
%        PP.cellttestloose(ii) = i_template;
%        ii=ii+1;
%     end
   
end

%% get the layers of template and rpv of each cells
% for V1
% for i_template = 1:length(psth.area)
%  
%     if  spikes.template_pos(i_template,2)<=LFP.layers(1)
%         PP.layers(i_template) = 1;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(1) && spikes.template_pos(i_template,2)<=LFP.layers(2)
%         PP.layers(i_template) = 3;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(2) && spikes.template_pos(i_template,2)<=LFP.layers(3)
%         PP.layers(i_template) = 4;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(3) && spikes.template_pos(i_template,2)<=LFP.layers(4)
%         PP.layers(i_template) = 5;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(4)
%         PP.layers(i_template) = 6;
%     end
%         
%     
%    
% end 

% for HVAs, 1: above the initial big sink
% 2: in the initial sink
% 3: in between the superficial and deep sink
% 4: the deep sink or below it
% for HVAs, 1: above the SGS
% 2:IN SGS / 3:in SO / 4: below SO
 % for LGN, 1: above the sink, 2:in the sink, 3:in the source, 4: below the
 % source
% for i_template = 1:length(psth.area)
%  
%     if  spikes.template_pos(i_template,2)<=LFP.layers(1)
%         PP.layers(i_template) = 1;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(1) && spikes.template_pos(i_template,2)<=LFP.layers(2)
%         PP.layers(i_template) = 2;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(2) && spikes.template_pos(i_template,2)<=LFP.layers(3)
%         PP.layers(i_template) = 3;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(3) 
%         PP.layers(i_template) = 4;
% 
%     end
%         
%     
%    
% end


% for layers accross shanks for HVAs
for i_template = 1:length(psth.area)
    i_shank = [];
    if round(spikes.template_pos(i_template,1))==0
        i_shank =1;
    elseif  round(spikes.template_pos(i_template,1)) ==400
        i_shank =2;
    elseif round(spikes.template_pos(i_template,1)) ==800
        i_shank = 3;
    elseif round(spikes.template_pos(i_template,1)) ==1200
        i_shank = 4;
    end 
 
    if  spikes.template_pos(i_template,2)<=LFP.layers(i_shank,1)
        PP.layers(i_template) = 1;
    elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,1) && spikes.template_pos(i_template,2)<=LFP.layers(i_shank,2)
        PP.layers(i_template) = 2;
    elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,2) && spikes.template_pos(i_template,2)<=LFP.layers(i_shank,3)
        PP.layers(i_template) = 3;
    elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,3) 
        PP.layers(i_template) = 4;

    end
        
    
   
end

% for layers accross shanks for V1
% for i_template = 1:length(psth.area)
%     i_shank = [];
%     if round(spikes.template_pos(i_template,1))==0
%         i_shank =1;
%     elseif  round(spikes.template_pos(i_template,1)) ==400
%         i_shank =2;
%     elseif round(spikes.template_pos(i_template,1)) ==800
%         i_shank = 3;
%     elseif round(spikes.template_pos(i_template,1)) ==1200
%         i_shank = 4;
%     end 
%  
%     if  spikes.template_pos(i_template,2)<=LFP.layers(i_shank,1)
%         PP.layers(i_template) = 1;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,1) && spikes.template_pos(i_template,2)<=LFP.layers(i_shank,2)
%         PP.layers(i_template) = 3;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,2) && spikes.template_pos(i_template,2)<=LFP.layers(i_shank,3)
%         PP.layers(i_template) = 4;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,3) && spikes.template_pos(i_template,2)<=LFP.layers(i_shank,4)
%         PP.layers(i_template) = 5;
%     elseif spikes.template_pos(i_template,2)>LFP.layers(i_shank,4) 
%         PP.layers(i_template) = 6;
%     end
%         
%     
%    
% end



% PP.rpv= spikes.isi_vilation;

PP.baseline = baseline;
PP.binwidth = binwidth;
PP.edges = edges;
PP.Offs = Offs;

%% check the first response 
% for i_template =PP.anovastrict%1:length(psth.area)
%   figure
%  plot(edges,squeeze(PP.Traces.raw_mean(i_template,6,:)))
%  hold on
%  vline([edges(base) edges(PP.responseIX(i_template,1))])
% %  ylim([0 1.2])
%  title(['template:' num2str(i_template)])
%  
% end


%% 
for i_template =PP.cellanova%PP.cellIX%1:length(psth.area)PP.cellanova
  figure
  h=errorbar([Offs,8000],PP.response.mean(i_template,:)./binwidth,PP.response.sem(i_template,:)./binwidth,'Color',[0 0 0]);
  h.LineStyle = 'none';

  hold on
  scatter([Offs,8000], PP.response.mean(i_template,:)./binwidth,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1])

 
%     ylim([0 120])
   set(gca,'Xtick',[0 2000 4000 8000])
   set(gca,'Xticklabel',{'0';'2';'4';'baseline'})
  ylabel('FR(Hz)')
  xlabel('Inter-stimulus interval (s)')
 title(['template:' num2str(i_template) 'layers' num2str(PP.layers(i_template))])
 
end
 %% plot for check 
% for i_template =PP.cellIX%1:length(psth.area)
% figure
% 
% for i_off = 1:length(Offs)
%     supertitle(['template: ' num2str(i_template)])
% subplot(2,3,i_off)
% plot(edges_full,squeeze(PP.full.Trace_mean(i_template,i_off,:)))
% xlim([-0.1 5])
% end
% end
% %% get the real plot for an example template
% for i_template =PP.cellanova
% % i_template = 11;
% y_lim = 20;
% y_low = -10;
% figure
% supertitle(['template:' num2str(i_template) '  Layers-' num2str(PP.layers(i_template))])
% for i_off = 1:length(Offs)
% 
% if i_off<=4
%  subplot(6,1,i_off)
% shadedErrorBar_ch(edges_full,squeeze(PP.full.Trace_mean(i_template,i_off,:))./binwidth,squeeze(PP.full.Trace_sem(i_template,i_off,:))./binwidth)
% xlim([-0.2 3])
% ylim([y_low y_lim])
% else
%  subplot(6,1,5)
% shadedErrorBar_ch(edges_full,squeeze(PP.full.Trace_mean(i_template,5,:))./binwidth,squeeze(PP.full.Trace_sem(i_template,5,:))./binwidth)
% xlim([-0.2 3])
% ylim([y_low y_lim])
% subplot(6,1,6)
% shadedErrorBar_ch(edges_full,squeeze(PP.full.Trace_mean(i_template,5,:))./binwidth,squeeze(PP.full.Trace_sem(i_template,5,:))./binwidth)
% xlim([1.8 5])
% ylim([y_low y_lim])
%     
% end
% end
% end
end

%% save it 
save('ppPM.mat','PP')


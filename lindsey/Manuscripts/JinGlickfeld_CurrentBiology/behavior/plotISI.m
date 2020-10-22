%% plot different led conditions
% note that when led condition, output is Output{1} non led output{2} LED
 %output = Output{1};
%% quickly plot
ou = output.target;
Orien = output.Infor.Orien;
Off = output.Infor.Off;
colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
figure
supertitle({'  ',['i' num2str(output.Infor.ID)]});
%a=21;
a=10;
%a=0;
subplot(1,2,1)


for i_off =[1 3]%1:length(Off)
    % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
    plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
    
    hold on
    
    
end
for i_off =[1 3]%1:length(Off)
    %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
    scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
    hold on
    for i_orien = 1:length(Orien)
        % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
        line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
        hold on
    end
    
end


axis([a 100 0 1])

set(gca,'xscale','log','TickDir','out')

xlabel('Orientation change (Deg)')
ylabel('Hit Rate')
axis square
subplot(1,2,2)
scatter(Off, output.FA.FA,'k')
hold on
for i_off=1:length(Off)
    
    line([Off(i_off) Off(i_off)],output.FA.FA_confi(i_off,1:2),'Color',[0 0 0])
    
end
xlim([200 800])
ylim([0 0.2])
ylabel('FA rate')
set(gca,'XTick',Off,'TickDir','out')
axis square


%% plot the hit rate confidence intervals over orientation change on different off times 
  ou = output.target;
  Orien = output.Infor.Orien;
  Off = output.Infor.Off;
  colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
  F=figure;
  set(F,'position',[50 50 550 660])
  supertitle({'  ',['i' num2str(output.Infor.ID) '-' 'target']});
 %a=21;
  a=10;
  %a=0;
  subplot(3,2,1)

  
  for i_off =1:length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   
  end
  for i_off =1:length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 100 0 1])
    
  set(gca,'xscale','log')
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-all trials')
 
  subplot(3,2,2)   
  for i_off = 1: length(Off)%[1 3]
      
   errorbar(Orien, ou.RTonHit_mean {i_off,1}(1,:),ou.RTonHit_mean {i_off,1}(2,:),'Color',colorchoice{i_off})
   hold on
  end
  ylim([100 500])
  xlim([a 93])
  set(gca,'xscale','log','XTick', round(Orien))
  xlabel('Orientation change degree')
  ylabel('RT on Hits')
  title('RT on hits-all trials')
  
  

  % plot the short cycles separate on two offs
   subplot(3,2,3)

  
  for i_off =[1 3]%1: length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   
   hold on
   
   
  end
  for i_off =[1 3]%1: length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  h=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  h.MarkerFaceColor = colorchoice{i_off};
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-2-4 cycles')
  
% plot the long cycles separate on two offs
   subplot(3,2,4)

  
  for i_off =[1 3]%1: length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   
   hold on
   
   
  end
  for i_off =[1 3]%1: length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-6-9 cycles')

% plot the short cycle vs long cycle on 250ms    
  subplot(3,2,5)
  
   for i_off =1%1: length(Off)
    %plot(Orien, L_hit{i_off,1},'Color',colorchoice{i_off})
%      plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
     hold on
    plot(Orien,ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   hold on
    plot(Orien,ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off},'LineStyle',':');
   
   
   
  end
  for i_off =1%1: length(Off)
     
%     scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   s=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   s.MarkerFaceColor = colorchoice{i_off}+0.5;
   hold on
   scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
%    line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   
   hold on
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off},'LineStyle',':');
   end
   
  end
  
  
  axis([a 93 0 1])
%   legend(legend_infor)  
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('250:Filled:2-4C;Open:6-9C')

% plot the short cycle vs long cycle on 750ms    
  subplot(3,2,6)
  
   for i_off =3%1: length(Off)
    %plot(Orien, L_hit{i_off,1},'Color',colorchoice{i_off})
%      plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
     hold on
    plot(Orien,ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   hold on
    plot(Orien,ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off},'LineStyle',':');
   
   
   
  end
  for i_off =3%1: length(Off)
     
%     scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   s=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   s.MarkerFaceColor = colorchoice{i_off};
   hold on
   scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
%    line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   
   hold on
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off},'LineStyle',':');
   end
   
  end
  
  
  axis([a 93 0 1])

%   legend(legend_infor)  
  set(gca,'xscale','log','XTick', round(Orien))
 
 xlabel('Orientation change degree')

  
  ylabel('Hit Rate')
  title('750:Filled:2-4C;Open:6-9C')   
  
  %% plot the RT on hit 
     
  
%     figure
%    
%     
%    for i_off = 1: length(Off)%[1 3]
%       
%    errorbar(Orien, ou.RTonHit_mean {i_off,1}(1,:),ou.RTonHit_mean {i_off,1}(2,:),'Color',colorchoice{i_off})
%    hold on
%   end
%   ylim([200 550])
%   xlim([10 100])
%   set(gca,'xscale','log','XTick', [10:10:100])
%   set(gca,'XTickLabel',{'10','','','', '50','','','','','100'})
%   set(gca,'TickDir','out')
%   xlabel('Orientation change degree')
%   ylabel('RT on Hits')
%%   
ou = output.pre;
 
  colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
  F=figure;
  set(F,'position',[50 50 550 660])
   h = supertitle({'  ',['i' num2str(output.Infor.ID) '-' 'pre']});
 
  subplot(3,2,1)

  
  for i_off =1: length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   
  end
  for i_off =1: length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-all trials')
 
  subplot(3,2,2)   
  for i_off = [1 3]
      
   errorbar(Orien, ou.RTonHit_mean {i_off,1}(1,:),ou.RTonHit_mean {i_off,1}(2,:),'Color',colorchoice{i_off})
   hold on
  end
  ylim([100 550])
  xlim([a 93])
  set(gca,'xscale','log','XTick', round(Orien))
  xlabel('Orientation change degree')
  ylabel('RT on Hits')
  title('RT on hits-all trials')
  
  

  % plot the short cycles separate on two offs
   subplot(3,2,3)

  
  for i_off =[1 3]%1: length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   
   hold on
   
   
  end
  for i_off =[1 3]%1: length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  h=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  h.MarkerFaceColor = colorchoice{i_off};
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-2-4 cycles')
  
% plot the long cycles separate on two offs
   subplot(3,2,4)

  
  for i_off =[1 3]%1: length(Off)
   % plot(Orien, S_hit{i_off,1},'Color',colorchoice{i_off})
   plot(Orien, ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   
   hold on
   
   
  end
  for i_off =[1 3]%1: length(Off)
      %scatter(Orien, S_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-6-9 cycles')

% plot the short cycle vs long cycle on 250ms    
  subplot(3,2,5)
  
   for i_off =1%1: length(Off)
    %plot(Orien, L_hit{i_off,1},'Color',colorchoice{i_off})
%      plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
     hold on
    plot(Orien,ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   hold on
    plot(Orien,ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off},'LineStyle',':');
   
   
   
  end
  for i_off =1%1: length(Off)
     
%     scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   s=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   s.MarkerFaceColor = colorchoice{i_off}+0.5;
   hold on
   scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
%    line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   
   hold on
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off},'LineStyle',':');
   end
   
  end
  
  
  axis([a 93 0 1])
%   legend(legend_infor)  
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('250:Filled:2-4C;Open:6-9C')

% plot the short cycle vs long cycle on 750ms    
  subplot(3,2,6)
  
   for i_off =3%1: length(Off)
    %plot(Orien, L_hit{i_off,1},'Color',colorchoice{i_off})
%      plot(Orien, ou.HT_rate{i_off,1}, 'Color',colorchoice{i_off});
     hold on
    plot(Orien,ou.S_HT_rate{i_off,1}, 'Color',colorchoice{i_off});
   hold on
    plot(Orien,ou.L_HT_rate{i_off,1}, 'Color',colorchoice{i_off},'LineStyle',':');
   
   
   
  end
  for i_off =3%1: length(Off)
     
%     scatter(Orien, ou.HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   hold on
   s=scatter(Orien, ou.S_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   s.MarkerFaceColor = colorchoice{i_off};
   hold on
   scatter(Orien, ou.L_HT_rate{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   for i_orien = 1:length(Orien)
   % line([Orien(i_orien) Orien(i_orien)],L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
%    line([Orien(i_orien) Orien(i_orien)],ou.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   line([Orien(i_orien) Orien(i_orien)],ou.S_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off});
   
   hold on
   line([Orien(i_orien) Orien(i_orien)],ou.L_confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off},'LineStyle',':');
   end
   
  end
  
  
  axis([a 93 0 1])
  V=axis;
%   legend(legend_infor)  
  set(gca,'xscale','log','XTick', round(Orien))
 
 xlabel('Orientation change degree')

  
  ylabel('Hit Rate')
  title('750:Filled:2-4C;Open:6-9C')     


%% plot FA rate
  
   F=figure; 
   set(F,'position',[50 50 800 600])
   subplot(2,2,1)
   scatter(Off, output.FA.FA,'k')
   hold on
   for i_off=1:length(Off)
   
   line([Off(i_off) Off(i_off)],output.FA.FA_confi(i_off,1:2),'Color',[0 0 0])
   
   end
   xlim([200 800])
   ylim([0 0.2])
   ylabel('FA rate')   
  set(gca,'XTick',Off)
  title('FA-all trials')
  
   subplot(2,2,2)
  
  for i_off = 1: length(Off)
      
   errorbar(Off(i_off),mean(output.FA.FA_RT{i_off,1}),std(double(output.FA.FA_RT{i_off,1}))./sqrt(length(output.FA.FA_RT{i_off,1})),'Color',[0 0 0])
   hold on
   scatter(Off(i_off),mean(output.FA.FA_RT{i_off,1}),'k')
  
  end
  ylim([200 550])
  xlim([200 800])
  set(gca,'XTick',Off)
  set(gca,'TickDir','out')
  
  ylabel('RT on FAs')
  title('RT on FAs-all trials')
  
  subplot(2,2,3)
  s=scatter(Off, output.FA.S_FA,'k');
  s.MarkerFaceColor = [0.8 0.8 0.8];
   hold on
  scatter(Off, output.FA.L_FA,'k')
   for i_off=1:length(Off)
   
   line([Off(i_off) Off(i_off)],output.FA.S_FA_confi(i_off,1:2),'Color',[0 0 0])
   line([Off(i_off) Off(i_off)],output.FA.L_FA_confi(i_off,1:2),'Color',[0 0 0])
   end
   xlim([200 800])
   ylim([0 0.2])
   ylabel('FA rate')   
  set(gca,'XTick',Off)
  title('FA-close:2-4C-open:6-8C')
  
  subplot(2,2,4)
  
  for i_off = 1: length(Off)
      
   errorbar(Off(i_off),mean(output.FA.S_FA_RT{i_off,1}),std(double(output.FA.S_FA_RT{i_off,1}))./sqrt(length(output.FA.S_FA_RT{i_off,1})),'Color',[0 0 0])
   hold on
   s=scatter(Off(i_off),mean(output.FA.S_FA_RT{i_off,1}),'k');
   s.MarkerFaceColor = [0.8 0.8 0.8];
    
   errorbar(Off(i_off),mean(output.FA.L_FA_RT{i_off,1}),std(double(output.FA.L_FA_RT{i_off,1}))./sqrt(length(output.FA.L_FA_RT{i_off,1})),'Color',[0 0 0])
  
   scatter(Off(i_off),mean(output.FA.L_FA_RT{i_off,1}),'k');
  
  end
  ylim([0 550])
  xlim([200 800])
  set(gca,'XTick',Off)
  
  ylabel('RT on FAs')
  title('RTonFA-close:2-4C-open:6-8C')
  
  supertitle(['i' num2str(output.Infor.ID) '-' 'FA'])
  
  %% plot release in 100ms visual presentation 

   F=figure; 
   set(F,'position',[50 50 800 600])
   subplot(2,2,1)
   scatter(Off, output.Re.Re,'k')
   hold on
   for i_off=1:length(Off)
   
   line([Off(i_off) Off(i_off)],output.Re.Re_confi(i_off,1:2),'Color',[0 0 0])
   
   end
   xlim([200 800])
   ylim([0 0.2])
   ylabel('Rate within 100ms')   
  set(gca,'XTick',Off)
  title('100ms-all trials')
  
   subplot(2,2,2)
  
  for i_off = 1: length(Off)
      
   errorbar(Off(i_off),mean(output.Re.Re_RT{i_off,1}),std(double(output.Re.Re_RT{i_off,1}))./sqrt(length(output.Re.Re_RT{i_off,1})),'Color',[0 0 0])
   hold on
   scatter(Off(i_off),mean(output.Re.Re_RT{i_off,1}),'k')
  
  end
  ylim([0 550])
  xlim([200 800])
  set(gca,'XTick',Off)
  
  ylabel('RT on 100ms')
  title('RT on 100ms-all trials')
  
  subplot(2,2,3)
  s=scatter(Off, output.Re.S_Re,'k');
  s.MarkerFaceColor = [0.8 0.8 0.8];
   hold on
  scatter(Off, output.Re.L_Re,'k')
   for i_off=1:length(Off)
   
   line([Off(i_off) Off(i_off)],output.Re.S_Re_confi(i_off,1:2),'Color',[0 0 0])
   line([Off(i_off) Off(i_off)],output.Re.L_Re_confi(i_off,1:2),'Color',[0 0 0])
   end
   xlim([200 800])
   ylim([0 0.2])
   ylabel('Rate 100 ms')   
  set(gca,'XTick',Off)
  title('response 100ms-close:2-4C-open:6-8C')
  
  subplot(2,2,4)
  
  for i_off = 1: length(Off)
      
   errorbar(Off(i_off),mean(output.Re.S_Re_RT{i_off,1}),std(double(output.Re.S_Re_RT{i_off,1}))./sqrt(length(output.Re.S_Re_RT{i_off,1})),'Color',[0 0 0])
   hold on
   s=scatter(Off(i_off),mean(output.Re.S_Re_RT{i_off,1}),'k');
   s.MarkerFaceColor = [0.8 0.8 0.8];
    
   errorbar(Off(i_off),mean(output.Re.L_Re_RT{i_off,1}),std(double(output.Re.L_Re_RT{i_off,1}))./sqrt(length(output.Re.L_Re_RT{i_off,1})),'Color',[0 0 0])
  
   scatter(Off(i_off),mean(output.Re.L_Re_RT{i_off,1}),'k');
  
  end
  ylim([0 550])
  xlim([200 800])
  set(gca,'XTick',Off)
  
  ylabel('RT on 100 ms')
  title('RTon100ms -close:2-4C-open:6-8C')
  
  supertitle(['i' num2str(output.Infor.ID) '-' 'FA'])
 %% plot dprime and c
  ou = output.sdt;
 
  colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
  F=figure;
  set(F,'position',[50 200 800 300])
  supertitle({'  ',['i' num2str(output.Infor.ID) '-' 'target']});
  subplot(1,2,1)
   
  for i_off =1: length(Off)
  
   plot(Orien, ou.dprime{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   scatter(Orien, ou.dprime{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   
  end 
  
 
    set(gca,'xscale','log','XTick',[10:10:100])
    set(gca,'XTickLabel',{'10','','','', '50','','','','','100'})
    set(gca,'TickDir','out')
   xlim([10 100])
   ylim([0 5])
  xlabel('Orientation change degree')
  ylabel('d prime')
  title('d prime-all trials')
  
    subplot(1,2,2)
   
  for i_off =1: length(Off)
  
   plot(Orien, ou.criterion{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   scatter(Orien, ou.criterion{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   
  end 
  
 
    set(gca,'xscale','log','XTick',[10:10:100])
    set(gca,'XTickLabel',{'10','','','', '50','','','','','100'})
    set(gca,'TickDir','out')
   xlim([10 100])
   ylim([-1 2])
   hline(0,'k:')
  xlabel('Orientation change degree')
  ylabel('d prime')
  title('d prime-all trials')
  
  
 
  %% plot d prime only
  ou = output.sdt;
 
  colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
  F=figure;
  set(F,'position',[50 200 800 300])
  supertitle({'  ',['i' num2str(output.Infor.ID) '-' 'target']});
 
  b=5;

  subplot(1,3,1)

  
  for i_off =1: length(Off)
  
   plot(Orien, ou.dprime{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   scatter(Orien, ou.dprime{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   
  end 
  
  axis([a 93 0 b])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('d prime')
  title('d prime-all trials')
 

  % plot the short cycles separate on two offs
   subplot(1,3,2)

  for i_off =[1 2]%1: length(Off)
  
   plot(Orien, ou.S_dprime{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
    scatter(Orien, ou.S_dprime{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   
  end 
 
  
  axis([a 93 0 b])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('d prime')
  title('d prime-2-4 cycles')
  
% plot the long cycles separate on two offs
   subplot(1,3,3)

   for i_off =[1 2]%1: length(Off)
   scatter(Orien, ou.L_dprime{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
    hold on
   plot(Orien, ou.L_dprime{i_off,1}, 'Color',colorchoice{i_off});
  
  
   
   
  end 
 
  
  axis([a 93 0 b])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('d prime')
  title('d prime-6-8 cycles')
 %% plot criterion change 




  figure
  for i_off =1: length(Off)
  
   plot(Orien, ou.criterion{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
   
   scatter(Orien, ou.criterion{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
   
  end 
  hold on 
  hline(0,'k:')
  
  axis([a 93 -1 2])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('c criterion')
  title('c criterion')
 
  %% plot RT distribution of ISIs
  skew = [];
  Off = output.Infor.Off;
  figure
subplot(2,2,1)
  for i_off = 1: length(Off)
       h = cdfplot( output.RT.B_RT{i_off});
       h.Color = colorchoice{i_off};
      
       hold on
       skew (i_off,1) = skewness(output.RT.B_RT{i_off}); 
%        for i_orien = 1:5
%            temp =[temp output.target.RTonHit{i_off,1}{i_orien,1}];
%        end
  end
  
  
  

 
  set(gca,'TickDir','out')
  set(gca,'XTick',0:100:600)
  legend('250ms-FA','500ms-FA','750ms-FA','Location','northwest')
  legend boxoff
  grid off
  title(['i' num2str(output.Infor.ID) '-' 'RT distribution at FAs']);
  ylabel('Fraction of RTs')
  xlabel('Reaction time (ms)')
  xlim([0 600])
  set(gca,'XTick',0:200:600,'TickDir','out')
  vline([200 550],'k:')
  axis square
  
  subplot(2,2,2)
  h = cdfplot(output.RT.B_RT{3});
  h.Color = colorchoice{3};
  hold on
  h=cdfplot(output.RT.T_orien{1}(output.RT.T_orien{1}<=550 ));
  skew (4,1) = skewness(output.RT.T_orien{1}(output.RT.T_orien{1}<=550));
  h.Color = [0 0 0]
  set(gca,'TickDir','out')
  set(gca,'XTick',0:100:600)
  legend('7500ms-FA','22.5 Deg-Hit','Location','northwest')
  legend boxoff
  grid off
  title(['i' num2str(output.Infor.ID) '-' 'RT distribution at FAs']);
  ylabel('Fraction of RTs')
  xlabel('Reaction time (ms)')
  xlim([0 600])
  vline([200 550],'k:')
  set(gca,'XTick',0:200:600,'TickDir','out')
  axis square
  
  subplot(2,2,3)
  scatter([250 500 750 1000], skew,'MarkerEdgeColor',[0 0 0],'SizeData',80)
%   ylim([-0.2 1])
  xlim([200 1050])
  set(gca,'TickDir','out')
  set(gca, 'XTick',[250 500 750 1000])
  set(gca,'XTickLabel',{'250','500','750','22.5 deg'})
  ylabel('skewness')
 % h=cdfplot([output.target.RTonHit{1,1}{1,1} output.target.RTonHit{2,1}{1,1} output.target.RTonHit{3,1}{1,1} ]);

  
 % plot the RT that exclude the first 200 ms 
 %% get the slope
%  figure
%  [f,x] = ecdf( output.RT.B_RT{3});
% subplot(2,1,1)
%  plot(x,f)
%  
%  subplot(2,1,2)
%  plot(x,movingslope(f,200,1,0.01));
  
% [ h,p] = kstest2(output.RT.T_orien{1}(output.RT.T_orien{1}<=550),output.RT.B_RT{3})
%  hold on

  %% plot RT distribution for target collapsed for all orientation
 
   figure
  subplot(1,3,1)
  for i_off = 1: length(Off)
       cdfplot(output.RT.T_RT{i_off})
       hold on
  end
  legend('250ms','500ms','750ms','Location','northwest')
  title(['i' num2str(output.Infor.ID) '-' 'RT distribution at all Hits']);
  ylabel('Fraction of RTs')
  xlabel('Reaction time at FAs')
   xlim([0 1000])

  subplot(1,3,2)
  for i_off = 1: length(Off)
       cdfplot(output.RT.B_RT{i_off,1})
       hold on
  end
  
  title(['i' num2str(output.Infor.ID) '-' 'RT distribution at current BS']);
  ylabel('Fraction of RTs')
  xlabel('Reaction time at current BS')
   xlim([0 1000])
  
  subplot(1,3,3)
   
   for i_off = 1: length(Off)
       cdfplot(output.RT.pre_B_RT{i_off,1})
       hold on
  end
  
  title(['i' num2str(output.Infor.ID) '-' 'RT distribution at pre BS']);
  ylabel('Fraction of RTs')
  xlabel('Reaction time at pre BS')
   xlim([0 1800])
  
  %% plot Target RT and FA RT overlay for different offs
   figure
  
  for i_off = 1: length(Off)
       subplot(1,3,i_off)
       h=cdfplot(output.RT.T_RT{i_off});
       set(h, 'Color','k')
       hold on
       h=cdfplot(output.RT.B_RT{i_off});
       set(h, 'Color','r')
       ylabel('Fraction of RTs')
       xlabel('RT at Hits(blk); at FAs(red)')
       title([num2str(Off(i_off)) 'ms'])
      
  end
  
  
 %% plot hit rate seperate by n offs
 
  F=figure;
  set(F,'position',[50 200 800 800])
  supertitle(['  ','i' num2str(output.Infor.ID) '-' 'Hit rate' ':  ' 'n offs: blk-250; gray-500; yellow-750']);
  subplot(2,2,1)
  for i_off =1: length(Off)
   
   plot(Orien, output.target.c_hit{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
  scatter(Orien, output.target.c_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  
   for i_orien = 1:length(Orien)
   
   line([Orien(i_orien) Orien(i_orien)],output.target.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-all n-1 trials')
  
 for i_off = 1:length(Off)
     subplot(2,2,i_off+1)
     for ii_off = 1:length(Off)
         plot(Orien, squeeze(output.con.c_hit(i_off,ii_off,:)) , 'Color',colorchoice{ii_off});
         hold on
         scatter(Orien, squeeze(output.con.HT_rate(i_off,ii_off,:)) ,'MarkerEdgeColor',colorchoice{ii_off});
         hold on
         for i_orien = 1:length(Orien)
             
             line([Orien(i_orien) Orien(i_orien)],squeeze(output.con.confi(i_off,ii_off,i_orien,1:2)), 'Color',colorchoice{ii_off})
             hold on
             
         end
         
     end
     title([num2str(Off(i_off)) ' ' 'n-1' ' ' 'trials'])
     xlabel('Orientation change degrees')
     ylabel('Hit Rate')
       axis([a 93 0 1])
       set(gca,'xscale','log','XTick', round(Orien)) 
 end

 %% plot the hit rate seperate by n-1 offs 
 F=figure;
  set(F,'position',[50 200 800 800])
  supertitle(['  ','i' num2str(output.Infor.ID) '-' 'Hit rate' ':  ' 'n-1 offs: blk-250; gray-500; yellow-750']);
  subplot(2,2,1)
  for i_off =1: length(Off)
   
   plot(Orien, output.pre.c_hit{i_off,1}, 'Color',colorchoice{i_off});
  
   hold on
  scatter(Orien, output.pre.c_hit{i_off,1},'MarkerEdgeColor',colorchoice{i_off});
  
   for i_orien = 1:length(Orien)
   
   line([Orien(i_orien) Orien(i_orien)],output.pre.confi{i_off,1}(i_orien,1:2), 'Color',colorchoice{i_off})
   hold on
   end
   
  end
  
  
  axis([a 93 0 1])
    
  set(gca,'xscale','log','XTick', round(Orien))
 
  xlabel('Orientation change degree')
  ylabel('Hit Rate')
  title('Hit Rate-all n trials')

 for ii_off = 1:length(Off)
     subplot(2,2,ii_off+1)
     for i_off = 1:length(Off)
         plot(Orien, squeeze(output.con.c_hit(i_off,ii_off,:)) , 'Color',colorchoice{i_off});
         hold on
         scatter(Orien, squeeze(output.con.HT_rate(i_off,ii_off,:)) ,'MarkerEdgeColor',colorchoice{i_off});
         hold on
         for i_orien = 1:length(Orien)
             
             line([Orien(i_orien) Orien(i_orien)],squeeze(output.con.confi(i_off,ii_off,i_orien,1:2)), 'Color',colorchoice{i_off})
             hold on
             
         end
         
     end
   title([num2str(Off(ii_off)) ' ' 'n' ' ' 'trials'])
     xlabel('Orientation change degrees')
     ylabel('Hit Rate')
       axis([a 93 0 1])
       set(gca,'xscale','log','XTick', round(Orien)) 
 end
 
 
 %% plot FA rate seperate by n offs
 F=figure;
 set(F,'position',[50 200 800 800])
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'FA rate' ':  ' 'n offs: in X axis']);
 subplot(2,2,1)
 scatter(Off, output.con.pre.FA,'k')
 hold on
 for i_off=1:length(Off)
     
     line([Off(i_off) Off(i_off)],output.con.pre.FA_confi(i_off,1:2),'Color',[0 0 0])
     
 end
 
 ylim([0 0.2])
 xlim([200 800])
 set(gca,'XTick',Off)
 ylabel('FA Rate')
 xlabel('n Offs (ms)')
 title('FA rate-All n-1 trials')
 
 for i_off = 1:length(Off)
     subplot(2,2,i_off+1)
     scatter(Off, output.con.FA(i_off,:),'k')
     for ii_off=1:length(Off)
         
         hold on
         line([Off(ii_off) Off(ii_off)],squeeze(output.con.FA_confi(i_off,ii_off,1:2)),'Color',[0 0 0])
         
     end
     title([num2str(Off(i_off)) ' ' 'n-1' ' ' 'trials'])
     ylabel('FA Rate')
     xlabel('n Offs (ms)')
     ylim([0 0.2])
     xlim([200 800])
 end
   %% plot FA rate seperate by n-1 offs
 F=figure;
 set(F,'position',[50 200 800 800])
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'FA rate' ':  ' 'n-1 offs: in X axis']);
 subplot(2,2,1)
 scatter(Off, output.con.pre2.FA,'k')
 hold on
 for i_off=1:length(Off)
     
     line([Off(i_off) Off(i_off)],output.con.pre2.FA_confi(i_off,1:2),'Color',[0 0 0])
     
 end
 
 ylim([0 0.2])
 xlim([200 800])
 set(gca,'XTick',Off)
 ylabel('FA Rate')
 xlabel('n-1 Offs (ms)')
 title('FA rate-All n trials')
 
 
 for ii_off = 1:length(Off)
      subplot(2,2,ii_off+1)
     scatter(Off, output.con.FA(:,ii_off)','k')
     for i_off=1:length(Off)
         
         hold on
         line([Off(i_off) Off(i_off)],squeeze(output.con.FA_confi(i_off,ii_off,1:2)),'Color',[0 0 0])
         
     end
     title([num2str(Off(ii_off)) ' ' 'n' ' ' 'trials'])
     ylabel('FA Rate')
     xlabel('n-1 Offs (ms)')
     ylim([0 0.2])
     xlim([200 800])
      set(gca,'XTick',Off)
 end
 %% plot the RT difference by n offs
 F=figure;
 set(F,'position',[50 200 800 800])
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'RT at FAs' ':  ' 'n offs: blk-250; gray-500; yellow-750']);
 subplot(2,2,1)
 for i_off = 1: length(Off)
 h=cdfplot(output.con.pre.FA_RT{i_off,1});
 set(h, 'Color',colorchoice{i_off})
 hold on
 
 end
 xlim([0 500])
 xlabel('RT at FAs (ms)')
 ylabel('Fraction of RT')
 title('RT at FAs-All n-1 trials')
 
 for i_off = 1:length(Off)
     subplot(2,2,i_off+1)
    
     for ii_off=1:length(Off)
         
         h=cdfplot(output.con.FA_RT{i_off,ii_off});
        set(h, 'Color',colorchoice{ii_off})
        
         hold on
     end
     title([num2str(Off(i_off)) ' ' 'n-1' ' ' 'trials'])
     ylabel('Fraction of RT')
     xlabel('RT at FAs (ms)')
     xlim([0 500])
 end

 %% plot the RT difference by n-1 offs
 F=figure;
 set(F,'position',[50 200 800 800])
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'RT at FAs' ':  ' 'n-1 offs: blk-250; gray-500; yellow-750']);
 subplot(2,2,1)
 for i_off = 1: length(Off)
 h=cdfplot(output.con.pre2.FA_RT{i_off,1});
 set(h, 'Color',colorchoice{i_off})
 hold on
 
 end
 xlim([0 500])
 xlabel('RT at FAs (ms)')
 ylabel('Fraction of RT')
 title('RT at FAs-All n trials')
 

 for ii_off = 1:length(Off)
     subplot(2,2,ii_off+1)
    
     for i_off=1:length(Off)
         
         h=cdfplot(output.con.FA_RT{i_off,ii_off});
        set(h, 'Color',colorchoice{i_off})
        
         hold on
     end
     title([num2str(Off(ii_off)) ' ' 'n' ' ' 'trials'])
     ylabel('Fraction of RT')
     xlabel('RT at FAs (ms)')
     xlim([0 500])
 end
 %% try to plot reaction time on hits seperate by n offs
 F=figure;
 set(F,'position',[50 200 800 800])
 i_orien=5;
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'RT at Hits' '-' num2str(round(Orien(i_orien))) 'deg' ':  ' 'n offs: blk-250; gray-500; yellow-750']);

 subplot(2,2,1)
 
 
 for i_off = 1: length(Off)
 h=cdfplot(output.target.RTonHit{i_off,1}{i_orien,1});
 set(h, 'Color',colorchoice{i_off})
 hold on
 
 end
 xlim([0 700])
 xlabel(['RT at Hits (ms)' 'target-' num2str(round(Orien(i_orien))) 'deg'])
 ylabel('Fraction of RT')
 title('RT at Hits-All n-1 trials')
 
 for i_off = 1:length(Off)
     subplot(2,2,i_off+1)
    
     for ii_off=1:length(Off)
        
         h=cdfplot(output.con.RTonHit{i_off,1}{ii_off,1}{i_orien,1});
        set(h, 'Color',colorchoice{ii_off})
        
         hold on
     end
     title([num2str(Off(i_off)) ' ' 'n-1' ' ' 'trials'])
     ylabel('Fraction of RT')
      xlabel(['RT at Hits (ms)' 'target-' num2str(round(Orien(i_orien))) 'deg'])
     xlim([0 700])
 end
 
 %% try to plot reaction time on hits seperate by n-1 offs
 F=figure;
 set(F,'position',[50 200 800 800])
 
 supertitle(['  ','i' num2str(output.Infor.ID) '-' 'RT at Hits' '-' num2str(round(Orien(i_orien))) 'deg' ':  ' 'n-1 offs: blk-250; gray-500; yellow-750']);

 subplot(2,2,1)
 
 
 for i_off = 1: length(Off)
 h=cdfplot(output.pre.RTonHit{i_off,1}{i_orien,1});
 set(h, 'Color',colorchoice{i_off})
 hold on
 
 end
 xlim([0 700])
 xlabel(['RT at Hits (ms)' 'target-' num2str(round(Orien(i_orien))) 'deg'])
 ylabel('Fraction of RT')
 title('RT at Hits-All n trials')
 
 for ii_off = 1:length(Off)
     subplot(2,2,ii_off+1)
    
     for i_off=1:length(Off)
        
         h=cdfplot(output.con.RTonHit{i_off,1}{ii_off,1}{i_orien,1});
        set(h, 'Color',colorchoice{i_off})
        
         hold on
     end
     title([num2str(Off(ii_off)) ' ' 'n' ' ' 'trials'])
     ylabel('Fraction of RT')
      xlabel(['RT at Hits (ms)' 'target-' num2str(round(Orien(i_orien))) 'deg'])
     xlim([0 700])
 end
 %% plot cycle distribution
 colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
 
 if isfield (output.FA,'cycle')
     figure
     subplot(1,2,1)
     for i = [3 2 1]
     h= histogram(output.FA.cycle{i,1});
     h.FaceColor = colorchoice{i};
     h.EdgeColor = [1 1 1];
     h.FaceAlpha = 1;
     hold on
     end
     xlim([1 11])
     xlabel('Cycle Number')
     ylabel('FA number')
     subplot(1,2,2)
     for i = [3 2 1]
         h= cdfplot(output.FA.cycle{i,1});
         h.Color =colorchoice{i}; 
         hold on
         
     end
     grid off
     xlim([1 11])
     xlabel('Cycle Number')
     ylabel('Fraction of FAs')
 end

 
  
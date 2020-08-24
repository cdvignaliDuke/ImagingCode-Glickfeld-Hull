function [output] = ISI_HR(StimOff,Off,Orien,Success,Miss,tCycleNum,tOrientation,RT_HM)

 for i_off = 1:length(Off)
      for i_orien = 1: length(Orien)
          % find all the corrects on orientation changes, ignore the cycle
          % number
             
          output.HT_num{i_off,1}(i_orien,1) = sum(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien)));
          output.Miss_num{i_off,1}(i_orien,1) = sum (Miss &(StimOff==Off(i_off))&(tOrientation==Orien(i_orien)));
          output.HT_rate{i_off,1}(i_orien,1) = output.HT_num{i_off,1}(i_orien,1) ./ (output.HT_num{i_off,1}(i_orien,1)+output.Miss_num{i_off,1}(i_orien,1));

          % use the newly corrected reaction time
          [a,b]= binofit(output.HT_num{i_off,1}(i_orien,1),(output.HT_num{i_off,1}(i_orien,1)+output.Miss_num{i_off,1}(i_orien,1)));
          output.confi{i_off,1}(i_orien,1:2) = b;
          output.c_hit{i_off,1}(i_orien,1) = a;
          % try to seperate by cycle number, short:2-4; long 6-10
          output.S_HT_num{i_off,1}(i_orien,1) = sum(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=2)&(tCycleNum<=4));
          output.S_Miss_num{i_off,1}(i_orien,1) = sum (Miss &(StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=2)&(tCycleNum<=4));
          output.S_HT_rate{i_off,1}(i_orien,1) = output.S_HT_num{i_off,1}(i_orien,1) ./ (output.S_HT_num{i_off,1}(i_orien,1)+output.S_Miss_num{i_off,1}(i_orien,1));
          
          output.L_HT_num{i_off,1}(i_orien,1) = sum(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=6)&(tCycleNum<=9));
          output.L_Miss_num{i_off,1}(i_orien,1) =  sum (Miss &(StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=6)&(tCycleNum<=9));
          output.L_HT_rate{i_off,1}(i_orien,1) = output.L_HT_num{i_off,1}(i_orien,1) ./ (output.L_HT_num{i_off,1}(i_orien,1)+output.L_Miss_num{i_off,1}(i_orien,1));
          
          
          [e,f] = binofit(output.S_HT_num{i_off,1}(i_orien,1),(output.S_HT_num{i_off,1}(i_orien,1)+ output.S_Miss_num{i_off,1}(i_orien,1)));
          
          output.S_confi{i_off,1}(i_orien,1:2) = f;
          output.S_hit{i_off,1}(i_orien,1) = e;
          [g,h] = binofit(output.L_HT_num{i_off,1}(i_orien,1),(output.L_HT_num{i_off,1}(i_orien,1)+ output.L_Miss_num{i_off,1}(i_orien,1)));
          
          output.L_confi{i_off,1}(i_orien,1:2) = h;
          output.L_hit{i_off,1}(i_orien,1) = g;
          % RT_on hit for all trials 
          output.RTonHit {i_off,1}{i_orien,1} = RT_HM(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien)));
          output.RTonHit_mean {i_off,1}(1,i_orien) = mean(output.RTonHit {i_off,1}{i_orien,1});
          output.RTonHit_mean {i_off,1}(2,i_orien) = std(output.RTonHit {i_off,1}{i_orien,1})./sqrt(output.HT_num{i_off,1}(i_orien,1));
          % RT_on hit for short trials
          output.S_RTonHit {i_off,1}{i_orien,1} = RT_HM(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=2)&(tCycleNum<=4));
          output.S_RTonHit_mean {i_off,1}(1,i_orien) = mean(output.S_RTonHit {i_off,1}{i_orien,1});
          output.S_RTonHit_mean {i_off,1}(2,i_orien) = std(output.S_RTonHit {i_off,1}{i_orien,1})./sqrt(output.S_HT_num{i_off,1}(i_orien,1));
          % RT_on hit for long trials
          output.L_RTonHit {i_off,1}{i_orien,1} = RT_HM(Success & (StimOff==Off(i_off))&(tOrientation==Orien(i_orien))&(tCycleNum>=6)&(tCycleNum<=9));
          output.L_RTonHit_mean {i_off,1}(1,i_orien) = mean(output.L_RTonHit {i_off,1}{i_orien,1});
          output.L_RTonHit_mean {i_off,1}(2,i_orien) = std(output.L_RTonHit {i_off,1}{i_orien,1})./sqrt(output.L_HT_num{i_off,1}(i_orien,1));
      end
      
      
  end

%% collapsed all the off times 

 for i_orien = 1: length(Orien)
          % find all the corrects on orientation changes, ignore the cycle
          % number
             
          output.all.HT_num(i_orien,1) = sum(Success &(tOrientation==Orien(i_orien)));
          output.all.Miss_num(i_orien,1) = sum (Miss &(tOrientation==Orien(i_orien)));
          output.all.HT_rate(i_orien,1) = output.all.HT_num(i_orien,1)./ (output.all.HT_num(i_orien,1)+output.all.Miss_num(i_orien,1));

          % use the newly corrected reaction time
          [a,b]= binofit(output.all.HT_num(i_orien,1),(output.all.HT_num(i_orien,1)+output.all.Miss_num(i_orien,1)));
          output.all.confi(i_orien,1:2) = b;
          output.all.c_hit(i_orien,1) = a;
          
          output.all.RTonHit{i_orien,1} = RT_HM(Success &(tOrientation==Orien(i_orien)));
          output.all.RTonHit_mean(1,i_orien) = mean(output.all.RTonHit{i_orien,1});
          output.all.RTonHit_mean(2,i_orien) = std(output.all.RTonHit{i_orien,1})./sqrt(output.all.HT_num(i_orien,1));
          
 end


end
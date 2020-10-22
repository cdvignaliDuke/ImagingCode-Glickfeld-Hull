function [out] = Fixed_sdt(Orien,output)
% calculate d prime and criterion for different orientations
%Adjust only the extreme values by replacing rates of 0 with 0.5/n  and rates
%of 1 with (n-0.5)/n where nn is the number of signal or noise trials (Macmillan & Kaplan, 1985)
% calculate d' and c for collapsed all conditions
for i_orien = 1:length(Orien)
    if output.target.all.c_hit(i_orien,1)<1
    [out.all.dprime(i_orien), out.all.criterion(i_orien)] = dprime_simple(output.target.all.c_hit(i_orien,1),output.FA.all.FA);
    else
    signaltrialN = output.target.all.HT_num(i_orien,1) + output.target.all.Miss_num(i_orien,1);
    [out.all.dprime(i_orien), out.all.criterion(i_orien)] = dprime_simple((signaltrialN-0.5)./signaltrialN,output.FA.all.FA);
    end 
 
end 
  
    
end
function [out] = Contrast_sdt(Contrast,output)
% calculate d prime and criterion for different orientations
%Adjust only the extreme values by replacing rates of 0 with 0.5/n  and rates
%of 1 with (n-0.5)/n where nn is the number of signal or noise trials (Macmillan & Kaplan, 1985)
% calculate d' and c for collapsed all conditions
for i_contrast = 1:length(Contrast)
    if output.target.c_hit(i_contrast,1)<1
    [out.dprime(i_contrast), out.criterion(i_contrast)] = dprime_simple(output.target.c_hit(i_contrast,1),output.FA.all.FA);
    else
    signaltrialN = output.target.HT_num(i_contrast,1) + output.target.Miss_num(i_contrast,1);
    [out.dprime(i_contrast), out.criterion(i_contrast)] = dprime_simple((signaltrialN-0.5)./signaltrialN,output.FA.all.FA);
    end 
 
end 
  
    
end
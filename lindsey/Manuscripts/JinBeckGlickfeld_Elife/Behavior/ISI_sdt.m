function [out] = ISI_sdt(Off,Orien,output)
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



for i_off = 1:length(Off)
    for  i_orien = 1: length(Orien)
       
        if output.target.c_hit{i_off,1}(i_orien,1)<1
            [out.dprime{i_off,1}(i_orien,1),out.criterion{i_off,1}(i_orien,1)] = dprime_simple(output.target.c_hit{i_off,1}(i_orien,1),output.FA.FA(i_off,1));
        else
            signaltrialN = output.target.HT_num{i_off,1}(i_orien,1) + output.target.Miss_num{i_off,1}(i_orien,1);
            [out.dprime{i_off,1}(i_orien,1),out.criterion{i_off,1}(i_orien,1)] = dprime_simple((signaltrialN-0.5)./signaltrialN,output.FA.FA(i_off,1));
        end
        
        if output.target.S_hit{i_off,1}(i_orien,1)<1
            [out.S_dprime{i_off,1}(i_orien,1),out.S_criterion{i_off,1}(i_orien,1)] = dprime_simple(output.target.S_hit{i_off,1}(i_orien,1),output.FA.S_FA(i_off,1));
        else
            signaltrialN = output.target.S_HT_num{i_off,1}(i_orien,1) + output.target.S_Miss_num{i_off,1}(i_orien,1);
            if signaltrialN ~=0
            [out.S_dprime{i_off,1}(i_orien,1),out.S_criterion{i_off,1}(i_orien,1)] = dprime_simple((signaltrialN-0.5)./signaltrialN,output.FA.S_FA(i_off,1));
            end
        end
        
        if output.target.L_hit{i_off,1}(i_orien,1)<1
            [out.L_dprime{i_off,1}(i_orien,1),out.L_criterion{i_off,1}(i_orien,1)] = dprime_simple(output.target.L_hit{i_off,1}(i_orien,1),output.FA.L_FA(i_off,1));
        else
            signaltrialN = output.target.L_HT_num{i_off,1}(i_orien,1) + output.target.L_Miss_num{i_off,1}(i_orien,1);
            if signaltrialN ~=0
            [out.L_dprime{i_off,1}(i_orien,1),out.L_criterion{i_off,1}(i_orien,1)] = dprime_simple((signaltrialN-0.5)./signaltrialN,output.FA.L_FA(i_off,1));
            end
        end
        
        
    end
    
    
    
    
    
end
end
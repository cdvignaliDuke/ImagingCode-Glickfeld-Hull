function [output] = Contrast_HR(Contrast,Success,Miss,tContrast,RT_HM)


      for i_contrast = 1: length(Contrast)
          % find all the corrects on Contrasttation changes, ignore the cycle
          % number
             
          output.HT_num(i_contrast,1) = sum(Success &(tContrast==Contrast(i_contrast)));
          output.Miss_num(i_contrast,1) = sum (Miss&(tContrast==Contrast(i_contrast)));
          output.HT_rate(i_contrast,1) = output.HT_num(i_contrast,1) ./ (output.HT_num(i_contrast,1)+output.Miss_num(i_contrast,1));

          % use the newly corrected reaction time
          [a,b]= binofit(output.HT_num(i_contrast,1),(output.HT_num(i_contrast,1)+output.Miss_num(i_contrast,1)));
          output.confi(i_contrast,1:2) = b;
          output.c_hit(i_contrast,1) = a;
        
          % RT_on hit for all trials 
          output.RTonHit{i_contrast,1} = RT_HM(Success &(tContrast==Contrast(i_contrast)));
          output.RTonHit_mean (1,i_contrast) = mean(output.RTonHit{i_contrast,1});
          output.RTonHit_mean (2,i_contrast) = std(output.RTonHit {i_contrast,1})./sqrt(output.HT_num(i_contrast,1));
          
          
         
      end
      
    
 


end
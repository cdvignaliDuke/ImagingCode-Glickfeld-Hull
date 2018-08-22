function [Y_yeses, Y_targets] = getStimAndBehaviorYs(trOut)

[~,n]=size(trOut);
   
   k=0;
   clear Y_yeses Y_targets
   for j=1:n
       k=k+1;
       if(trOut{j}=='h')
           Y_yeses(k,1)=1;
           Y_targets(k,1)=1;
       elseif(trOut{j}=='fa')
           Y_yeses(k,1)=1;
           Y_targets(k,1)=0;
       elseif(trOut{j}=='m')
           Y_yeses(k,1)=0;
           Y_targets(k,1)=1;
       elseif(trOut{j}=='cr')
           Y_yeses(k,1)=0;
           Y_targets(k,1)=0;
       end
   end
end
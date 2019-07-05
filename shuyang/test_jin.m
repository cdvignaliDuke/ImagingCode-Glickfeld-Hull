avedfOvFrun = zeros(1,size(deriv,2));
for i = 1: size(deriv,2) %for each cell
avedfOvFrun(i) = mean(dfOvF(frm_run,i));
end 

avedfOvFstay = zeros(1,size(deriv,2));
for i = 1: size(deriv,2) %for each cell
avedfOvFstay(i) = mean(dfOvF(frm_stay,i));
end 
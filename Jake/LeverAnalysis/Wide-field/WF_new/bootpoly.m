
% fita = bootstrp(1000, @bootfun, x, y);
% 
% cia = bootci(1000,@bootfun,x, y);
diff_fit = []; nROI = length(totCorrS);
for i = 1:1000
    ind = randsample(nROI, 30, 'true');

    diff_fit = [diff_fit; polyfit(totCorrS(ind,2), totCorrS(ind,1), 1) - polyfit(totCorrF(ind,2), totCorrF(ind,1), 1)];
end
x = diff_fit(:,2);
SEM = std(x)/sqrt(length(x));               % Standard Error
ts = tinv([0.025  0.975],length(x)-1);      % T-Score
CI = mean(x) + ts*SEM;                      % Confidence Intervals


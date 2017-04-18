function [A, tau, R_square] = singExpDecayFit(timeCourse);

t = 1:length(timeCourse);
A_guess = mean(timeCourse(1:5,:),1);
ss = mean(timeCourse(length(t)-5:length(t),:),1);
tau_guess = -1*find(timeCourse<(ss+ (0.4*(A_guess-ss))),1,'first');

options = optimset('MaxFunEvals',inf,'MaxIter',100000);
[out,fval,success] = fminsearch(@singExpDecay_sse, [A_guess, tau_guess],options);
if success == 1
    A=out(1);
    tau=out(2);
    sse = fval;
    sse_tot = sum((timeCourse - mean(timeCourse)).^2);
    R_square = 1-(sse/sse_tot);
else
    A=NaN;
    tau=NaN;
    sse = NaN;
    sse_tot = NaN;
    R_square = NaN;
end
    
    %nested function
    function temp_sse = singExpDecay_sse(in)
        A_tmp = in(1);
        tau_tmp = in(2);
       
        y_fit = A_tmp.*exp(t'./tau_tmp);
        residuals = timeCourse - y_fit;
        temp_sse = sum(residuals.^2);
        
    end
end
function [A1, A2, tau1, tau2, R_square] = doubExpDecayFit(timeCourse);

t = length(timeCourse);
A1_guess = mean(timeCourse(1:5,:),1);
tau1_guess = -1*find(timeCourse<(0.4*A1_guess),1,'first');
A2_guess = mean(timeCourse(t-5:t,:),1);
tau2_guess = -.5*t;

options = optimset('MaxFunEvals',inf,'MaxIter',inf);
[out,fval,success] = fminsearch(@singExpDecay_sse, [A1_guess, tau1_guess, A2_guess, tau2_guess],options);
if success == 1
    A1=out(1);
    tau1=out(2);
    A2=out(3);
    tau2=out(4);
    sse = fval;
    sse_tot = sum((timeCourse - mean(timeCourse)).^2);
    R_square = 1-(sse/sse_tot);
end
    
    %nested function
    function temp_sse = singExpDecay_sse(in)
        A1_tmp = in(1);
        tau1_tmp = in(2);
        A2_tmp = in(3);
        tau2_tmp = in(4);
       
        y_fit = A1_tmp.*exp((1:t)'./tau1_tmp) + A2_tmp.*exp((1:t)'./tau1_tmp);
        residuals = timeCourse - y_fit;
        temp_sse = sum(residuals.^2);
        
    end
end
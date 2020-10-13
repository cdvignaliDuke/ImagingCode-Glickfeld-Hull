function [tao,A,sse,R_square] = SingleExponFit(time,Response)


tao_guess=100;


A_guess = 1;





options = optimset('MaxFunEvals',1e8,'MaxIter',1e8);
[out,fval,success] = fminsearch(@Expon_ssen, [tao_guess,A_guess],options);
if success == 1
    
    tao=out(1);
    
    
    A = out(2);
   
    sse=fval;
    sse_tot = sum((Response - mean(Response)).^2);
    R_square = 1-(sse/sse_tot);
    
end


    % define the nested subfunction
    function miaosse = Expon_ssen(in)
        % pull out the slope and intercept
        
        tao_tmp = in(1);       
        A_tmp = in(2);
    
        
       
       
        y_fit = 1-A_tmp.*exp(time./(-tao_tmp));
        
        y_fit_zero = 1-A_tmp;
        
        residuals = Response - y_fit;
        miaosse = sum(residuals.^2);
        
        if tao_tmp < 0 || A_tmp<0 ||  y_fit_zero<0
            miaosse = inf;
        end
        
    end



end
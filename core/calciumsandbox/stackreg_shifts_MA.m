%STACKREG_SHIFTS_MA Stack registration using interp2 
%
%stack_reg = stackreg_shifts_MA(stack, shifts);
%
%where stack is a 3D FastLoadStackAK-loaded tiffstack, and shifts is the
%output of RegMulticore_shifts.m (Nframes x 2);
%080503 MA

    function stack_reg = stackreg_shifts_MA(stack, shifts);
    Lx = size(stack,1);
    Ly = size(stack,2);
    Nfr = size(stack,3);
    X = [1:Lx]'*ones(1,Ly);
    Y = ([1:Ly]'*ones(1,Lx))';
    
    stack_reg = zeros(size(stack));
    for count = 1:Nfr
        if round(count/1000)*1000 == count
            fprintf('%s\n',num2str(count));
        end
        X2 = X + shifts(count,2);
        Y2 = Y + shifts(count,1);
       stack_reg(:,:,count) = interp2(Y,X,stack(:,:,count),Y2,X2,'*linear');
    end
    
    
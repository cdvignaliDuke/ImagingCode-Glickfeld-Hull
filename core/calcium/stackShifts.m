% STACKSHIFTS Stack registration using interp2
%  STACK_REG = STACKSHIFTS(STACK, SHIFTS);
%  where STACK is 3d stack array, 
%  and SHIFTS is the output of REGMULTICORE_SHIFTS (Nframes x 2);
%
% To register using shifts from STACKREGISTER:
%  OUTS = STACKREGISTER(STACK,TARGET)
%  STACK_REG = STACKSHIFTS(STACK,-OUTS(:,4:-1:3))
%
%  STACK_REG = STACKSHIFTS(STACK, SHIFTS, INTFLAG) if INTFLAG == 1
%   constrains shifts to integers.  
% 
% by Mark Andermann 2008
%
% see STACKREGISTER, REGMULTICORE_SHIFTS
%
% 09/11/13 K Ohki added INTFLAG
% 10/05/04 V Bonin values outside range (extrapval) set to 0

function stack_reg = stackShifts(stack, shifts, int_shift);

if nargin <3
    int_shift=0;
end

[Lx, Ly, Nfr] = size(stack);

if int_shift==1
    for count = 1:Nfr
        stack_reg(:,:,count)=circshift(stack(:,:,count),round(shifts(count,:)));
    end
    return
end

extrapval = 0; % pad with zeros outside interpolation range

X = [1:Lx]'*ones(1,Ly);
Y = ([1:Ly]'*ones(1,Lx))';

stack_reg = zeros(size(stack));
for count = 1:Nfr
    if round(count/1000)*1000 == count
        fprintf('%s\n',num2str(count));
    end
    X2 = X + shifts(count,2);
    Y2 = Y + shifts(count,1);
    stack_reg(:,:,count) = interp2(Y,X,stack(:,:,count),Y2,X2,'*linear',extrapval);
end

return;

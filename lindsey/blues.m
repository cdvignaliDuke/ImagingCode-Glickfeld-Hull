function x = blues(num);
% create a color map of blue: pale to dark
x = zeros(num,3);
x(:,3) = 1;
x(:,1) = 1-logspace(log(1/num),log(1-(1/num)),num);
x(:,2) = 1-logspace(log(1/num),log(1-(1/num)),num);

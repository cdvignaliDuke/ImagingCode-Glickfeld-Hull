function x = grays(num);
% create a color map of gray: white to black
x = zeros(num,3);
x(:,3) = 1-logspace(log(1/num),log(1-(1/num)),num);
x(:,1) = 1-logspace(log(1/num),log(1-(1/num)),num);
x(:,2) = 1-logspace(log(1/num),log(1-(1/num)),num);

function x = reds(num);
% create a color map of reds: pink to red
x = zeros(num,3);
x(:,1) = 1;
x(:,2) = 1-logspace(log(1/num),log(1-(1/num)),num);
x(:,3) = 1-logspace(log(1/num),log(1-(1/num)),num);


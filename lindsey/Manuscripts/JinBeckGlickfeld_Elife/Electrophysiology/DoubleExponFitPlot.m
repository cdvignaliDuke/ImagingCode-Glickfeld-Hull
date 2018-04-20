function y=DoubleExponFitPlot(time,tao1,tao2,ft,A1,A2,color)

y = ft-A1.*exp(time./(-tao1))-A2.*exp(time./-tao2);

plot(time,y,'color',color,'linewidth',1.5)

end    
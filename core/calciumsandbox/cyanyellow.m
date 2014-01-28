function c = cyanyellow

top = [zeros(32,1),1*linspace(1,0,32)',1*linspace(1,0,32)'];
bottom = [1*linspace(0,1,32)',1*linspace(0,1,32)',zeros(32,1)];
bottom(1,:)=[];
c = [top;bottom];

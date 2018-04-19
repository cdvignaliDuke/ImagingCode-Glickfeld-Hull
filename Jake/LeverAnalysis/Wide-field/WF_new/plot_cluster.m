 [class,type] = dbscan(data,3,3);
for k = 1:length(data)
    if(type(k)==1)
       plot(data(k,1),data(k,2),'ks');hold on;
    elseif(type(k) ==0)
            plot(data(k,1),data(k,2),'bo');hold on;
    else
            plot(data(k,1),data(k,2),'rx');hold on;
    end
    text(data(k,1),data(k,2),num2str(class(k)))
   
end
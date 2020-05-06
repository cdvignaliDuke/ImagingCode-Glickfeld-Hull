
function [dunns]=ioclustervalidity2_min(allclusters, eventsmat, combclusters);


numclusters=size(allclusters,1);
cluster = cell(1,numclusters);
if combclusters~=0;
    a=allclusters{combclusters(1),1};
    b=allclusters{combclusters(2),1};
    c=[a b];
    c=sort(c);
    e = unique(c);    
    cluster{combclusters(1)} = e;
end

cnt=1;
for i=1:numclusters;
    if combclusters~=0;
        if cnt==combclusters(1);
            cnt=cnt+1;
        end
        if i~=combclusters(1) && i~=combclusters(2);
            cluster{cnt} = allclusters{i,1};
            cnt=cnt+1; 
        end
    else
        cluster{cnt} = allclusters{i,1};
        cnt=cnt+1; 
    end          
end

if combclusters~=0;
    numclusters=numclusters-1;
end


maxdiam=0;

for p=1:numclusters;
    clmembers=cluster{p};
    numpnts=length(clmembers);
    if numpnts~=1;
        for r=1:numpnts;
            for s=1:numpnts;
                %pairdist=(sum((eventsmat(clmembers(r),:)-eventsmat(clmembers(s),:)).^2))^0.5;
                pairdist=1-corr2(eventsmat(clmembers(r),:),eventsmat(clmembers(s),:));
                maxdiam=max(maxdiam,pairdist);
            end
        end
    end
end

mind=10000;

for p=1:numclusters;
    for r=1:numclusters;
        if p~=r;
             cl1members=cluster{p};
             cl2members=cluster{r};
             numpnts1=size(cl1members,2);
             numpnts2=size(cl2members,2);
             for s=1:numpnts1;
                for tt=1:numpnts2;
                    pairdist=1-corr2(eventsmat(cl1members(s),:),eventsmat(cl2members(tt),:));
                    if pairdist~=0;  
                        mind=min(mind,pairdist);
                    end
                end
             end
        end
     end
end

dunns = mind/maxdiam;






                            
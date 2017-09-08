
function [newclusters, dunns, centroidcorr, dendmem]=ioclustervalidity2(allclusters, eventsmat, combclusters)

numclusters = size(allclusters,1);

if combclusters~=0;
    a=allclusters{combclusters(1),1};
    b=allclusters{combclusters(2),1};
    newclusters{combclusters(1), 1} = unique([a b]);
%     eval(['cluster' num2str(combclusters(1)) '=e;']);
end



cnt=1;
for i=1:numclusters;
    if combclusters~=0;
        if cnt==combclusters(1);
            cnt=cnt+1;
        end
        if i~=combclusters(1) && i~=combclusters(2);
%             eval(['cluster' num2str(cnt) '=allclusters{' num2str(i) ',1};']);
            newclusters{cnt, 1} = allclusters{i, 1};
            cnt=cnt+1; 
        end
    else
%         eval(['cluster' num2str(cnt) '=allclusters{' num2str(i) ',1};']);
        newclusters{cnt, 1} = allclusters{i, 1};
        cnt=cnt+1; 
    end          
end

if combclusters~=0;
    numclusters=numclusters-1;
end


centroidmat=[];
dendmem=zeros(size(eventsmat,1),numclusters);

% for i=1:numclusters;
%     eval(['aa=cluster' num2str(i) ';']);
%     nummembers=length(aa);
%     centroid=zeros(1,numtimes);
%     for j=1:nummembers;
%         centroid(1,1:numtimes)=eventsmat(aa(j),1:numtimes)+centroid(1,1:numtimes);
%     end
%     eval(['centroid' num2str(i) '=centroid/nummembers;']);    
%     centroidmat=[centroidmat; centroid/nummembers];
% %     parfor ii=1:size(eventsmat,1);
% %         dendmem(ii,i)=corr2(eventsmat(ii,:),centroid/nummembers);
% %     end
%         
% end

for i = 1:numclusters
    centroid = sum(eventsmat(allclusters{i,1},:),1) / length(allclusters{i,1});
    newclusters{i,2} = centroid;
    centroidmat=[centroidmat; centroid];
end

centroidcorr = corrcoef(centroidmat');

maxdiam=0;

mind=10000;


parfor p=1:numclusters;
    
    cl1members = newclusters{p,1};
    numpnts=length(cl1members);
    if numpnts~=1;
        pairdist = 1 - corr(eventsmat(cl1members,:)',eventsmat(cl1members,:)');
        maxdiam = max(maxdiam, max(max(pairdist)));
    end
    for r=1:numclusters;
        if p~=r;
%              eval(['cl1members=cluster' num2str(p) ';']);
%              eval(['cl2members=cluster' num2str(r) ';']);
             cl2members = newclusters{r,1};
             pairdist = 1 - corr(eventsmat(cl1members,:)',eventsmat(cl2members,:)');
             mind = min(mind, min(pairdist(pairdist > 0)));
        end
     end
end

dunns = mind/maxdiam;


% newclusters=cell(numclusters,2);
% 
% for i=1:numclusters;
%     eval(['newclusters{' num2str(i) ',1}=cluster' num2str(i) ';']);
%     eval(['newclusters{' num2str(i) ',2}=centroid' num2str(i) ';']);
% end

                            
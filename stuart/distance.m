% calculate distance of cells in mask using mean position

mx=zeros(nCells,1);
my=mx;
for ii=1:nCells
    [x,y]=find(mask_cell{1}==ii);
    mx(ii)=mean(x);
    my(ii)=mean(y);
end

distance=size(nCells);
for ii=1:nCells
    for jj=1:nCells
        distance(ii,jj)=sqrt((mx(jj)-mx(ii))^2+(my(jj)-my(ii))^2);
    end
end


%%
distance(distance==0)=nan;
figure; hold on
for ii=1:nCells
    scatter(distance(ii,:),heatu(ii,:),'b')
end
title('Unsuccessful')
xlabel('Distance (pixels)')
ylabel('Correlation Strength')
ylim([-.3 .6])

%%
figure; hold on
for ii=1:nCells
    for jj=1:nCells
        scatter(distance(ii,jj),mean([auclist(ii),auclist(jj)]))
    end
end
ylabel('AUC')
xlabel('Distance (pixels)')

    
    
    
    
    
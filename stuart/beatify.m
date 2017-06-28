function [aresp,leftresp,rightresp]=beatify(nCells,tLeftTrial,data_stim_dfof)
% variance for cells across trials to gauge responsiveness, used to
% generate typical response, or canon
nTrials=length(tLeftTrial);
vc=cell(1,nCells);
leftresp=zeros(1,nCells);
rightresp=zeros(1,nCells);
aresp=zeros(1,nCells);
for ii=1:nCells
    vc{ii}=zeros(1,sum(tLeftTrial));
    nn=1;
    for jj=find(tLeftTrial) %for left trials
        vc{ii}(nn)=std(data_stim_dfof(:,ii,jj)); %get standard deviation of each left trial for each cell
        nn=nn+1;
    end
    if std(vc{ii})>.0135 %threshold of responsive or not picked out from below
        leftresp(ii)=1;
    end
    %
    vc{ii}=zeros(1,sum(~tLeftTrial)); 
    nn=1;
    for jj=find(~tLeftTrial) %for right trials, otherwise identical as above
        vc{ii}(nn)=std(data_stim_dfof(:,ii,jj));
        nn=nn+1;
    end
    if std(vc{ii})>.0135
        rightresp(ii)=1;
    end
    if and(~leftresp(ii),~rightresp(ii)) %if responsive to neither side check if general non specific trend
        vc{ii}=zeros(1,nTrials);
        for jj=nTrials
            vc{ii}(jj)=std(data_stim_dfof(:,ii,jj));
        end
        if std(vc{ii})>.0135
            aresp(ii)=1;
        end
    elseif and(leftresp(ii),rightresp(ii)) %if responsive to both mark for averaging response
        aresp(ii)=1;
    end
    %}
end
%% determine threshold for if cell is responsive, troubleshooting only
%{
for ii=1:nCells
    temp(ii,1)=mean(vc{ii});
    temp(ii,2)=std(vc{ii});
    %temp(ii,3)=std(mean(data_stim_dfof(:,ii,:),3));
end
figure; hold on
scatter(temp(temp(:,2)>.0135,1),temp(temp(:,2)>.0135,2),'b')
scatter(temp(temp(:,2)<.0135,1),temp(temp(:,2)<.0135,2),'r')
  
for ii=[24 25 29 32 52 55 57 91]
    figure; hold on
    for jj=1:nTrials
        plot(tt,data_stim_dfof(:,ii,jj))
    end
    title(num2str(ii))
    plot(tt,mean(data_stim_dfof(:,ii,:),3),'k')
end
%}
end     
%% determine discriminability between different conditions of stimulus presentations
clear N P can%clear any residual values
LvR=[0,6,12]; %index of row beginnings for different contrasts
auclist=cell(3,1);
negated=cell(3,1);
group=newhighLvR6;%1:nCells; %choose group of cells to run here
Trialside=tLeftTrial;
%generate typical trace/canon
[aresp,leftresp,rightresp]=beatify(nCells,tLeftTrial,data_stim_dfof); %function to determine response
%figure;hold on
cc=1; %counter of appropriate conditions
for jj=1:nCells
    if aresp(jj)
        can(:,cc)=mean(data_stim_dfof(:,jj,tGratingContrast==1),3); %if cell responds to both, take all high contrast
        % plot(tt,can(:,cc),'b')
        cc=cc+1;
    elseif leftresp(jj)
        can(:,cc)=mean(data_stim_dfof(:,jj,and(tGratingContrast==1,tLeftTrial)),3); %if one side only take just that side
        % plot(tt,can(:,cc),'r')
        cc=cc+1;
    elseif rightresp(jj)
        can(:,cc)=mean(data_stim_dfof(:,jj,and(tGratingContrast==1,~tLeftTrial)),3);
        %  plot(tt,can(:,cc),'g')
        cc=cc+1;
    end
end
canon=mean(can,2); %all together the responses above average to the typical response
% plot(tt,canon,'k')
% plot(tt,mean(can(:,setdiff(1:79,[40,53])),2),'k')
for kk=LvR(1:2)
    auclist{kk/6+1}=zeros(1,length(group));
    negated{kk/6+1}=zeros(1,length(group));
    %figure
    P=cell(1,length(group));
    N=cell(1,length(group));
    cond1=zeros(50,1);
    cond2=zeros(50,1);
    Trialside=Trialside*-1+1; %flips side the binary is true for
    for jj=group
        %figure; hold on
        %cond1=zeros(50,length(idx{kk+1})+length(idx{kk+2})+length(idx{kk+3}));
        %cond2=zeros(50,length(idx{kk+4})+length(idx{kk+5})+length(idx{kk+6}));
        ll=1; rr=1;
        for ii=1:nTrials
            %if and(tLeftTrial(ii),tGratingContrast(ii)==cons(kk/6+1)) %left/ipsi at increaseing contrasts -or-
            if and(Trialside(ii),tGratingContrast(ii)==cons(1)) %First if contra and 6p, then if ipsi and 6p
                cond1(:,ll)=data_stim_dfof(:,jj,ii);
                ll=ll+1;
            end
        end
        for ii=1:nTrials
            %if and(~tLeftTrial(ii),tGratingContrast(ii)==cons(kk/6+1)) %right/contra at increasing contrasts -or-
            if and(Trialside(ii),tGratingContrast(ii)==cons(2)) %first if contra and 10p, then if ipsi and 10p
                cond2(:,rr)=data_stim_dfof(:,jj,ii);
                rr=rr+1;
            end
        end
        %compare canon to traces
        for ii=1:size(cond1,2)
            trace=cond1(:,ii); %take each trial one at a time 
            [c, lgs]=xcov(canon(21-buf:21+buf), trace(21-buf:21+buf), 'coeff'); %simple dot product here instead of whole correlation?
            N{find(group==jj)}(ii)=mean(c(find(lgs==0)-10:find(lgs==0)+10));
        end
        for ii=1:size(cond2,2)
            trace=cond2(:,ii);
            [c,lgs]=xcov(canon(21-buf:21+buf), trace(21-buf:21+buf), 'coeff');
            P{find(group==jj)}(ii)=mean(c(find(lgs==0)-10:find(lgs==0)+10));
        end
        %
        %generate PDF %ROC curve
        criterion= -1:.01:1;
        TP=zeros(size(criterion));
        FP=TP;
        for ii=1:length(criterion)
            TP(ii)=sum(P{find(group==jj)}>criterion(ii))/length(P{find(group==jj)}); %acurate condition id based on criterion
            FP(ii)=sum(N{find(group==jj)}>criterion(ii))/length(N{find(group==jj)}); %inaccurate " 
        end
        
        %AUC (if >.5 flip classifier sign)
        tbl=table(TP',FP');
        tbl=unique(tbl); %eliminate duplicates to graph and get area
        Y=table2array(tbl(:,1));
        X=table2array(tbl(:,2));
        auc=trapz(X,Y); %area under curve
        if auc<0.5
            auc=1-auc; %if lesss than .5 still has information just applied wrong, so essentially switch to true and false negatives
            %plot(TP,FP,'r')
            negated{kk/6+1}(jj)=1; 
        else
            %plot(FP,TP,'b')
        end
        
        auclist{kk/6+1}(jj)=auc;
        %{
        lst=min([P{jj},N{jj}]):min([range(P{jj}),range(N{jj}),.06])/3:max([P{jj},N{jj}])+.05; %set histogram width 
        figure
        hold on
        histogram(P{jj},lst)
        histogram(N{jj},lst)
        title([num2str(jj),' ',num2str(cons(kk*1/6+1))])
        %}
    end
    %std of each cell
    %each trial one response weighted with the cells' std
    %histogram of weighted correlations from pop
    %auc from that, save pop auc for different contrasts
    %
    %P{1}=P{48};
    %N{1}=N{48};
    Pwi=zeros(1,length(group));
    Nwi=zeros(1,length(group));
    for nc=1:length(group)
        Pwi(nc)=1/std(P{nc})^2; %1/standard deviation of correlations squared
        Nwi(nc)=1/std(N{nc})^2;
    end
    waP=-1*ones(1,length(P{1}));
    waN=-1*ones(1,length(N{1}));
    for ii=1:length(P{1})
        wa=0;
        for jj=1:length(group)
            wa=wa+P{jj}(ii)*Pwi(jj); %sum of weighted responses
        end
        waP(ii)=wa/sum(Pwi); %sum divided by sum of weights
    end
    for ii=1:length(N{1})
        wa=0;
        for jj=1:length(group)
            wa=wa+N{jj}(ii)*Nwi(jj);
        end
        waN(ii)=wa/sum(Nwi);
    end
    %lst=min([waP,waN]):min([range(waP),range(waN),.15])/3:max([waP,waN])+.05;
    figure; hold on
    histogram(waP,lst)
    histogram(waN,lst)
    plot([mean(waP) mean(waP)],[0 25],'b')
    plot([mean(waN) mean(waN)],[0 25],'r')
    title(num2str(cons(kk/6+1)))
    
    TP=zeros(size(criterion));
    FP=TP;
    for ii=1:length(criterion)
        TP(ii)=sum(waP>criterion(ii))/length(waP);
        FP(ii)=sum(waN>criterion(ii))/length(waN);
    end
    tbl=table(TP',FP');
    tbl=unique(tbl);
    Y=table2array(tbl(:,1));
    X=table2array(tbl(:,2));
    auc=trapz(X,Y);
    if auc<.5
        auc=1-auc;
        subp_negated(kk/6+1)=1;
    else
        subp_negated(kk/6+1)=0;
    end
    subp_auclist(kk/6+1)=auc;
    %}
    
    %{
    plot([0 1],[0 1],'k','LineWidth',3)
    title('ROC Curve')
    xlabel('False Positive (Blue) True Negative (Red)')
    ylabel('True Positive (Blue) False Negative (Red)')
%
%rg=range(N)/3;
lst=min([P,N]):min([range(P),range(N),.06])/3:max([P,N]);
figure
hold on
histogram(P,lst)
histogram(N,lst)
xlabel('')
title([num2str(jj),' ',num2str(cons(kk*1/6+1))])
xlim([-.2 .7]);ylim([0 25])
xlabel('Orange==Ipsi     Blue==Contra')
    %}
    %{
   if kk<=1
      %t12=P;
      t22=N;%(23:29)
      %t32=waP;
      t42=waN;
      %t52=Pwi;
      t62=Nwi;
   end
    %}
end
subp_auclist
figure
scatter([1 2 3],subp_auclist)
xlim([0.5 3.5])
ylim([.4 1])
ax=gca;
ax.XTick=[1 2 3];
ax.XTickLabel={'6%','10%','100%'};
xlabel('Contrast')
ylabel('AUC')
title('Contra vs Ipsi')

xlabel('side')
xlim([.5 2.5])
ax.XTickLabel={'Contra', 'Ipsi'};
title('10 v 6')
%%
figure
hold on
for ii=1:3
    yup=ii*ones(nCells,1);
    scatter(yup,auclist{ii},'b')
    scatter(yup(1), mean(auclist{ii}),500,'+k')
    scatter(yup(negated{ii}==1),auclist{ii}(negated{ii}==1),'r')
end
xlim([0.5 3.5])
ylim([.4 1])
ax=gca;
ax.XTick=[1 2 3];
ax.XTickLabel={'6%','10%','100%'};
xlabel('Contrast')
ylabel('AUC')
title('Contra vs Ipsi')

[h,p]=ttest(auclist{3},auclist{1}, 'tail', 'right');

%% go through cells at lowest contrast and see who can discriminate, then maybe look at how correlated they are vs the ones who can't discriminate as well, correlation strength/discrimination strength plot, maybe also distance vs both
% options for calculating discriminabiity of pair: average each one;
% treat/combine both as single unit and try and discriminate, essentially
% twice as many trials to plot
%{
UvS=[0,6,12];
auclist=cell(3,1);
negated=cell(3,1);

buf=10;

for kk=UvS
    auclist{kk*1/6+1}=zeros(1,nCells);
    negated{kk*1/6+1}=zeros(1,nCells);
    figure; hold on
    for jj=1:nCells %
        %get low contrast combined contra and ipsi success vs fail
        suc6=zeros(50,length(idx{kk+1})+length(idx{kk+4}));
        un6=zeros(50,length(idx{kk+2})+length(idx{kk+3})+length(idx{kk+5})+length(idx{kk+6}));
        s6=1; u6=1;
        for ii=1:nTrials
            if and(SIx(ii),tGratingContrast(ii)==cons(kk*1/6+1))
                suc6(:,s6)=data_stim_dfof(:,jj,ii);
                s6=s6+1;
            end
        end
        for ii=1:nTrials
            if and(~SIx(ii),tGratingContrast(ii)==cons(kk*1/6+1))
                un6(:,u6)=data_stim_dfof(:,jj,ii);
                u6=u6+1;
            end
        end
        
        %compare canon to traces
        clear N P
        
        for ii=1:nTrials
            can(:,ii)=data_stim_dfof(:,jj,ii);
        end
        canon=mean(can,2);
        
        for ii=1:size(suc6,2)
            trace=suc6(:,ii);
            [c, lgs]=xcov(canon(21-buf:21+buf), trace(21-buf:21+buf), 'coeff');
            P(ii)=mean(c(find(lgs==0)-10:find(lgs==0)+10));
        end
        %canon=mean(un6,2);
        for ii=1:size(un6,2)
            trace=un6(:,ii);
            [c,lgs]=xcov(canon(21-buf:21+buf), trace(21-buf:21+buf), 'coeff');
            N(ii)=mean(c(find(lgs==0)-10:find(lgs==0)+10));
        end
        
        %generate PDF %ROC curve
        criterion= -1:.01:1;
        TP=zeros(size(criterion));
        FP=TP;
        for ii=1:length(criterion)
            TP(ii)=sum(P>criterion(ii))/length(P);
            FP(ii)=sum(N>criterion(ii))/length(N);
        end
        
        %AUC (if >.5 flip classifier sign)
        tbl=table(TP',FP');
        tbl=unique(tbl);
        Y=table2array(tbl(:,1));
        X=table2array(tbl(:,2));
        auc=trapz(X,Y);
        if auc<0.5
            auc=1-auc;
            plot(TP,FP,'r')
            negated{kk*1/6+1}(jj)=1;
        else
            plot(FP,TP,'b')
        end
        
        auclist{kk*1/6+1}(jj)=auc;
        
    end
    plot([0 1],[0 1],'k','LineWidth',3)
    title('ROC Curve')
    xlabel('False Positive (Blue) True Negative (Red)')
    ylabel('True Positive (Blue) False Negative (Red)')
    
%{
%rg=range(N)/3;
lst=min([P,N]):min([range(P),range(N),.15])/3:max([P,N]);
figure
hold on
histogram(P,lst)
histogram(N,lst)
xlabel('')
title([num2str(jj),' ',num2str(cons(kk*1/6+1))])
%xlim([-.2 .7])
%ylim([0 9])
%}
end
%%
%nCells=nCells;
%highcells=highcells;
figure
hold on
for ii=1:3
    yup=ii*ones(nCells,1);
    scatter(yup,auclist{ii},'b')
    scatter(yup(1), mean(auclist{ii}),500,'+k')
    
    scatter(yup(negated{ii}==1),auclist{ii}(negated{ii}==1),'r')
end
xlim([0.5 3.5])
ylim([.4 1])
ax=gca;
ax.XTick=[1 2 3];
ax.XTickLabel={'6%','10%','100%'};
xlabel('Contrast')
ylabel('AUC')
title('Successful vs Unsuccessful')

[h,p]=ttest(auclist{2},auclist{3}, 'tail', 'right');
%}

%% hist PDF
lst=min([P,N]):min([range(P),range(N),.15])/3:max([P,N]);
figure
hold on
histogram(P,lst)
histogram(N,lst)
xlabel('')

%% heatmaps, combine success and fail, get rid of autos for distance/auc plots
%{
all=[];
%all low contrast trials dfof
%xcov
%assign value to heat array
%figure


%% heatmaps of low contrast combined ipsi and contra, success vs fail
%also look at differences
heats=nan(nCells);
heatu=heats;

suc=[];
for ii=[idx{1},idx{4}] %success low contrast left and right
    suc=[suc;data_stim_dfof(:,:,ii)];
end
[c,lgs]=xcov(suc,25,'coeff');

figure
hold on
[ui, iu]=find(and(c(26,:)>.2,c(26,:)<.9));
for ii=unique(iu)
    plot(lgs,c(:,ii))
end

for ii=1:nCells
    for jj=1:nCells
        heats(ii,jj)=mean(c(24:28,nCells*(ii-1)+jj));%take 2 frames around time 0 and 0 frame
    end
end
figure
heatmap(heats,1:nCells,1:nCells,[],'Colorbar',true,'MinColorValue',-.6,'MaxColorValue',1)
xlabel('Cell');ylabel('Cell');
title('Successful Trials')

figure
un=[];
for ii=[idx{2},idx{3},idx{5},idx{6}] %all fails both sides
    un=[un;data_stim_dfof(:,:,ii)];
end
[c,lgs]=xcov(un,25,'coeff');
for ii=1:nCells
    for jj=1:nCells
        heatu(ii,jj)=mean(c(24:28,nCells*(ii-1)+jj));
    end
end
heatmap(heatu,1:nCells,1:nCells,[],'Colorbar',true,'MinColorValue',-.6,'MaxColorValue',1)
xlabel('Cell');ylabel('Cell');
title('Unsuccessful Trials')

figure
heatmap((heats-heatu),1:nCells,1:nCells,[],'Colorbar',true,'MinColorValue',-.6,'MaxColorValue',1)
xlabel('Cell');ylabel('Cell');
title('Difference between Successful and Unsucessful')

figure
for ii=1:nCells
    for jj=1:nCells
        subplot(1,3,1); hold on
        scatter(mean([auclist(ii),auclist(jj)]),heats(ii,jj),pointsize,distance(ii,jj))
        subplot(1,3,2); hold on
        scatter(mean([auclist(ii),auclist(jj)]),heatu(ii,jj),pointsize,distance(ii,jj))
        subplot(1,3,3); hold on
        scatter(mean([auclist(ii),auclist(jj)]),heats(ii,jj)-heatu(ii,jj),pointsize,distance(ii,jj))
    end
end

subplot(1,3,1)
title('Success Responses Colored by Distance')
xlabel('AUC')
ylabel('Correlation Strength')
subplot(1,3,2)
title('Unsuccessful Responses Colored by Distance')
xlabel('AUC')
ylabel('Correlation Strength')
subplot(1,3,3)
title('Response Difference Colored by Distance')
xlabel('AUC')
ylabel('Correlation Strength')
ylim([-.4 1])

%}



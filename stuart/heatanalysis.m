%% heat analysis of population correlations, run stu_singleChannelTC_2AFC.m first 

tLeftTrial = celleqel2mat_padded(input.tLeftTrial); %binary of if left trial in vector form, 143
%remember frames are taken at 30Hz, 50 frames = 1665 ms
data_stim_dfof; %might have to regenerate from stu_singleChannelTC_2AFC
data_choice_dfof;

SIx = strcmp(input.trialOutcomeCell, 'success');%success index, 182
MIx = strcmp(input.trialOutcomeCell, 'ignore');%miss/ignore index, 11
FIx = strcmp(input.trialOutcomeCell, 'incorrect');%fail index, 51

tGratingContrast = celleqel2mat_padded(input.tGratingContrast);%vector of contrasts
cons = unique(tGratingContrast);%possible contrast values
ncon = length(cons);

idx=cell(3,6);
for ii=1:3 %cells of conditions in a 3 by 6 grid for contrast, outcome, and side
idx{(1+(ii-1)*6)}=intersect(intersect(find(tLeftTrial),find(SIx)),find(tGratingContrast ==cons(ii)));
idx{(2+(ii-1)*6)}=intersect(intersect(find(tLeftTrial),find(MIx)),find(tGratingContrast ==cons(ii)));
idx{(3+(ii-1)*6)}=intersect(intersect(find(tLeftTrial),find(FIx)),find(tGratingContrast ==cons(ii)));
idx{(4+(ii-1)*6)}=intersect(intersect(find(~tLeftTrial),find(SIx)),find(tGratingContrast ==cons(ii)));
idx{(5+(ii-1)*6)}=intersect(intersect(find(~tLeftTrial),find(MIx)),find(tGratingContrast ==cons(ii)));
idx{(6+(ii-1)*6)}=intersect(intersect(find(~tLeftTrial),find(FIx)),find(tGratingContrast ==cons(ii)));
end
ttl{1}='Success Left';
ttl{2}='Ignore Left';
ttl{3}='Incorrect Left';
ttl{4}='Success Right';
ttl{5}='Ignore Right';
ttl{6}='Incorrect Right';
tcon=cell(3,1);
for ii=1:3
    tcon{ii}=num2str(cons(ii));
end
lth=length(data_stim_dfof(1,:,1));
heat=nan(lth);
%tst=nan(length(idx{nn})*50,lth);

for ii=1:6 %just to see trial numbers, need to do three times 
lol(3,ii)=length(idx{ii+12});
end

%% plots
%for after stim on
%use last 11 frames for response (40-50) (referenced to stim)
figure
for nn=1:18 %loop through conditions
    if idx{nn}
        tst=[];
        for ii=idx{nn}
            tst=[tst;data_stim_dfof(:,:,ii)];
        end
        
        [c,lgs]=xcov(tst,25,'coeff');
      
        for ii=1:lth
            for jj=1:lth
                heat(ii,jj)=mean(c(24:28,lth*(ii-1)+jj));%take 2 frames after time 0 and 0 frame
            end
        end
        subplot(3,6,nn)
        heatmap(heat,1:lth,1:lth,[],'Colorbar',true,'MinColorValue',-.6,'MaxColorValue',1)
        %maxs=max(maxs,max(max(heat)));
        %mins=min(mins,min(min(heat)));
    else continue
    end
end
%suptitle('stim')

%for after choice
figure
for nn=1:18
    if idx{nn}
        tst=[];
        for ii=idx{nn}
            tst=[tst;data_choice_dfof(:,:,ii)];
        end
        
        [c,lgs]=xcov(tst,25,'coeff');
        
        for ii=1:lth
            for jj=1:lth
                heat(ii,jj)=mean(c(24:28,lth*(ii-1)+jj));
            end
        end
        subplot(3,6,nn)
        heatmap(heat,1:lth,1:lth,[],'Colorbar',true,'MinColorValue',-.6,'MaxColorValue',1)
        %maxc=max(maxc,max(max(heat)));
        %minc=min(minc,min(min(heat)));
    else continue
    end
end
%suptitle('choice')
%% labels
for ii=1:6
    subplot(3,6,ii)
    title(ttl{ii})
end
xx=[1,7,13];
for ii=1:3
    subplot(3,6,xx(ii))
    ylabel([tcon{ii},'           '],'FontWeight','bold','rot',0)
end


%% sanity plots of correlograms
figure
for ii=1:16
    subplot(4,4,ii)
    if ii<=4
        plot(lgs,c(:,ii))
    elseif ii<=8
        plot(lgs,c(:,lth+ii-4))
    elseif ii<=12
        plot(lgs,c(:,lth*2+ii-8))
    elseif ii<=16
        plot(lgs,c(:,lth*3+ii-12))
    end
    ylim([0 1])
end

%% cross correlation function (CCF)
%plot all ccfs together with autos red and crosses blue 
stt=cell(18,1);
cht=cell(18,1);

figure %for stims
for nn=1:18 %loop through conditions
    if idx{nn}
        tst=[];
        for ii=idx{nn}
            tst=[tst;data_stim_dfof(:,:,ii)];
        end
        
        [c,lgs]=xcov(tst,25,'coeff');
        
        subplot(3,6,nn); hold on
        auto = c(26,:)>.99;
        shadedErrorBar(lgs,mean(c(:,auto),2),std(c(:,auto),[],2),'r')
        shadedErrorBar(lgs,mean(c(:,~auto),2),std(c(:,auto),[],2),'b')
        ylim([-.2,1])
        ylabel('correlation')
        xlabel('frames')
        
        stt{nn}=c(:,~auto);
        
    end
end

figure %for choice
for nn=1:18 %loop through conditions
    if idx{nn}
        tst=[];
        for ii=idx{nn}
            tst=[tst;data_choice_dfof(:,:,ii)];
        end
        
        [c,lgs]=xcov(tst,25,'coeff');
        
        subplot(3,6,nn); hold on
        auto = c(26,:)>.99;
        shadedErrorBar(lgs,mean(c(:,auto),2),std(c(:,auto),[],2),'r')
        shadedErrorBar(lgs,mean(c(:,~auto),2),std(c(:,auto),[],2),'b')
        ylim([-.2,1])
        ylabel('correlation')
        xlabel('frames')
        
        cht{nn}=c(:,~auto);
        
    end
end

%% calculate ccf significance
[S1_lr_h_stim, S1_lr_p_stim]=ttest(stt{13}',stt{16}','tail','right');
%3422 pairs, how correlated each pair, mean and std
%[1]comparing left and right high contrast success
%[2]comparing left high contrast success and ignore (no fails)
%[3]comparing left success with high contrast and middle, [4]then lower contrast
%[5]comparing right succcess and fail (no ignore) at high contrast
%[6][7]compare different contrasts for right success
first=[13 13 13 13 18 10 4];
second=[16 14 7 1 16 16 16];
p=-1*ones(7,1);
for nn=1:7
    avg=zeros(3422,1);
    for ii=1:3422
        c1=xcov(stt{first(nn)}(:,ii),stt{second(nn)}(:,ii),25,'coeff');
        avg(ii) = mean(c1(26:35));
    end
    u=mean(avg);
    o=std(avg);
    z=(u-1)/(o/sqrt(3422));
    p(nn)=normcdf(z);
    %figure; hist(avg)
end

%% test shuffle comparison
%everything still rejects null that they're correlated across conditions,
%takes a long time to run too
%{
for nn=1:1%7
    avg=zeros(3422,1);
    for jj=1:3422
    for ii=1:3422
        c1=xcov(stt{first(nn)}(:,jj),stt{second(nn)}(:,ii),25,'coeff');
        avg(ii) = mean(c1(26:35));
    end
    u=mean(avg);
    o=std(avg);
    z=(u-1)/(o/sqrt(3422));
    p(jj)=normcdf(z);
    %figure; hist(avg)
    end
end
%}

%% choice p vals
%[1]comparing left and right high contrast success
%[2]comparing left high contrast success and ignore (no fails)
%[3]comparing left success with high contrast and middle, [4]then lower contrast
%[5]comparing right succcess and fail (no ignore) at high contrast
%[6][7]compare different contrasts for right success
for nn=1:7
    avg=zeros(3422,1);
    for ii=1:3422
        c1=xcov(cht{first(nn)}(:,ii),cht{second(nn)}(:,ii),25,'coeff');
        avg(ii) = mean(c1(26:35));
    end
    u=mean(avg);
    o=std(avg);
    z=(u-1)/(o/sqrt(3422));
    p(nn)=normcdf(z);
    %figure; hist(avg)
end

%% raw traces of dfof for the cells, coded by responders to stim, choice, or nothing

h_stim;
h_choice;
nCells;

tt = (-20:29)*1000/frame_rate;%frames into msec for the 50 frames
%for stim
figure
for ii=1:6
    subplot(3,2,ii)
    hold on                                                          %)intersect(find(SIx,
    indS = intersect(find(tGratingContrast == cons(4-(round(ii/2)))), find(tLeftTrial  == mod(ii,2)));%high contrast and succesful and left, then right
    xlabel('msec')
    ylabel('dF/F')
    ylim([-.2,.4])
    for cc = 1:nCells%highLvR100 %
        if find(find(h_stim(1+mod(ii,2),:)) == cc)
            clr='b';
        else
            clr='r';
        end
        plot(tt,nanmean(data_stim_dfof(:,cc,indS),3),clr);%each cell's dF/F across time for all selected trials
        rt(cc,1)=max(nanmean(data_stim_dfof(:,cc,indS),3));
        %cell 38
    end
end

tcon=cell(3,1);
for ii=1:3
    tcon{ii}=num2str(cons(4-ii));
end

xx=[1,3,5];
for ii=1:3
    subplot(3,2,xx(ii))
    ylabel([tcon{ii},'           '],'FontWeight','bold','rot',0)
end

subplot(3,2,1)
title('Ipsilateral Stim')
subplot(3,2,2)
title('Contralateral Stim')
subplot(3,2,5)
xlabel('Blue responds to Ipsilateral Stim High Contrast, Red does not')
subplot(3,2,6)
xlabel('Blue responds to Contralateral Stim High Contrast, Red does not')

%for choice
figure
for ii=1:6
    subplot(3,2,ii)
    hold on
    indS = intersect(find(tGratingContrast == cons(4-(round(ii/2)))), intersect(find(SIx), find(tLeftTrial == mod(ii,2))));%high contrast and succesful and right, then left
    xlabel('msec')
    ylabel('dF/F')
    ylim([-.2,.8])
    for cc = 1:nCells
        if find(find(h_choice(1+mod(ii,2),:)) == cc)
            clr='b';
        else
            clr='r';
        end
        plot(tt,nanmean(data_choice_dfof(:,cc,indS),3),clr);%each cell's dF/F across time for all selected trials
        rt(cc,1)=max(nanmean(data_choice_dfof(:,cc,indS),3));
        %cell 38
    end
end

tcon=cell(3,1);
for ii=1:3
    tcon{ii}=num2str(cons(4-ii));
end

xx=[1,3,5];
for ii=1:3
    subplot(3,2,xx(ii))
    ylabel([tcon{ii},'           '],'FontWeight','bold','rot',0)
end

subplot(3,2,1)
title('Ipsilateral Stim')
subplot(3,2,2)
title('Contralateral Stim')
subplot(3,2,5)
xlabel('Blue responds to Choice for Ipsilateral Stim High Contrast, Red does not')
subplot(3,2,6)
xlabel('Blue responds to Choice for Contralateral Stim High Contrast, Red does not')

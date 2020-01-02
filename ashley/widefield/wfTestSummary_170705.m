load(fullfile(rc.structOutput, [mouse '_concatTC.mat']))
load(fullfile(rc.structOutput, [mouse '_roi_masks.mat']))
data_f = mean(data_stim_tc_all(:,1:10,:),2);
data_stim_dfof = dfof(data_stim_tc_all, data_f);
data_dec_dfof = dfof(data_dec_tc_all, data_f);

indR = intersect(intersect(find(tGratingContrast_all == 1), find(SIx_all)),find(tLeftTrial_all==0));
indL = intersect(intersect(find(tGratingContrast_all == 1), find(SIx_all)),find(tLeftTrial_all));
stimR = mean(data_stim_dfof(:,:,indR),3);
stimL = mean(data_stim_dfof(:,:,indL),3);
decR = mean(data_dec_dfof(:,:,indR),3);
decL = mean(data_dec_dfof(:,:,indL),3);
dateStr = [];
for i = 1:size(dateStr_all,1)
    dateStr = [dateStr ' ' dateStr_all(i,:)];
end
frameRateHz = 10;
tt = (1-10:1+29).*(1000/frameRateHz);
figure;
for i = 1:4
    subplot(2,2,i)
    plot(tt,stimR(i,:))
    hold on
    plot(tt,stimL(i,:))
    title(area_list(i,:))
end
suptitle([mouse dateStr '- Stim align: Blue- Contra; Red- Ipsi'])
print(fullfile(rc.structOutput, [mouse '_stim_IvC_byAreas.pdf']), '-dpdf', '-bestfit')

figure;
for i = 1:4
    subplot(2,2,i)
    plot(tt,decR(i,:))
    hold on
    plot(tt,decL(i,:))
    title(area_list(i,:))
end
suptitle([mouse dateStr '- Decision align: Blue- Contra; Red- Ipsi'])
print(fullfile(rc.structOutput, [mouse '_dec_IvC_byAreas.pdf']), '-dpdf', '-bestfit')

probs = unique(tProbLeft_all);
nprob = length(probs);
cons = unique(tGratingContrast_all);
ncon = length(cons);
figure;
[n n2] = subplotn(ncon);
for icon = 1:ncon
    indR = intersect(intersect(find(con == cons(icon)), find(SIx_all)),find(tLeftTrial_all==0));
    indL = intersect(intersect(find(con == cons(icon)), find(SIx_all)),find(tLeftTrial_all));
    subplot(n,n2,icon)
    plot(tt,mean(data_stim_dfof(1,:,indR),3));
    hold on
    plot(tt,mean(data_stim_dfof(1,:,indL),3));
    title(num2str(chop(cons(icon),2)*100))
    ylim([-0.01 0.03])
end
suptitle([mouse dateStr 'Blue = Contra; Red = Ipsi'])
print(fullfile(rc.structOutput, [mouse '_V1_IvC_byCon.pdf']), '-dpdf', '-bestfit')

figure;
[n n2] = subplotn(ncon);
for i = 1:4
    for icon = 1:ncon
        indR = intersect(intersect(find(con == cons(icon)), find(SIx_all)),find(tLeftTrial_all==0));
        subplot(2,2,i)
        plot(tt,mean(data_stim_dfof(i,:,indR),3));
        hold on
        title(area_list(i,:))
        ylim([-0.01 0.03])
    end
end
suptitle([mouse dateStr ' Blue = ' num2str(chop(cons(1),2)*100) '; Red = ' num2str(chop(cons(2),2)*100) '; Yellow = ' num2str(chop(cons(3),2)*100) '; Purple = ' num2str(chop(cons(4),2)*100)])
print(fullfile(rc.structOutput, [mouse '_cons_byArea.pdf']), '-dpdf', '-bestfit')

figure;
for icon = 1:ncon
    for i = 2:nprob-1
        indR = intersect(intersect(find(con == cons(icon)), find(SIx_all)),find(tLeftTrial_all==0));        
        indP = intersect(indR, find(tProbLeft_all == probs(i)));
        subplot(n,n2,icon)
        plot(tt,mean(data_stim_dfof(1,:,indP),3));
        hold on
    end
    ylim([-0.01 0.03])
    title(num2str(chop(cons(icon),2)*100))
end
suptitle([mouse dateStr ' Blue = ' num2str(1-probs(2)) ' ; Red = '  num2str(1-probs(3)) ' ; Yellow = '  num2str(1-probs(4))])
print(fullfile(rc.structOutput, [mouse '_probs_byCon.pdf']), '-dpdf', '-bestfit')

tGratingContrast = tGratingContrast_all;
tGratingContrast(find(tLeftTrial_all)) = tGratingContrast(find(tLeftTrial_all)).*-1;
cons = unique(tGratingContrast);
ncon = length(cons);
nSIx = zeros(ncon,nprob-2);
nFIx = zeros(ncon,nprob-2);
pC = zeros(ncon,nprob-2);
ciC = zeros(2,ncon,nprob-2);
for icon = 1:ncon
    indC = find(tGratingContrast == cons(icon));        
    for i = 2:nprob-1
        iprob = i-1;
        indP = intersect(indC, find(tProbLeft_all == probs(i)));
        nSIx(icon,iprob) = sum(SIx_all(:,indP),2);
        nFIx(icon,iprob) = sum(FIx_all(:,indP),2);
        [pC(icon,iprob) ciC(:,icon,iprob)] = binofit(nSIx(icon,iprob), nSIx(icon,iprob)+nFIx(icon,iprob));
    end
end

col_mat = strvcat('b', 'k', 'g');
figure;
for i = 1:nprob-2
    for icon = 1:ncon
        if cons(icon) < 0
            errorbar(cons(icon), 1-pC(icon,i), ciC(2,icon,i)-pC(icon,i), pC(icon,i)-ciC(1,icon,i),['o' col_mat(i,:)])
        else
            errorbar(cons(icon), pC(icon,i), pC(icon,i)-ciC(1,icon,i),ciC(2,icon,i)-pC(icon,i), ['o' col_mat(i,:)])
        end
        hold on
    end
end
ylabel('Prob right choice')
xlabel('Contrast (R-L)')
title([mouse dateStr 'Blue = ' num2str(probs(2)*100) '; Black = ' num2str(probs(3)*100) '; Green = ' num2str(probs(4)*100)])
print(fullfile(rc.structOutput, [mouse '_behavior_byProb.pdf']), '-dpdf', '-bestfit')

%% reshuffle data
analdir = 'G:\users\lindsey\analysisLG\active mice\DR7\110425\';
seqfile = '110425_DR7_run5_dec_Seqposition.mat';
load(fullfile(analdir,'analysis',seqfile));

nON = 144;
nOFF = 144;
nStim = 45;

for stim = 1:nStim
    stim_matrix_down(1+((stim-1)*nON):nON+((stim-1)*nON))= (nOFF+1+((stim-1)*(nON+nOFF)):(nON+nOFF)+((stim-1)*(nON+nOFF)));
    baseline_matrix_down(1+((stim-1)*(nOFF/3)):(nOFF/3)+((stim-1)*(nOFF/3))) =1+(2*nOFF/3) + ((stim-1)*(nON+nOFF)):nOFF +((stim-1)*(nON+nOFF));
end

for iS = 1:length(Seqposition)
    oris(iS)=Seqposition(iS).TFSFetc(3);
    cons(iS)=Seqposition(iS).TFSFetc(4);
end

uoris = unique(oris);
noris = length(uoris);

ic_oris = zeros(((nON+nOFF)),5, noris,nIC);
ic_oris_avg = zeros(((nON+nOFF)),1,noris,nIC);
for ic = 1:6
    figure;
    for ori = 1:noris
        subplot(1,noris,ori)
        for rep = 1:5;
            stim = Seqposition(ori).ind(rep);
            ic_oris(:,rep,ori,ic) = ica_sig(ic, 1+((stim-1)*((nON+nOFF))):((nON+nOFF))+((stim-1)*((nON+nOFF))));
            plot(ic_oris(:,rep,ori,ic), 'c');
            hold on; 
            axis([0 30 -.01 .1]);
            axis off
        end
        ic_oris_avg(:,:,ori,ic) = mean(ic_oris(:,:,ori,ic),2);
        plot(ic_oris_avg(:,:,ori,ic), 'k'); 
    end
end

dFoverF = zeros(5, noris, nIC);
for ic = 1:nIC;
    for ori = 1:noris
        for trial = 1:5;
            stim_avg = mean(ic_oris(1+(nON/10):(nOFF+nON)/10,trial, ori, ic),1);
            baseline_avg = mean(ic_oris(1+(2*nOFF/30):nOFF/10,trial, ori, ic),1);
            dFoverF(trial,ori,ic)= (stim_avg-abs(baseline_avg))./abs(baseline_avg);
        end
    end
end

dFoverF_avg = mean(dFoverF,1);
dFoverF_stdev = nanstd(dFoverF,[],1);
dFoverF_sem = dFoverF_stdev/sqrt(5); 

figure;
for ic = 1:nIC;
    subplot(5,5,ic)
    errorbar((uoris/pi)*180, dFoverF_avg(:,:,ic), dFoverF_sem(:,:,ic));
    xlim([0 360]);
    xlabel('Orientation (Deg)');
    ylabel('dF/F');
    hold on;
end

for ic = 5;


subplot(2,3,2);
errorbar((uoris/pi)*180, dFoverF_avg(:,:,ic), dFoverF_sem(:,:,ic));
xlim([0 360]);
xlabel('Orientation (Deg)');
ylabel('dF/F');
hold on;

ic = 5;
figure
for ori = 1:noris
    subplot(1,noris+1,ori)
    for rep = 1:5;
        stim = Seqposition(ori).ind(rep);
        ic_oris(:,rep,ori,ic) = ica_sig_down(ic, 1+((stim-1)*((nON+nOFF)/10)):((nON+nOFF)/10)+((stim-1)*((nON+nOFF)/10)));
        plot(ic_oris(:,rep,ori,ic), 'c');
        hold on; 
        axis([0 30 -.01 .1]);
        axis off
    end
    ic_oris_avg(:,:,ori,ic) = mean(ic_oris(:,:,ori,ic),2);
    plot(ic_oris_avg(:,:,ori,ic), 'k'); 
end
subplot(1,noris+1,noris +1);
imstretch(sm(:,:,ic),[.1 .99],1.5);
text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
title(sprintf('IC %i',ic))




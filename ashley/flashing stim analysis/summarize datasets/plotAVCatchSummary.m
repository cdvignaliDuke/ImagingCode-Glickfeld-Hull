function plotAVCatchSummary(datasetStr,cellsInd)
% takes mouse FSAV Ca structure and plots responses to anticipation phase
% of task
close all
av = behavParamsAV;
eval(['awFSAVdatasets' datasetStr])
titleStr = datasetStr;
if strcmp(titleStr, '')
    titleStr = 'V1';
else
    titleStr = titleStr(2:end);
end
rc = behavConstsAV;
if strcmp(rc.name,'ashley')
    dataGroup = ['awFSAVdatasets' datasetStr];
else
    dataGroup = [];
end
str = unique({expt.SubNum});
values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
mouse_str = ['i' strjoin(str,'_i')];
mouse_ind = find(intersect(cell2mat({av.mouse}),values));
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' datasetStr '.mat']));
pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
cycTime = mouse(1).expt(1).info(1).cyc_time;
cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name]; 
iav = 1; % for visual trials
% parameters for accessing data in mouse struct
catchAlign = 3;
targetAlign = 2;
cFA = 3;
cCR = 4;
hits = 1;
misses = 2;
hitmatch = 5;
missmatch = 6;
famatch = 7;
crmatch = 10;

% cycTime = mouse(1).expt(1).info(1).cyc_time;

fnout = fullfile('Z:\Analysis\temp figs\160331', [titleStr '_' mouse_str]); %% maybe lose mouse_str

nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
n = ceil(sqrt(nexp+1));
if (n^2)-n > nexp+1
    n2 = n-1;
else
    n2 = n;
end

%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% set(0,'DefaultaxesFontSize', 16)

tt = -pre_event_frames:post_event_frames-1;
ttMs = tt/(cycTime/cycTimeMs);
baseStimFrames = -(floor(pre_event_frames/cycTime)*cycTime):cycTime:0;

%% avg resp accross all direction for success,miss,catch-FA,catch-CR - target responsive cells
% 
% respTvsCFig = figure;
% suptitle(titleStr)
% i = 1;
% resp_tar_hits = [];
% resp_tar_misses = [];
% resp_catch_fa = [];
% resp_catch_cr = [];
% % resp_FA = [];
% % resp_CR = [];
% n_cFA = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
% n_cCR = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
% n_hits = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
% n_misses = zeros(1,size(mouse,2)+size(mouse(imouse).expt,2));
% 
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         figure(respTvsCFig);
%         subplot(n,n2,i)
%         respTrace = [];
%         if mouse(imouse).expt(iexp).info.isCatch
%             cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
%             cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
%             % invalid hits matched to invalid miss
%             respTrace(1) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),2),3), 'c');
%             n_cFA(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),3);
%             hold on
%             %invalid miss matched to invalid hits
%             respTrace(2) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),2),3), 'b');
%             n_cCR(i) = size(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),3);
%             hold on
%             %valid hits matched to valid miss
%             respTrace(3) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),2),3), 'k');
%             n_hits(i) = size(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),3);
%             hold on
%             %valid miss match to valid hits
%             respTrace(4) = plot(tt,mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),2),3), 'r');
%             n_misses(i) = size(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),3);
%             hold on
%             vline(baseStimFrames,':k')
%             hold on
%             vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%             hold on
%             vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
%             hold on
%             xlim([-10 20])
%             legend(respTrace,{[num2str(n_cFA(i)) ' invalid hit'];[num2str(n_cCR(i)) ' invalid miss'];[num2str(n_hits(i)) ' hit'];[num2str(n_misses(i)) ' miss']},'Location','northwest');
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells']})
%         
%             resp_tar_hits = cat(2, resp_tar_hits, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hitmatch).resp(:,cell_ind,:),3));
%             resp_tar_misses  = cat(2, resp_tar_misses, mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(missmatch).resp(:,cell_ind,:),3));
%             resp_catch_fa = cat(2, resp_catch_fa, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(famatch).resp(:,cell_ind,:),3));
%             resp_catch_cr = cat(2, resp_catch_cr, mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crmatch).resp(:,cell_ind,:),3));
%         end
%          i = i+1;       
%     end
% end
% 
% n_cFA_all = sum(n_cFA);
% n_cCR_all = sum(n_cCR);
% n_misses_all = sum(n_misses);
% n_hits_all = sum(n_hits);
% 
% respTrace = [];
% subplot(n,n2,i)
% respTrace(3) = plot(tt, nanmean(resp_tar_hits,2), 'k');
% hold on
% respTrace(4) = plot(tt, nanmean(resp_tar_misses,2), 'r');
% hold on
% respTrace(1) = plot(tt, nanmean(resp_catch_fa,2), 'c');
% hold on
% respTrace(2) = plot(tt, nanmean(resp_catch_cr,2), 'b');
% vline(baseStimFrames,':k')
% hold on
% vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
% hold on
% vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
% xlim([-10 20])
% legend(respTrace,{[num2str(n_cFA_all) ' cFA'];[num2str(n_cCR_all) ' cCR'];[num2str(n_hits_all) ' hits'];[num2str(n_misses_all) ' misses']},'Location','northwest');
% title(['All cells; n = ' num2str(size(resp_tar_hits,2))])
% figure(respTvsCFig);
% print([fnout 'catch_align_TCs' datasetStr '.pdf'], '-dpdf')

%% scatter transient response C-FA vs. T-suc; C-CR vs. T-miss;C-CR vs T-suc; T-miss vs T-suc
scatFAvsH = figure; %oucomes 1h 5fa
scatCRvsM = figure;%oucomes 4m 9cr
scatCRvsH = figure;%oucomes 2h 8cr
scatMvsH = figure; % align = 2; oucomes 5h 6m
scatCRvsFA = figure; %oucomes 7fa 10cr
i = 1;

resp_hvsfa_all = [];
resp_favsh_all = [];
resp_mvscr_all = [];
resp_crvsm_all = [];
resp_hvscr_all = [];
resp_crvsh_all = [];
resp_hvsm_all = [];
resp_mvsh_all = [];
resp_favscr_all = [];
resp_crvsfa_all = [];
n_hvsfa = zeros(1,size(expt,2));
n_mvscr = zeros(1,size(expt,2));
n_hvscr = zeros(1,size(expt,2));
n_hvsm = zeros(1,size(expt,2));
n_favscr = zeros(1,size(expt,2));
n_all = zeros(1,size(expt,2));

tc_hvsfa_all = [];
tc_favsh_all = [];
tc_mvscr_all = [];
tc_crvsm_all = [];
tc_hvscr_all = [];
tc_crvsh_all = [];
tc_hvsm_all = [];
tc_mvsh_all = [];
tc_favscr_all = [];
tc_crvsfa_all = [];

resp_val_all = [];
resp_inv_all = [];
tc_val_all = [];
tc_inv_all = [];

hits = 1;
misses = 2;
fas = 3;
crs = 4;

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            cell_ind = intersect(mouse(imouse).expt(iexp).cells(1).ind,cell_ind);
        % find number of each trial outcome type
            cdirs = mouse(imouse).expt(iexp).info.cDirs;
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp,'unif',false)),3,length(cdirs));
            nDirs_h = sz(3,:);
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp,'unif',false)),3,length(cdirs));
            nDirs_m = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_m = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp,'unif',false)),3,length(cdirs));
            nDirs_fa = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_fa = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp,'unif',false)),3,length(cdirs));
            nDirs_cr = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_cr = nT(ind3);
            end
            
        % fa vs hits matched
            ind = find(lt(nDirs_fa,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvsfa = nDirs_fa;
            elseif isempty(ind)
                nDirs_hvsfa = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_fa(ind) nDirs_h(ind2)]
                nDirs_hvsfa = nT(ind3);
            end
            if any(nDirs_hvsfa > 1)
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            tc_hvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
            tc_favsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
            for iboot = 1:100
                resp_rTh = [];
                resp_rTfa = [];
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_fa(:,:,:,iboot) = resp_rTfa;
            if iboot == 1
               tc_hvsfa = resp_rTh;
               tc_favsh = resp_rTfa;
            end
            end
            resp_h = mean(resp_h,4);
            resp_fa = mean(resp_fa,4);
            else
            resp_h = [];
            resp_fa = [];
            tc_hvsfa = [];
            tc_favsh = [];
            end
            
            tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(tc_hvsfa,3)); 
            tc_favsh_all = cat(2,tc_favsh_all,mean(tc_favsh,3));  
            
            if any(nDirs_hvsfa > 1)
            resp_hvsfa = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
            resp_favsh = squeeze(mean(mean(resp_fa(trans_win,:,:),3),1))-squeeze(mean(mean(resp_fa(pre_win,:,:),3),1));
            else
                resp_hvsfa = [];
                resp_favsh = [];
            end
            resp_val_all = cat(2,resp_val_all,resp_hvsfa);
            resp_inv_all = cat(2,resp_inv_all,resp_favsh);
            tc_val_all = cat(2,tc_val_all,mean(tc_hvsfa,3));
            tc_inv_all = cat(2,tc_inv_all,mean(tc_favsh,3));
            
%             resp_hvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_favsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
            resp_hvsfa_all = cat(2,resp_hvsfa_all,resp_hvsfa);
            resp_favsh_all = cat(2,resp_favsh_all,resp_favsh);
            n_hvsfa(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1));
            figure(scatFAvsH);
            subplot(n,n2,i)
            scatter(resp_hvsfa,resp_favsh,50,'k.')
            hold on
            errorbarxy(nanmean(resp_hvsfa),nanmean(resp_favsh),nanstd(resp_hvsfa)/length(resp_hvsfa),nanstd(resp_favsh)/length(resp_favsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('invalid hits (fa)')
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_hvsfa(i)) ' valid hits;' num2str(n_favsh(i)) ' invalid hits']})
                    

        % cr vs miss matched
%             misses = 4;
%             crs = 9;

            ind = find(lt(nDirs_cr,nDirs_m));
            if length(ind) == length(cdirs)
                nDirs_mvscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_mvscr = nDirs_m;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_m(ind2)];
                nDirs_mvscr = nT(ind3);
            end
            if any(nDirs_mvscr > 1)
            resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
            tc_mvscr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)));
            tc_crvsm = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)));
            for iboot = 1:100
                resp_rTm = [];
                resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                end
            end
            resp_m(:,:,:,iboot) = resp_rTm;
            resp_cr(:,:,:,iboot) = resp_rTcr;            
            if iboot == 1
               tc_mvscr = resp_rTm;
               tc_crvsm = resp_rTcr;
            end
            end
            resp_m = mean(resp_m,4);
            resp_cr = mean(resp_cr,4);
            else
                resp_m = [];
                resp_cr = [];
                tc_mvscr = [];
                tc_crvsm = [];
            end
            
            tc_mvscr_all = cat(2,tc_mvscr_all,mean(tc_mvscr,3)); 
            tc_crvsm_all = cat(2,tc_crvsm_all,mean(tc_crvsm,3));  
            
            if any(nDirs_mvscr > 1)
                    resp_mvscr = squeeze(mean(mean(resp_m(trans_win,:,:),3),1))-squeeze(mean(mean(resp_m(pre_win,:,:),3),1));
                resp_crvsm = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_mvscr = [];
                resp_crvsm = [];
            end
            resp_val_all = cat(2,resp_val_all,resp_mvscr);
            resp_inv_all = cat(2,resp_inv_all,resp_crvsm);
            tc_val_all = cat(2,tc_val_all,mean(tc_mvscr,3));
            tc_inv_all = cat(2,tc_inv_all,mean(tc_crvsm,3));
            n_all(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1))+sum(nDirs_mvscr(nDirs_mvscr > 1));
            
%             resp_mvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_mvscr_all = cat(2,resp_mvscr_all,resp_mvscr);
            resp_crvsm_all = cat(2,resp_crvsm_all,resp_crvsm);
            n_mvscr(i) = sum(nDirs_mvscr(nDirs_mvscr > 1));
            figure(scatCRvsM);
            subplot(n,n2,i)
            scatter(resp_mvscr,resp_crvsm,50,'k.')
            hold on
            errorbarxy(nanmean(resp_mvscr),nanmean(resp_crvsm),nanstd(resp_mvscr)/length(resp_mvscr),nanstd(resp_crvsm)/length(resp_crvsm),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid miss')
            ylabel('invalid miss')
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsm(i)) ' valid miss;' num2str(n_crvsm(i)) ' invalid miss']})
                  
        % cr vs hits matched
%             hits = 2;
%             crs = 8;

            ind = find(lt(nDirs_cr,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_hvscr = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_h(ind2)];
                nDirs_hvscr = nT(ind3);
            end
            if any(nDirs_hvscr > 1)
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)),100);
            tc_hvscr = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
            tc_crvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
            for iboot = 1:100
            resp_rTh = [];
            resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_hvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_cr(:,:,:,iboot) = resp_rTcr;
            if iboot == 1
               tc_hvscr = resp_rTh;
               tc_crvsh = resp_rTcr;
            end
            end
            resp_h = mean(resp_h,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_h = [];
            resp_cr = [];
            tc_hvscr = [];
            tc_crvsh = [];
            end   
            
            tc_hvscr_all = cat(2,tc_hvscr_all,mean(tc_hvscr,3)); 
            tc_crvsh_all = cat(2,tc_crvsh_all,mean(tc_crvsh,3));  
            
            if any(nDirs_hvscr > 1)
                resp_hvscr = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
                resp_crvsh = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_hvscr = [];
                resp_crvsh = [];
            end

%             resp_hvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_hvscr_all = cat(2,resp_hvscr_all,resp_hvscr);
            resp_crvsh_all = cat(2,resp_crvsh_all,resp_crvsh);
            n_hvscr(i) = sum(nDirs_hvscr(nDirs_hvscr > 1));
            figure(scatCRvsH);
            subplot(n,n2,i)
            scatter(resp_hvscr,resp_crvsh,50,'k.')
            hold on
            errorbarxy(nanmean(resp_hvscr),nanmean(resp_crvsh),nanstd(resp_hvscr)/length(resp_hvscr),nanstd(resp_crvsh)/length(resp_crvsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('invalid miss')
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsh(i)) ' valid hits;' num2str(n_crvsh(i)) ' invalid miss']})
        
        % hits vs miss matched
%             hits = 5;
%             misses = 6;

            ind = find(lt(nDirs_m,nDirs_h));
            if length(ind) == length(cdirs)
                nDirs_hvsm = nDirs_m;
            elseif isempty(ind)
                nDirs_hvsm = nDirs_h;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_m(ind) nDirs_h(ind2)];
                nDirs_hvsm = nT(ind3);
            end
            if any(nDirs_hvsm > 1)
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)),100);
            resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)),100);
            tc_hvsm = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
            tc_mvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
            for iboot = 1:100
                resp_rTh = [];
                resp_rTm = [];
            for idir = 1:length(cdirs)
                if nDirs_hvsm(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_m(:,:,:,iboot) = resp_rTm;
            if iboot == 1
               tc_hvsm = resp_rTh;
               tc_mvsh = resp_rTm;
            end
            end
            resp_h = mean(resp_h,4);
            resp_m = mean(resp_m,4);
            else
            resp_h = [];
            resp_m = [];
            tc_hvsm = [];
            tc_mvsh = [];
            end
            
            
            tc_hvsm_all = cat(2,tc_hvsm_all,mean(tc_hvsm,3)); 
            tc_mvsh_all = cat(2,tc_mvsh_all,mean(tc_mvsh,3));  
            
            if any(nDirs_hvsm > 1)
                resp_hvsm = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
                resp_mvsh = squeeze(mean(mean(resp_m(trans_win,:,:),3),1))-squeeze(mean(mean(resp_m(pre_win,:,:),3),1));
            else
                resp_hvsm = [];
                resp_mvsh = [];
            end

%             resp_hvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_mvsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(targetAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
            resp_hvsm_all = cat(2,resp_hvsm_all,resp_hvsm);
            resp_mvsh_all = cat(2,resp_mvsh_all,resp_mvsh);
            n_hvsm(i) = sum(nDirs_hvsm(nDirs_hvsm > 1));
            figure(scatMvsH);
            subplot(n,n2,i)
            scatter(resp_hvsm,resp_mvsh,50,'k.')
            hold on
            errorbarxy(nanmean(resp_hvsm),nanmean(resp_mvsh),nanstd(resp_hvsm)/length(resp_hvsm),nanstd(resp_mvsh)/length(resp_mvsh),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('valid hits')
            ylabel('valid miss')
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_mvsh(i)) ' valid hits;' num2str(n_mvsh(i)) ' valid miss']})
            
        % fa vs cr matched
%             fas = 7;
%             crs = 10;

            ind = find(lt(nDirs_cr,nDirs_fa));
            if length(ind) == length(cdirs)
                nDirs_favscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_favscr = nDirs_fa;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_fa(ind2)];
                nDirs_favscr = nT(ind3);
            end
            if any(nDirs_favscr > 1)
            resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);           
            tc_favscr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
            tc_crvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
            for iboot = 1:100
                resp_rTfa = [];
                resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_favscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                end
            end
            resp_fa(:,:,:,iboot) = resp_rTfa;
            resp_cr(:,:,:,iboot) = resp_rTcr;
            if iboot == 1
               tc_favscr = resp_rTfa;
               tc_crvsfa = resp_rTcr;
            end
            end
            resp_fa = mean(resp_fa,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_fa = [];
            resp_cr = [];
            tc_favscr = [];
            tc_crvsfa = [];
            end
            
            tc_favscr_all = cat(2,tc_favscr_all,mean(tc_favscr,3)); 
            tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(tc_crvsfa,3));  
            
            if any(nDirs_favscr > 1)
                resp_favscr = squeeze(mean(mean(resp_fa(trans_win,:,:),3),1))-squeeze(mean(mean(resp_fa(pre_win,:,:),3),1));
                resp_crvsfa = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_favscr = [];
                resp_crvsfa = [];
            end

%             resp_favscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_favscr_all = cat(2,resp_favscr_all,resp_favscr);
            resp_crvsfa_all = cat(2,resp_crvsfa_all,resp_crvsfa);
            n_favscr(i) = sum(nDirs_favscr(nDirs_favscr > 1));
            figure(scatCRvsFA);
            subplot(n,n2,i)
            scatter(resp_favscr,resp_crvsfa,50,'k.')
            hold on
            errorbarxy(nanmean(resp_favscr),nanmean(resp_crvsfa),nanstd(resp_favscr)/length(resp_favscr),nanstd(resp_crvsfa)/length(resp_crvsfa),{'ro','r','r'});
            xlim([-0.05 0.1])
            ylim([-0.05 0.1])
            hold on
            plot([-10:0.1:20],[-10:0.1:20],'k--')
            axis square
            xlabel('invalid hits')
            ylabel('invalid miss')
%             title({mouse(imouse).expt(iexp).date, [' n = ' num2str(length(cell_ind)) ' cells'],[num2str(n_crvsfa(i)) ' invalid hits;' num2str(n_crvsfa(i)) ' invalid miss']})
            
%combo vaild h+m, invalid h+m (fa+cr)
        resp_v = resp_h;
        resp_inv = resp_fa;
        tc_v = tc_hvsfa;
        tc_inv = tc_favsh;
        ind = find(lt(nDirs_cr,nDirs_fa));
            if length(ind) == length(cdirs)
                nDirs_favscr = nDirs_cr;
            elseif isempty(ind)
                nDirs_favscr = nDirs_fa;
            else
                ind2 = setdiff(1:length(cdirs),ind);
                [j ind3] = sort([ind ind2]);
                nT = [nDirs_cr(ind) nDirs_fa(ind2)];
                nDirs_favscr = nT(ind3);
            end
            if any(nDirs_favscr > 1)
            resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);           
            tc_favscr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
            tc_crvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
            for iboot = 1:100
                resp_rTfa = [];
                resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_favscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                end
            end
            resp_fa(:,:,:,iboot) = resp_rTfa;
            resp_cr(:,:,:,iboot) = resp_rTcr;
            if iboot == 1
               tc_favscr = resp_rTfa;
               tc_crvsfa = resp_rTcr;
            end
            end
            resp_fa = mean(resp_fa,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_fa = [];
            resp_cr = [];
            tc_favscr = [];
            tc_crvsfa = [];
            end
            
            tc_favscr_all = cat(2,tc_favscr_all,mean(tc_favscr,3)); 
            tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(tc_crvsfa,3));  
            
            if any(nDirs_favscr > 1)
                resp_favscr = squeeze(mean(mean(resp_fa(trans_win,:,:),3),1))-squeeze(mean(mean(resp_fa(pre_win,:,:),3),1));
                resp_crvsfa = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
            else
                resp_favscr = [];
                resp_crvsfa = [];
            end

%             resp_favscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_favscr_all = cat(2,resp_favscr_all,resp_favscr);
            resp_crvsfa_all = cat(2,resp_crvsfa_all,resp_crvsfa);
            n_favscr(i) = sum(nDirs_favscr(nDirs_favscr > 1));
            
        end
        i = i+1;
    end
end
        
figure(scatFAvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvsfa_all,resp_favsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvsfa_all),mean(resp_favsh_all),std(resp_hvsfa_all)/length(resp_hvsfa_all),std(resp_favsh_all)/length(resp_favsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('invalid hits')
% title({[' All cells = ' num2str(length(resp_hvsfa_all)) ' cells'],[num2str(sum(n_hvsfa)) ' valid hits;' num2str(sum(n_favsh)) ' invalid hits']})

cdfCatchResp = figure;
subplot(3,2,1)
hHvsFA = cdfplot(resp_hvsfa_all);
hHvsFA.Color = 'k';
hold on
hFAvsH = cdfplot(resp_favsh_all);
hFAvsH.Color = 'c';
[h p] = kstest2(resp_hvsfa_all,resp_favsh_all);
xlim([-0.05 0.05])
title(['H-v vs H-inv ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

tcCatchResp = figure;
suptitle([mouse(1).expt(1).cells(cellsInd).name ' pref'])
subplot(3,2,1)
plot(ttMs,mean(tc_hvsfa_all,2)-mean(mean(tc_hvsfa_all(pre_event_frames,:),2)),'k')
hold on
plot(ttMs,mean(tc_favsh_all,2)-mean(mean(tc_favsh_all(pre_event_frames,:),2)),'c')
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({[num2str(sum(n_hvsfa)) ' H-v, H-inv; ' num2str(size(resp_hvsfa_all,2)) ' cells']})

figure(scatCRvsM);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_mvscr_all,resp_crvsm_all,50,'k.')
hold on
errorbarxy(mean(resp_mvscr_all),mean(resp_crvsm_all),std(resp_mvscr_all)/length(resp_mvscr_all),std(resp_crvsm_all)/length(resp_crvsm_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid miss')
ylabel('invalid miss')
% title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsm)) ' valid miss;' num2str(sum(n_crvsm)) ' invalid miss']})

figure(cdfCatchResp)
subplot(3,2,2)
hHvsFA = cdfplot(resp_mvscr_all);
hHvsFA.Color = 'r';
hold on
hFAvsH = cdfplot(resp_crvsm_all);
hFAvsH.Color = 'b';
xlim([-0.05 0.05])
[h p] = kstest2(resp_mvscr_all,resp_crvsm_all);
title(['M-v vs M-inv ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

figure(tcCatchResp);
subplot(3,2,2)
plot(ttMs,mean(tc_mvscr_all,2)-mean(mean(tc_mvscr_all(pre_event_frames,:),2)),'r')
hold on
plot(ttMs,mean(tc_crvsm_all,2)-mean(mean(tc_crvsm_all(pre_event_frames,:),2)),'b')
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({[num2str(sum(n_mvscr)) ' M-v, M-inv; ' num2str(size(resp_mvscr_all,2)) ' cells']})

figure(scatCRvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvscr_all,resp_crvsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvscr_all),mean(resp_crvsh_all),std(resp_hvscr_all)/length(resp_hvscr_all),std(resp_crvsh_all)/length(resp_crvsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('invalid miss')
% title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsh)) 'valid hits;' num2str(sum(n_crvsh)) ' invalid miss']})

figure(cdfCatchResp)
subplot(3,2,3)
hHvsFA = cdfplot(resp_hvscr_all);
hHvsFA.Color = 'k';
hold on
hFAvsH = cdfplot(resp_crvsh_all);
hFAvsH.Color = 'b';
[h p] = kstest2(resp_hvscr_all,resp_crvsh_all);
xlim([-0.05 0.05])
title(['H-v vs M-inv ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

figure(tcCatchResp);
subplot(3,2,3)
plot(ttMs,mean(tc_hvscr_all,2)-mean(mean(tc_hvscr_all(pre_event_frames,:),2)),'k')
hold on
plot(ttMs,mean(tc_crvsh_all,2)-mean(mean(tc_crvsh_all(pre_event_frames,:),2)),'b')
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({[num2str(sum(n_hvscr)) ' H-v, M-inv; ' num2str(size(resp_hvscr_all,2)) ' cells']})

figure(scatMvsH);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_hvsm_all,resp_mvsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvsm_all),mean(resp_mvsh_all),std(resp_hvsm_all)/length(resp_hvsm_all),std(resp_mvsh_all)/length(resp_mvsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('valid miss')
% title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_mvsh)) ' valid hits;' num2str(sum(n_mvsh)) ' valid miss']})

figure(cdfCatchResp)
subplot(3,2,4)
hHvsFA = cdfplot(resp_hvsm_all);
hHvsFA.Color = 'k';
hold on
hFAvsH = cdfplot(resp_mvsh_all);
hFAvsH.Color = 'r';
[h p] = kstest2(resp_hvsm_all,resp_mvsh_all);
xlim([-0.05 0.05])
title(['H-v vs M-v ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

figure(tcCatchResp);
subplot(3,2,4)
plot(ttMs,mean(tc_hvsm_all,2)-mean(mean(tc_hvsm_all(pre_event_frames,:),2)),'k')
hold on
plot(ttMs,mean(tc_mvsh_all,2)-mean(mean(tc_mvsh_all(pre_event_frames,:),2)),'r')
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({[num2str(sum(n_hvsm)) ' H-v, M-v; ' num2str(size(resp_hvsm_all,2)) ' cells']})

figure(scatCRvsFA);
suptitle(titleStr)
subplot(n,n2,i)
scatter(resp_favscr_all,resp_crvsfa_all,50,'k.')
hold on
errorbarxy(mean(resp_favscr_all),mean(resp_crvsfa_all),std(resp_favscr_all)/length(resp_favscr_all),std(resp_crvsfa_all)/length(resp_crvsfa_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('invalid hit')
ylabel('invalid miss')
% title({[' All cells = ' num2str(length(cell_ind)) ' cells'],[num2str(sum(n_crvsfa)) ' invalid hit;' num2str(sum(n_crvsfa)) ' invalid miss']})

figure(cdfCatchResp)
subplot(3,2,5)
hHvsFA = cdfplot(resp_favscr_all);
hHvsFA.Color = 'c';
hold on
hFAvsH = cdfplot(resp_crvsfa_all);
hFAvsH.Color = 'b';
[h p] = kstest2(resp_favscr_all,resp_crvsfa_all);
xlim([-0.05 0.05])
title(['H-inv vs M-inv ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

figure(tcCatchResp);
subplot(3,2,5)
plot(ttMs,mean(tc_favscr_all,2)-mean(mean(tc_favscr_all(pre_event_frames,:),2)),'c')
hold on
plot(ttMs,mean(tc_crvsfa_all,2)-mean(mean(tc_crvsfa_all(pre_event_frames,:),2)),'b')
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title({[num2str(sum(n_favscr)) ' H-inv, M-inv; ' num2str(size(resp_favscr_all,2)) ' cells']})

% combo valid h+invalid h, valid m+invalid m
figure(cdfCatchResp)
subplot(3,2,6)
hVvsInv = cdfplot(resp_val_all);
hVvsInv.Color = 'k';
hold on
hInvvsV = cdfplot(resp_inv_all);
hInvvsV.Color = 'c';
[h p] = kstest2(resp_val_all,resp_inv_all);
xlim([-0.05 0.05])
title(['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)])

figure(tcCatchResp)
subplot(3,2,6)
plot(ttMs,mean(tc_val_all,2)-mean(mean(tc_val_all(pre_win,:),2)),'k');
hold on
plot(ttMs,mean(tc_inv_all,2)-mean(mean(tc_inv_all(pre_win,:),2)),'c');
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
hold on
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title([num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']);




figure(scatFAvsH);
print([fnout 'catch_align_FAvsH' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsM);
print([fnout 'catch_align_CRvsM' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsH);
print([fnout 'catch_align_CRvsH' datasetStr '.pdf'], '-dpdf')
figure(scatMvsH);
print([fnout 'catch_align_MvsH' datasetStr '.pdf'], '-dpdf')
figure(scatCRvsFA);
print([fnout 'catch_align_CRvsFA' datasetStr '.pdf'], '-dpdf')
figure(cdfCatchResp)
print([fnout 'catch_align_cdf' datasetStr '.pdf'], '-dpdf')
figure(tcCatchResp)
print([fnout 'catch_align_TC' datasetStr '.pdf'], '-dpdf')

% %% plot average 
% %  plot avg target trace for random subset of cells, all directions
% Dirs = mouse(end).info.allDirs; 
% colorsT = brewermap(length(Dirs)+15,'YlGn');
% colorindT = [3:2:length(Dirs)+12];
% colorsT = colorsT(colorindT(1:length(Dirs)),:);
% 
% 
% for imouse = 1:size(mouse,2)
%     for iexp = 1:size(mouse(imouse).expt,2)
%         figure;
%         suptitle({mouse(imouse).expt(iexp).date, mouse(imouse).expt(iexp).mouse_name})
%         cellTarResp = [];
%         dirLegend = [];
%         cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
%         if length(cell_ind) > 9
%             cell_ind = sort(randsample(cell_ind,9));
%         end
%         for icell = 1:length(cell_ind)
%             subplot(3,3,icell)
%             exptDirs = mouse(imouse).expt(iexp).visTargets;
%             for idir = 1:length(exptDirs)  
%                 dirColInd = find(Dirs == exptDirs(idir));
%                 cellTarResp(idir) = plot(tt,mean(mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).stimResp{idir}(:,cell_ind(icell),:),3),'color',colorsT(dirColInd,:));
%                 dirLegend{idir} = num2str(exptDirs(idir));
%                 hold on                
%             end
%             legend(cellTarResp,dirLegend,'Location','northwest')
%             vline(baseStimFrames,':k')
%             hold on
%             vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames], '--r')
%             hold on
%             vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames], '--k')
%             hold on
%             xlim([-10 10])
%             ylim([-0.1 0.1])
%             title(['cell# ' num2str(cell_ind(icell))])
%         end
%     end
% end

end


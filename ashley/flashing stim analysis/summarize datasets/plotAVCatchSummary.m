function plotAVCatchSummary(datasetStr,cellsInd)
% takes mouse FSAV Ca structure and plots responses to target phase
% of task, comparing valid and invalid trials
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

fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' mouse_str]); %% maybe lose mouse_str

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
resp_all_all = [];
tc_val_all = [];
tc_inv_all = [];
tc_all = [];

%find all catch dirs used
cds = [];
trLs = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)   
        cds = [cds mouse(imouse).expt(iexp).info.cDirs];
        trLs = [trLs mouse(imouse).expt(iexp).align(targetAlign).av(1).outcome(1).tcyc*mouse(imouse).expt(iexp).info.cyc_time_ms];
    end
end
cds = unique(cds);
trLs = unique(trLs);

resp_h_dir = cell(length(cds),1);
resp_fa_dir = cell(length(cds),1);
resp_m_dir = cell(length(cds),1);
resp_cr_dir = cell(length(cds),1);
val_all_dir = cell(length(cds),1);
inv_all_dir = cell(length(cds),1);
n_h_dir = zeros(length(cds),size(expt,2));
n_m_dir = zeros(length(cds),size(expt,2));

hits = 1;
misses = 2;
fas = 3;
crs = 4;

oriTuningResp_avg = [];
oriTuningResp_tc = [];

for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            val_temp = [];
            inv_temp = [];
            all_temp = [];
            cell_ind = mouse(imouse).expt(iexp).cells(cellsInd).ind;
            cell_ind = intersect(mouse(imouse).expt(iexp).cells(13).ind,cell_ind);
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
            
            if sum(nDirs_fa > 1) > 0 & sum(nDirs_h > 1) > 0  
                mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
                dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
                load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
                oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
                oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
            elseif sum(nDirs_m > 1) > 0 & sum(nDirs_cr > 1) > 0 
                mName = ['AW' mouse(imouse).expt(iexp).mouse_name(2:end)];
                dirtuning = expt(intersect( find(strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)) ,find(strcmp({expt.date},mouse(imouse).expt(iexp).date)) ) ).dirtuning;
                load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
                oriTuningResp_avg = cat(2,oriTuningResp_avg,dFoverF_OriResp_avg_rect(:,cell_ind));
                oriTuningResp_tc = cat(3,oriTuningResp_tc,dFoverF_OriResp_TC(:,:,cell_ind));
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
            resp_h_dir_temp = cell(length(cdirs),1);
            resp_fa_dir_temp = cell(length(cdirs),1);
            if any(nDirs_hvsfa > 1)
%             resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa),100);
%             resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa),100);           
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            for iboot = 1:100
                resp_rTh = [];
                resp_rTfa = [];
                rH_dir_temp = cell(length(cdirs),1);
                rFA_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rH_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_h_dir_temp{idir} = cat(4,resp_h_dir_temp{idir},rH_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rFA_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_fa_dir_temp{idir} = cat(4,resp_fa_dir_temp{idir},rFA_dir_temp{idir});
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_fa(:,:,:,iboot) = resp_rTfa;
%             if iboot == 1
%                tc_hvsfa = resp_rTh;
%                tc_favsh = resp_rTfa;
%             end
            end
            resp_h = nanmean(resp_h,4);
            resp_fa = nanmean(resp_fa,4);
            
%             tc_hvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
%             tc_favsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa >1)));
            
            
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                ind = find(cds == cdirs(idir));
                resp_h_dir{ind,:} = cat(2,resp_h_dir{ind,:},squeeze(nanmean(nanmean(resp_h_dir_temp{idir},4),3)));
                resp_fa_dir{ind,:} = cat(2,resp_fa_dir{ind,:},squeeze(nanmean(nanmean(resp_fa_dir_temp{idir},4),3)));
                end
            end
            else
            resp_h = [];
            resp_fa = [];
%             tc_hvsfa = [];
%             tc_favsh = [];
            end
            
            
            tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(resp_h,3));%tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(tc_hvsfa,3)); 
            tc_favsh_all = cat(2,tc_favsh_all,mean(resp_fa,3));%tc_favsh_all = cat(2,tc_favsh_all,mean(tc_favsh,3));  
            
            if any(nDirs_hvsfa > 1)
            resp_hvsfa = squeeze(mean(mean(resp_h(trans_win,:,:),3),1))-squeeze(mean(mean(resp_h(pre_win,:,:),3),1));
            resp_favsh = squeeze(mean(mean(resp_fa(trans_win,:,:),3),1))-squeeze(mean(mean(resp_fa(pre_win,:,:),3),1));
            val_temp = resp_h;
            inv_temp = resp_fa;
            all_temp = cat(3,resp_h,resp_fa);
            else
                resp_hvsfa = [];
                resp_favsh = [];
                val_temp = [];
                inv_temp = [];
                all_temp = [];
            end

%             resp_val_all = cat(2,resp_val_all,resp_hvsfa);
%             resp_inv_all = cat(2,resp_inv_all,resp_favsh);
%             tc_val_all = cat(2,tc_val_all,mean(resp_h,3));%tc_val_all = cat(2,tc_val_all,mean(tc_hvsfa,3));
%             tc_inv_all = cat(2,tc_inv_all,mean(resp_fa,3));%tc_inv_all = cat(2,tc_inv_all,mean(tc_favsh,3));
%             
%             resp_hvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).resp(pre_win,cell_ind,:),3),1));
%             resp_favsh = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
            resp_hvsfa_all = cat(2,resp_hvsfa_all,resp_hvsfa);
            resp_favsh_all = cat(2,resp_favsh_all,resp_favsh);
            n_hvsfa(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1));
            for idir = find(nDirs_hvsfa > 1)
                ind = find(cds == cdirs(idir));
                n_h_dir(ind,i) = nDirs_hvsfa(idir);
            end
            figure(scatFAvsH);
            subplot(n,n2,i)
            scatter(resp_hvsfa,resp_favsh,50,'k.')
            hold on
            cc(nanmean(resp_hvsfa),nanmean(resp_favsh),nanstd(resp_hvsfa)/length(resp_hvsfa),nanstd(resp_favsh)/length(resp_favsh),{'ro','r','r'});
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
            resp_m_dir_temp = cell(length(cdirs),1);
            resp_cr_dir_temp = cell(length(cdirs),1);
            if any(nDirs_mvscr > 1)
            resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
            resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr(nDirs_mvscr > 1)),100);
%             resp_m = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr),100);
%             resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_mvscr),100);
            for iboot = 1:100
                resp_rTm = [];
                resp_rTcr = [];
                rM_dir_temp = cell(length(cdirs),1);
                rCR_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rM_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_m_dir_temp{idir} = cat(4,resp_m_dir_temp{idir},rM_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rCR_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_cr_dir_temp{idir} = cat(4,resp_cr_dir_temp{idir},rCR_dir_temp{idir});
                end
            end
            end
            resp_m(:,:,:,iboot) = resp_rTm;
            resp_cr(:,:,:,iboot) = resp_rTcr;            
%             if iboot == 1
%                tc_mvscr = resp_rTm;
%                tc_crvsm = resp_rTcr;
%             end
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                ind = find(cds == cdirs(idir));
                resp_m_dir{ind,:} = cat(2,resp_m_dir{ind,:},squeeze(nanmean(nanmean(resp_m_dir_temp{idir},4),3)));
                resp_cr_dir{ind,:} = cat(2,resp_cr_dir{ind,:},squeeze(nanmean(nanmean(resp_cr_dir_temp{idir},4),3)));
                end
            end
            
            resp_m = nanmean(resp_m,4);
            resp_cr = nanmean(resp_cr,4);
            else
                resp_m = [];
                resp_cr = [];
%                 tc_mvscr = [];
%                 tc_crvsm = [];
            end
            
            tc_mvscr_all = cat(2,tc_mvscr_all,mean(resp_m,3)); 
            tc_crvsm_all = cat(2,tc_crvsm_all,mean(resp_cr,3)); 
            
            if any(nDirs_mvscr > 1)
                    resp_mvscr = squeeze(mean(mean(resp_m(trans_win,:,:),3),1))-squeeze(mean(mean(resp_m(pre_win,:,:),3),1));
                resp_crvsm = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
                val_temp = cat(3,val_temp,resp_m);
                inv_temp = cat(3,inv_temp,resp_cr);
                all_temp = cat(3,all_temp,resp_m,resp_cr);
            else
                resp_mvscr = [];
                resp_crvsm = [];
                val_temp = val_temp;
                inv_temp = inv_temp;
                all_temp = all_temp;
            end
            
            resp_val_all = cat(2,resp_val_all,squeeze(mean(mean(val_temp(trans_win,:,:),3),1))-squeeze(mean(mean(val_temp(pre_win,:,:),3),1)));%cat(2,resp_val_all,resp_mvscr);
            resp_inv_all = cat(2,resp_inv_all,squeeze(mean(mean(inv_temp(trans_win,:,:),3),1))-squeeze(mean(mean(inv_temp(pre_win,:,:),3),1)));%cat(2,resp_inv_all,resp_crvsm);
            resp_all_all = cat(2,resp_all_all,squeeze(mean(mean(all_temp(trans_win,:,:),3),1))-squeeze(mean(mean(all_temp(pre_win,:,:),3),1)));
            tc_val_all = cat(2,tc_val_all,mean(val_temp,3));%mean(resp_m,3));
            tc_inv_all = cat(2,tc_inv_all,mean(inv_temp,3));%mean(resp_cr,3));
            tc_all = cat(2,tc_all,mean(all_temp,3));
            n_all(i) = sum(nDirs_hvsfa(nDirs_hvsfa > 1))+sum(nDirs_mvscr(nDirs_mvscr > 1));
            
%             resp_mvscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).resp(pre_win,cell_ind,:),3),1));
%             resp_crvsm = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
            resp_mvscr_all = cat(2,resp_mvscr_all,resp_mvscr);
            resp_crvsm_all = cat(2,resp_crvsm_all,resp_crvsm);
            n_mvscr(i) = sum(nDirs_mvscr(nDirs_mvscr > 1));
            for idir = find(nDirs_mvscr > 1)
                ind = find(cds == cdirs(idir));
                n_m_dir(ind,i) = nDirs_mvscr(idir);
            end
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
      
        % all valid and invalid matched
        tH = cellfun(@(resp_h_dir_temp) nanmean(resp_h_dir_temp,4),resp_h_dir_temp,'unif',false);
        tM = cellfun(@(resp_m_dir_temp) nanmean(resp_m_dir_temp,4),resp_m_dir_temp,'unif',false);
        tFA = cellfun(@(resp_fa_dir_temp) nanmean(resp_fa_dir_temp,4),resp_fa_dir_temp,'unif',false);
        tCR = cellfun(@(resp_cr_dir_temp) nanmean(resp_cr_dir_temp,4),resp_cr_dir_temp,'unif',false);
        
        for idir = 1:length(cdirs)
            tVal = mean(cat(3,tH{idir},tM{idir}),3);
            tInv = mean(cat(3,tFA{idir},tCR{idir}),3);
            ind = find(cds == cdirs(idir));
            val_all_dir{ind} = cat(2,val_all_dir{ind},tVal);
            inv_all_dir{ind} = cat(2,inv_all_dir{ind},tInv);
        end

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
%             tc_hvscr = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
%             tc_crvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvscr(nDirs_hvscr > 1)));
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
%             if iboot == 1
%                tc_hvscr = resp_rTh;
%                tc_crvsh = resp_rTcr;
%             end
            end
            resp_h = mean(resp_h,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_h = [];
            resp_cr = [];
%             tc_hvscr = [];
%             tc_crvsh = [];
            end   
            
            tc_hvscr_all = cat(2,tc_hvscr_all,mean(resp_h,3)); 
            tc_crvsh_all = cat(2,tc_crvsh_all,mean(resp_cr,3));  
            
            
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
%             tc_hvsm = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
%             tc_mvsh = NaN(length(tt),length(cell_ind),sum(nDirs_hvsm(nDirs_hvsm > 1)));
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
%             if iboot == 1
%                tc_hvsm = resp_rTh;
%                tc_mvsh = resp_rTm;
%             end
            end
            resp_h = mean(resp_h,4);
            resp_m = mean(resp_m,4);
            else
            resp_h = [];
            resp_m = [];
%             tc_hvsm = [];
%             tc_mvsh = [];
            end
            
            
            tc_hvsm_all = cat(2,tc_hvsm_all,mean(resp_h,3)); 
            tc_mvsh_all = cat(2,tc_mvsh_all,mean(resp_m,3)); 
            
            
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
%             tc_favscr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
%             tc_crvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
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
%             if iboot == 1
%                tc_favscr = resp_rTfa;
%                tc_crvsfa = resp_rTcr;
%             end
            end
            resp_fa = mean(resp_fa,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_fa = [];
            resp_cr = [];
%             tc_favscr = [];
%             tc_crvsfa = [];
            end
            
            tc_favscr_all = cat(2,tc_favscr_all,mean(resp_fa,3)); 
            tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(resp_cr,3)); 
            
            
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
            
%         ind = find(lt(nDirs_cr,nDirs_fa));
%             if length(ind) == length(cdirs)
%                 nDirs_favscr = nDirs_cr;
%             elseif isempty(ind)
%                 nDirs_favscr = nDirs_fa;
%             else
%                 ind2 = setdiff(1:length(cdirs),ind);
%                 [j ind3] = sort([ind ind2]);
%                 nT = [nDirs_cr(ind) nDirs_fa(ind2)];
%                 nDirs_favscr = nT(ind3);
%             end
%             if any(nDirs_favscr > 1)
%             resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);
%             resp_cr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);           
% %             tc_favscr = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
% %             tc_crvsfa = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)));
%             for iboot = 1:100
%                 resp_rTfa = [];
%                 resp_rTcr = [];
%             for idir = 1:length(cdirs)
%                 if nDirs_favscr(idir) > 1
%                     rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
%                     resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
%                     rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
%                     resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
%                 end
%             end
%             resp_fa(:,:,:,iboot) = resp_rTfa;
%             resp_cr(:,:,:,iboot) = resp_rTcr;
% %             if iboot == 1
% %                tc_favscr = resp_rTfa;
% %                tc_crvsfa = resp_rTcr;
% %             end
%             end
%             resp_fa = mean(resp_fa,4);
%             resp_cr = mean(resp_cr,4);
%             else
%             resp_fa = [];
%             resp_cr = [];
% %             tc_favscr = [];
% %             tc_crvsfa = [];
%             end
%             
%             tc_favscr_all = cat(2,tc_favscr_all,mean(resp_fa,3)); 
%             tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(resp_cr,3));  
%             
%             if any(nDirs_favscr > 1)
%                 resp_favscr = squeeze(mean(mean(resp_fa(trans_win,:,:),3),1))-squeeze(mean(mean(resp_fa(pre_win,:,:),3),1));
%                 resp_crvsfa = squeeze(mean(mean(resp_cr(trans_win,:,:),3),1))-squeeze(mean(mean(resp_cr(pre_win,:,:),3),1));
%             else
%                 resp_favscr = [];
%                 resp_crvsfa = [];
%             end
% 
% %             resp_favscr = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).resp(pre_win,cell_ind,:),3),1));
% %             resp_crvsfa = squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(trans_win,cell_ind,:),3),1))-squeeze(mean(mean(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).resp(pre_win,cell_ind,:),3),1));
%             resp_favscr_all = cat(2,resp_favscr_all,resp_favscr);
%             resp_crvsfa_all = cat(2,resp_crvsfa_all,resp_crvsfa);
%             n_favscr(i) = sum(nDirs_favscr(nDirs_favscr > 1));
            
        end
        i = i+1;
    end
end

%% target response heatmaps
avgCellRespHeatMap_Val_Inv

%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% valid vs invalid hits, by direction (TC)
mincells = 5;
ncells_h_dir = cellfun(@(resp_h_dir) size(resp_h_dir,2),resp_h_dir,'unif',false);
ncells_m_dir = cellfun(@(resp_m_dir) size(resp_m_dir,2),resp_m_dir,'unif',false);
ncells_all = cellfun(@(val_all_dir) size(val_all_dir,2),val_all_dir,'unif',false);

dir_str = cellfun(@num2str,num2cell(cds),'unif',false);
ncells_h_str = cellfun(@num2str,ncells_h_dir,'unif',false)';
ncells_m_str = cellfun(@num2str,ncells_m_dir,'unif',false)';
ncells_all_str = cellfun(@num2str,ncells_all,'unif',false)';
ntrials_h_str = cellfun(@num2str,num2cell(sum(n_h_dir,2)),'unif',false)';
ntrials_m_str = cellfun(@num2str,num2cell(sum(n_m_dir,2)),'unif',false)';
ntrials_all_str = cellfun(@num2str,num2cell(sum(n_h_dir,2)+sum(n_m_dir,2)),'unif',false)';

colors_val = brewermap(length(cds)+2,'Greys');
colors_inv = brewermap(length(cds)+2,'Blues');
tarRespTraceDir = figure;
subplot(3,2,1)
for idir = 1:length(cds)
    if size(resp_h_dir{idir},2) > mincells
    tc = nanmean(resp_h_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_val(idir+1,:);
    hold on
    end
end
title('valid hits')
ind = find(cell2mat(cellfun(@(resp_h_dir) size(resp_h_dir,2) > 5 ,resp_h_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_h_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_h_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
subplot(3,2,2)
for idir = 1:length(cds)
    if size(resp_fa_dir{idir},2) > mincells
    tc = nanmean(resp_fa_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_inv(idir+1,:);
    hold on
    end
end
title('invalid hits')
ind = find(cell2mat(cellfun(@(resp_fa_dir) size(resp_fa_dir,2) > 5,resp_fa_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_h_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_h_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
subplot(3,2,3)
for idir = 1:length(cds)
    if size(resp_m_dir{idir},2) > mincells
    tc = nanmean(resp_m_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_val(idir+1,:);
    hold on
    end
end
title('valid miss')
ind = find(cell2mat(cellfun(@(resp_m_dir) size(resp_m_dir,2) > 5,resp_m_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_m_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_m_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
subplot(3,2,4)
for idir = 1:length(cds)
    if size(resp_cr_dir{idir},2) > mincells
    tc = nanmean(resp_cr_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_inv(idir+1,:);
    hold on
    end
end
title('invalid miss')
ind = find(cell2mat(cellfun(@(resp_cr_dir) size(resp_cr_dir,2) > 5,resp_cr_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_m_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_m_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
%combo hits and misses
subplot(3,2,5)
for idir = 1:length(cds)
    if size(val_all_dir{idir},2) > mincells
    tc = nanmean(val_all_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_val(idir+1,:);
    hold on
    end
end
title('all valid')
ind = find(cell2mat(cellfun(@(val_all_dir) size(val_all_dir,2) > 5,val_all_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_all_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_all_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
subplot(3,2,6)
for idir = 1:length(cds)
    if size(inv_all_dir{idir},2) > mincells
    tc = nanmean(inv_all_dir{idir,:},2);
    H = plot(ttMs,tc-mean(tc(pre_event_frames,:)));
    H.Color = colors_inv(idir+1,:);
    hold on
    end
end
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
title('all valid')
ind = find(cell2mat(cellfun(@(inv_all_dir) size(inv_all_dir,2) > 5,inv_all_dir,'unif',false)));
legend(strcat(char(dir_str(ind)),char(repmat({'-'},1,length(ind))),char(ncells_all_str(ind)),char(repmat({' c;'},1,length(ind))),char(ntrials_all_str(ind)),char(repmat({' tr'},1,length(ind)))),'location','northwest')
for iplot = 1:6
figure(tarRespTraceDir)
subplot(3,2,iplot)
hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
end

%% tuning curves for valid and invalid trials.
nC_h = cell2mat(ncells_h_dir);
nC_m = cell2mat(ncells_m_dir);
nC_all = cell2mat(ncells_all);

ind = nC_h  > mincells;
tcH = resp_h_dir(ind);
rH = cellfun(@(tcH) mean(tcH(trans_win,:),1),tcH,'unif',false);
bH = cellfun(@(tcH) mean(tcH(pre_win,:),1),tcH,'unif',false);
errH = cell2mat(cellfun(@(rH) std(rH),rH,'unif',false))./sqrt(nC_h(nC_h > mincells));
tcFA = resp_fa_dir(ind);
rFA = cellfun(@(tcFA) mean(tcFA(trans_win,:),1),tcFA,'unif',false);
bFA = cellfun(@(tcFA) mean(tcFA(pre_win,:),1),tcFA,'unif',false);
errFA = cell2mat(cellfun(@(rFA) std(rFA),rFA,'unif',false))./sqrt(nC_h(nC_h > mincells));

ind = nC_m > mincells;
tcM = resp_m_dir(ind);
rM = cellfun(@(tcM) mean(tcM(trans_win,:),1),tcM,'unif',false);
bM = cellfun(@(tcM) mean(tcM(pre_win,:),1),tcM,'unif',false);
errM = cell2mat(cellfun(@(rM) std(rM),rM,'unif',false))./sqrt(nC_m(nC_m > mincells));
tcCR = resp_cr_dir(ind);
rCR = cellfun(@(tcCR) mean(tcCR(trans_win,:),1),tcCR,'unif',false);
bCR = cellfun(@(tcCR) mean(tcCR(pre_win,:),1),tcCR,'unif',false);
errCR = cell2mat(cellfun(@(rCR) std(rCR),rCR,'unif',false))./sqrt(nC_m(nC_m > mincells));

ind = nC_all > mincells;
tcVal = val_all_dir(ind);
rVal = cellfun(@(tcVal) mean(tcVal(trans_win,:),1),tcVal,'unif',false);
bVal = cellfun(@(tcVal) mean(tcVal(pre_win,:),1),tcVal,'unif',false);
errVal = cell2mat(cellfun(@(rVal) std(rVal),rVal,'unif',false))./sqrt(nC_all(nC_all > mincells));
tcInv = inv_all_dir(ind);
rInv = cellfun(@(tcInv) mean(tcInv(trans_win,:),1),tcInv,'unif',false);
bInv = cellfun(@(tcInv) mean(tcInv(pre_win,:),1),tcInv,'unif',false);
errInv = cell2mat(cellfun(@(rInv) std(rInv),rInv,'unif',false))./sqrt(nC_all(nC_all > mincells));

tarRespTuning = figure;
subplot(3,2,1)
ind = nC_h > mincells;
r = cell2mat(cellfun(@mean,rH,'unif',false))-cell2mat(cellfun(@mean,bH,'unif',false));
errorbar(cds(ind),r,errH,'ko');
hold on
r = cell2mat(cellfun(@mean,rFA,'unif',false))-cell2mat(cellfun(@mean,bFA,'unif',false));
errorbar(cds(ind),r,errFA,'co')
title('hits')
subplot(3,2,3)
ind = nC_m > mincells;
r = cell2mat(cellfun(@mean,rM,'unif',false))-cell2mat(cellfun(@mean,bM,'unif',false));
errorbar(cds(ind),r,errM,'ko');
hold on
r = cell2mat(cellfun(@mean,rCR,'unif',false))-cell2mat(cellfun(@mean,bCR,'unif',false));
errorbar(cds(ind),r,errCR,'co')
title('misses')
subplot(3,2,5)
ind = nC_all > mincells;
r = cell2mat(cellfun(@mean,rVal,'unif',false))-cell2mat(cellfun(@mean,bVal,'unif',false));
errorbar(cds(ind),r,errVal,'ko');
hold on
r = cell2mat(cellfun(@mean,rInv,'unif',false))-cell2mat(cellfun(@mean,bInv,'unif',false));
errorbar(cds(ind),r,errInv,'co')
title('all')

for iplot = [1 3 5]
figure(tarRespTuning)
subplot(3,2,iplot)
hold on
xlim([0 100])
ylim([-0.01 0.03])
ax = gca;
ax.XTick = cds;
end


%% scatters with all cells
scatAllCells = figure;
suptitle([titleStr ', pairwise ttest'])

[h,p] = ttest(resp_hvsfa_all,resp_favsh_all,'alpha',0.05/(size(resp_favsh_all,2)));
subplot(3,2,1)
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
title(['H-v,H-inv; p=' num2str(p) ', ' num2str(size(resp_favsh_all,2)) '  cells'])

[h,p] = ttest(resp_mvscr_all,resp_crvsm_all,'alpha',0.05/(size(resp_crvsm_all,2)));
subplot(3,2,2)
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
title(['M-v,M-inv; p=' num2str(p) ', ' num2str(size(resp_crvsm_all,2)) '  cells'])

[h,p] = ttest(resp_hvscr_all,resp_crvsh_all,'alpha',0.05/(size(resp_crvsh_all,2)));
subplot(3,2,3)
scatter(resp_hvscr_all,resp_crvsh_all,50,'k.')
hold on
errorbarxy(mean(resp_hvscr_all),mean(resp_crvsh_all),std(resp_hvscr_all)/length(resp_hvscr_all),std(resp_crvsh_all)/length(resp_crvsh_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid hits')
ylabel('invalid hits')
title(['H-v,M-inv; p=' num2str(p) ', ' num2str(size(resp_crvsh_all,2)) '  cells'])

[h,p] = ttest(resp_hvsm_all,resp_mvsh_all,'alpha',0.05/(size(resp_mvsh_all,2)));
subplot(3,2,4)
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
title(['H-v,M-v; p=' num2str(p) ', ' num2str(size(resp_mvsh_all,2)) '  cells'])

[h,p] = ttest(resp_favscr_all,resp_crvsfa_all,'alpha',0.05/(size(resp_crvsfa_all,2)));
subplot(3,2,5)
scatter(resp_favscr_all,resp_crvsfa_all,50,'k.')
hold on
errorbarxy(mean(resp_favscr_all),mean(resp_crvsfa_all),std(resp_favscr_all)/length(resp_favscr_all),std(resp_crvsfa_all)/length(resp_crvsfa_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('invalid hits')
ylabel('invalid miss')
title(['H-inv,M-inv; p=' num2str(p) ', ' num2str(size(resp_crvsfa_all,2)) '  cells'])

[h,p] = ttest(resp_val_all,resp_inv_all,'alpha',0.05/(size(resp_inv_all,2)));
subplot(3,2,6)
scatter(resp_val_all,resp_inv_all,50,'k.')
hold on
errorbarxy(mean(resp_val_all),mean(resp_inv_all),std(resp_val_all)/length(resp_val_all),std(resp_inv_all)/length(resp_inv_all),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid all')
ylabel('invalid all')
title(['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells'])

figure(scatAllCells)
print([fnout 'catch_align_scatter' datasetStr '.pdf'], '-dpdf')
%% all directions -  scatter, cdfs, tc
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

figure(tarRespTuning)
subplot(3,2,2)
r = mean(resp_hvsfa_all);
errVal = std(resp_hvsfa_all)/sqrt(length(resp_hvsfa_all));
errorbar(mean(cds),r,errVal,'ko');
hold on
r = mean(resp_favsh_all);
errInv = std(resp_favsh_all)/sqrt(length(resp_favsh_all));
errorbar(mean(cds),r,errInv,'co')
xlim([0 100])
ylim([-0.01 0.03])


tcCatchResp = figure;
suptitle([mouse(1).expt(1).cells(cellsInd).name ' pref'])
subplot(3,2,1)
tc1 = mean(tc_hvsfa_all,2)-mean(mean(tc_hvsfa_all(pre_event_frames,:),2));
err1 = std(tc_hvsfa_all,[],2)/sqrt(size(tc_hvsfa_all,2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = mean(tc_favsh_all,2)-mean(mean(tc_favsh_all(pre_event_frames,:),2));
err2 = std(tc_favsh_all,[],2)/sqrt(size(tc_favsh_all,2));
shadedErrorBar(ttMs,tc2,err2,'c')
hold on
% plot(ttMs,mean(tc_hvsfa_all,2)-mean(mean(tc_hvsfa_all(pre_event_frames,:),2)),'k')
% hold on
% plot(ttMs,mean(tc_favsh_all,2)-mean(mean(tc_favsh_all(pre_event_frames,:),2)),'c')
% hold on
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

figure(tarRespTuning)
subplot(3,2,4)
r = mean(resp_mvscr_all);
errVal = std(resp_mvscr_all)/sqrt(length(resp_mvscr_all));
errorbar(mean(cds),r,errVal,'ko');
hold on
r = mean(resp_crvsm_all);
errInv = std(resp_crvsm_all)/sqrt(length(resp_crvsm_all));
errorbar(mean(cds),r,errInv,'co')
xlim([0 100])
ylim([-0.01 0.03])

figure(tcCatchResp);
subplot(3,2,2)
tc1 = mean(tc_mvscr_all,2)-mean(mean(tc_mvscr_all(pre_event_frames,:),2));
err1 = std(tc_mvscr_all,[],2)/sqrt(size(tc_mvscr_all,2));
shadedErrorBar(ttMs,tc1,err1,'r');
hold on
tc2 = mean(tc_crvsm_all,2)-mean(mean(tc_crvsm_all(pre_event_frames,:),2));
err2 = std(tc_crvsm_all,[],2)/sqrt(size(tc_crvsm_all,2));
shadedErrorBar(ttMs,tc2,err2,'b');
hold on
% plot(ttMs,mean(tc_mvscr_all,2)-mean(mean(tc_mvscr_all(pre_event_frames,:),2)),'r')
% hold on
% plot(ttMs,mean(tc_crvsm_all,2)-mean(mean(tc_crvsm_all(pre_event_frames,:),2)),'b')
% hold on
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
tc1 = mean(tc_hvscr_all,2)-mean(mean(tc_hvscr_all(pre_event_frames,:),2));
err1 = std(tc_hvscr_all,[],2)/sqrt(size(tc_hvscr_all,2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = mean(tc_crvsh_all,2)-mean(mean(tc_crvsh_all(pre_event_frames,:),2));
err2 = std(tc_crvsh_all,[],2)/sqrt(size(tc_crvsh_all,2));
shadedErrorBar(ttMs,tc2,err2,'b');
hold on
% plot(ttMs,mean(tc_hvscr_all,2)-mean(mean(tc_hvscr_all(pre_event_frames,:),2)),'k')
% hold on
% plot(ttMs,mean(tc_crvsh_all,2)-mean(mean(tc_crvsh_all(pre_event_frames,:),2)),'b')
% hold on
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
tc1 = mean(tc_hvsm_all,2)-mean(mean(tc_hvsm_all(pre_event_frames,:),2));
err1 = std(tc_hvsm_all,[],2)/sqrt(size(tc_hvsm_all,2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = mean(tc_mvsh_all,2)-mean(mean(tc_mvsh_all(pre_event_frames,:),2));
err2 = std(tc_mvsh_all,[],2)/sqrt(size(tc_mvsh_all,2));
shadedErrorBar(ttMs,tc2,err2,'r');
hold on
% plot(ttMs,mean(tc_hvsm_all,2)-mean(mean(tc_hvsm_all(pre_event_frames,:),2)),'k')
% hold on
% plot(ttMs,mean(tc_mvsh_all,2)-mean(mean(tc_mvsh_all(pre_event_frames,:),2)),'r')
% hold on
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
tc1 = mean(tc_favscr_all,2)-mean(mean(tc_favscr_all(pre_event_frames,:),2));
err1 = std(tc_favscr_all,[],2)/sqrt(size(tc_favscr_all,2));
shadedErrorBar(ttMs,tc1,err1,'c');
hold on
tc2 = mean(tc_crvsfa_all,2)-mean(mean(tc_crvsfa_all(pre_event_frames,:),2));
err2 = std(tc_crvsfa_all,[],2)/sqrt(size(tc_crvsfa_all,2));
shadedErrorBar(ttMs,tc2,err2,'b')
hold on
% plot(ttMs,mean(tc_favscr_all,2)-mean(mean(tc_favscr_all(pre_event_frames,:),2)),'c')
% hold on
% plot(ttMs,mean(tc_crvsfa_all,2)-mean(mean(tc_crvsfa_all(pre_event_frames,:),2)),'b')
% hold on
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

figure(tarRespTuning)
subplot(3,2,6)
r = mean(resp_val_all);
errVal = std(resp_val_all)/sqrt(length(resp_val_all));
errorbar(mean(cds),r,errVal,'ko');
hold on
r = mean(resp_inv_all);
errInv = std(resp_inv_all)/sqrt(length(resp_inv_all));
errorbar(mean(cds),r,errInv,'co')
xlim([0 100])
ylim([-0.01 0.03])


figure(tcCatchResp)
subplot(3,2,6)
tc1 = mean(tc_val_all,2)-mean(mean(tc_val_all(pre_win,:),2));
err1 = std(tc_val_all,[],2)/sqrt(size(tc_val_all,2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = mean(tc_inv_all,2)-mean(mean(tc_inv_all(pre_win,:),2));
err2 = std(tc_inv_all,[],2)/sqrt(size(tc_inv_all,2));
shadedErrorBar(ttMs,tc2,err2,'c');
hold on
% plot(ttMs,mean(tc_val_all,2)-mean(mean(tc_val_all(pre_win,:),2)),'k');
% hold on
% plot(ttMs,mean(tc_inv_all,2)-mean(mean(tc_inv_all(pre_win,:),2)),'c');
% hold on
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames,':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title([num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']);

%% resp size by trial length - 90 hits only
L = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        if mouse(imouse).expt(iexp).info.isCatch 
        cycTimeMs = mouse(imouse).expt(iexp).info.cyc_time_ms;
        L = [L mouse(imouse).expt(iexp).align(catchAlign).av(1).outcome(1).tcyc{end}*cycTimeMs];
        end
    end
end
L = unique(L);
Lshort = L(L < 2500);
Llong = L(L > 2500);
        
rVal_short = [];
rVal_long = [];
rInv_short = [];
rInv_long = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2) 
        if mouse(imouse).expt(iexp).info.isCatch
        % find 90 target index
        cdirs = mouse(imouse).expt(iexp).info.cDirs;
        ind = find(cdirs == 90);
        cycTimeMs = mouse(imouse).expt(iexp).info.cyc_time_ms;
        if ~isempty(ind)
            rVal = mouse(imouse).expt(iexp).align(catchAlign).av(1).outcome(hits).stimResp{ind};
            rInv = mouse(imouse).expt(iexp).align(catchAlign).av(1).outcome(fas).stimResp{ind};
            trL_val = mouse(imouse).expt(iexp).align(catchAlign).av(1).outcome(hits).tcyc{ind}*cycTimeMs;
            trL_inv = mouse(imouse).expt(iexp).align(catchAlign).av(1).outcome(fas).tcyc{ind}*cycTimeMs;
            
            if sum(sum(isnan(mean(rInv(:,:,find(trL_inv < 2500)),3)))) == 0 & sum(sum(isnan(mean(rInv(:,:,find(trL_inv > 2500)),3)))) == 0
                
                figure;
                hold on
                plot(mean(mean(rInv(:,:,find(trL_inv > 2500)),3),2))
                
            rVal_short = cat(2,rVal_short,mean(rVal(:,:,find(trL_val < 2500)),3));
            rVal_long = cat(2,rVal_short,mean(rVal(:,:,find(trL_val > 2500)),3));
            rInv_short = cat(2,rInv_short,mean(rInv(:,:,find(trL_inv < 2500)),3));
            rInv_long = cat(2,rInv_short,mean(rInv(:,:,find(trL_inv > 2500)),3));
            end
        end
        end
    end
end

%plot tc
trialLength_tarTC = figure;
subplot(3,2,1)
plot(ttMs,nanmean(rVal_short,2),'k')
hold on
plot(ttMs,nanmean(rInv_short,2),'c')
hold on
xlabel('time (ms)')
ylabel('dF/F')
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title('trials < 2.5s; 90deg hits')

subplot(3,2,2)
plot(ttMs,nanmean(rVal_long,2),'k')
hold on
plot(ttMs,nanmean(rInv_long,2),'c')
hold on
xlabel('time (ms)')
ylabel('dF/F')
xlim([-10 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
title('trials > 2.5s')

subplot(3,2,3)
r = mean(nanmean(rVal_short(trans_win,:),2),1)-mean(nanmean(rVal_short(pre_win,:),2),1);
plot(mean(Lshort),r,'ko')
hold on
r = mean(nanmean(rInv_short(trans_win,:),2),1)-mean(nanmean(rInv_short(pre_win,:),2),1);
plot(mean(Lshort),r,'co')
hold on
r = mean(nanmean(rVal_long(trans_win,:),2),1)-mean(nanmean(rVal_long(pre_win,:),2),1);
plot(mean(Llong),r,'ko')
hold on
r = mean(nanmean(rInv_long(trans_win,:),2),1)-mean(nanmean(rInv_long(pre_win,:),2),1);
plot(mean(Llong),r,'co')


%% 
% tuning2catch
%%
% 
% figure(scatFAvsH);
% print([fnout 'catch_align_FAvsH' datasetStr '.pdf'], '-dpdf')
% figure(scatCRvsM);
% print([fnout 'catch_align_CRvsM' datasetStr '.pdf'], '-dpdf')
% figure(scatCRvsH);
% print([fnout 'catch_align_CRvsH' datasetStr '.pdf'], '-dpdf')
% figure(scatMvsH);
% print([fnout 'catch_align_MvsH' datasetStr '.pdf'], '-dpdf')
% figure(scatCRvsFA);
% print([fnout 'catch_align_CRvsFA' datasetStr '.pdf'], '-dpdf')
figure(cdfCatchResp)
print([fnout 'catch_align_cdf' datasetStr '.pdf'], '-dpdf')
figure(tcCatchResp)
print([fnout 'catch_align_TC' datasetStr '.pdf'], '-dpdf')
figure(tarRespTraceDir)
print([fnout 'catch_align_TC_dirs' datasetStr '.pdf'], '-dpdf')
figure(tarRespTuning)
print([fnout 'tarRespTuning_val_inv' datasetStr '.pdf'], '-dpdf')

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


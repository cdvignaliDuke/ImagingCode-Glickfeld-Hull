
cellsCatchAlign = 13;

nbins = 3;
compInd = 1;
aucColors = brewermap(nbins,'*RdBu');
tc_val_lastBaseStim = [];
tc_inv_lastBaseStim = [];
dscCellsCount = 0;
start = 1;

auc_fig = figure;
colormap(brewermap([],'*RdBu'));
auc_fig_press_early = figure;
colormap(brewermap([],'*RdBu'));
auc_fig_press_late = figure;
colormap(brewermap([],'*RdBu'));
heatmap_fig = figure;
colormap(brewermap([],'*RdBu'));
valOverlay = figure;
colormap(brewermap([],'*RdBu'));
clear legData;

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

%pre-allocate press align tcs

tc_vis_late = [];
tc_aud_late = [];
tc_vis_early = [];
tc_aud_early = [];

oriTuningResp_avg = [];
oriTuningResp_tc = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)        
        if mouse(imouse).expt(iexp).info.isCatch
            val_temp = [];
            inv_temp = [];
            all_temp = [];
            cell_ind = mouse(imouse).expt(iexp).cells(cellsCatchAlign).ind;

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
%             resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa),100);           
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            for iboot = 1:100
                resp_rTh = NaN(length(tt),length(cell_ind), sum(nDirs_hvsfa(nDirs_hvsfa > 1)));
                resp_rTfa = NaN(length(tt),length(cell_ind), sum(nDirs_hvsfa(nDirs_hvsfa > 1)));
%                 rH_dir_temp = cell(length(cdirs),1);
%                 rFA_dir_temp = cell(length(cdirs),1);
                begin = 1;
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                dirsInd = begin:begin-1+nDirs_hvsfa(idir);
                begin = begin+nDirs_hvsfa(idir);
                    if lt(nDirs_fa(idir),nDirs_h(idir))
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                    cT = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).tcyc{idir}));
                    ncT = histc(double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).tcyc{idir}),cT);
                    for ict = 1:length(cT)
                        cycInd = find(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{idir} == cT(ict));
                        if ~isempty(cycInd) & ncT(ict) > 1 
                            if ncT(ict) < nDirs_hvsfa(idir)
                            resp_rTh(:,:,dirsInd(1:ncT(ict))) = rT(:,cell_ind,randsample(cycInd,ncT(ict),1));
                            else
                            resp_rTh(:,:,dirsInd) = rT(:,cell_ind,randsample(cycInd,nDirs_hvsfa(idir),1));
                            end
                        end
                    end
                        rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                        resp_rTfa(:,:,dirsInd) = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    else
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).stimResp{idir};
                    cT = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{idir}));
                    ncT = histc(double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).tcyc{idir}),cT);
                    for ict = 1:length(cT)
                        cycInd = find(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(fas).tcyc{idir} == cT(ict));
                        if ~isempty(cycInd) & ncT(ict) > 1  
                            if ncT(ict) < nDirs_hvsfa(idir)
                            resp_rTfa(:,:,dirsInd(1:ncT(ict))) = rT(:,cell_ind,randsample(cycInd,ncT(ict),1));
                            else
                            resp_rTfa(:,:,dirsInd) = rT(:,cell_ind,randsample(cycInd,nDirs_hvsfa(idir),1));
                            end
                        end
                    end
                        rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(hits).stimResp{idir};
                        resp_rTh(:,:,dirsInd) = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    end
%                     rH_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
%                     resp_h_dir_temp{idir} = cat(4,resp_h_dir_temp{idir},rH_dir_temp{idir});
%                     rFA_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
%                     resp_fa_dir_temp{idir} = cat(4,resp_fa_dir_temp{idir},rFA_dir_temp{idir});
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
%             resp(fas).all(:,:,:,iboot) = resp_rTfa;
%             if iboot == 1
%                tc_hvsfa = resp_rTh;
%                tc_favsh = resp_rTfa;
%             end
            end
            resp_h = nanmean(resp_h,4);
            resp(fas).all = nanmean(resp(fas).all,4);
            
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
            resp(fas).all = [];
            if iAuc == 3
            dscCellsCount = dscCellsCount+length(cell_ind);
            disp([mouse(imouse).expt(iexp).date '-h' num2str(length(mouse(imouse).expt(iexp).cells(cellInd_auROC).ind))])
            end
%             tc_hvsfa = [];
%             tc_favsh = [];
            end
            
            
            tc_hvsfa_all = cat(2,tc_hvsfa_all,nanmean(resp_h,3));%tc_hvsfa_all = cat(2,tc_hvsfa_all,mean(tc_hvsfa,3)); 
            tc_favsh_all = cat(2,tc_favsh_all,nanmean(resp(fas).all,3));%tc_favsh_all = cat(2,tc_favsh_all,mean(tc_favsh,3));  
            
            if sum(sum(isnan(resp_h))) > 0 
                disp('resp_h nan after hvsfa')
            end
            
            if any(nDirs_hvsfa > 1)
            resp_hvsfa = squeeze(nanmean(nanmean(resp_h(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(resp_h(pre_win,:,:),3),1));
            resp_favsh = squeeze(nanmean(nanmean(resp(fas).all(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(resp(fas).all(pre_win,:,:),3),1));
            val_temp = resp_h;
            inv_temp = resp(fas).all;
            all_temp = cat(3,resp_h,resp(fas).all);
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
%             tc_inv_all = cat(2,tc_inv_all,mean(resp(fas).all,3));%tc_inv_all = cat(2,tc_inv_all,mean(tc_favsh,3));
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
                resp_rTm = NaN(length(tt),length(cell_ind), sum(nDirs_mvscr(nDirs_mvscr > 1)));
                resp_rTcr = NaN(length(tt),length(cell_ind), sum(nDirs_mvscr(nDirs_mvscr > 1)));
%                 rM_dir_temp = cell(length(cdirs),1);
%                 rCR_dir_temp = cell(length(cdirs),1);
                begin = 1;
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                dirsInd = begin:begin-1+nDirs_mvscr(idir);
                begin = begin+nDirs_mvscr(idir);
                    if lt(nDirs_cr(idir),nDirs_m(idir))
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                    cT = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{idir}));
                    ncT = histc(double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{idir}),cT);
                    for ict = 1:length(cT)
                        cycInd = find(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).tcyc{idir} == cT(ict));
                        if ~isempty(cycInd) & ncT(ict) > 1 
                            if ncT(ict) < nDirs_mvscr(idir)
                            resp_rTm(:,:,dirsInd(1:ncT(ict))) = rT(:,cell_ind,randsample(cycInd,ncT(ict),1));
                            else
                            resp_rTm(:,:,dirsInd) = rT(:,cell_ind,randsample(cycInd,nDirs_mvscr(idir),1));
                            end
                        end
                    end
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    resp_rTcr(:,:,dirsInd) = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    else
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).stimResp{idir};
                    cT = double(unique(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).tcyc{idir}));
                    ncT = histc(double(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).tcyc{idir}),cT);
                    for ict = 1:length(cT)
                        cycInd = find(mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(crs).tcyc{idir} == cT(ict));
                        if ~isempty(cycInd) & ncT(ict) > 1
                            if ncT(ict) < nDirs_mvscr(idir)
                            resp_rTcr(:,:,dirsInd(1:ncT(ict))) = rT(:,cell_ind,randsample(cycInd,ncT(ict),1));
                            else
                            resp_rTcr(:,:,dirsInd) = rT(:,cell_ind,randsample(cycInd,nDirs_mvscr(idir),1));
                            end
                        end
                    end
                        rT = mouse(imouse).expt(iexp).align(catchAlign).av(iav).outcome(misses).stimResp{idir};
                        resp_rTm(:,:,dirsInd) = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));   
                    end
%                     rM_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
%                     resp_m_dir_temp{idir} = cat(4,resp_m_dir_temp{idir},rM_dir_temp{idir});
%                     rCR_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
%                     resp_cr_dir_temp{idir} = cat(4,resp_cr_dir_temp{idir},rCR_dir_temp{idir});
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
            if iAuc == 3 & any(nDirs_hvsfa > 1)
                dscCellsCount = dscCellsCount+length(cell_ind);
                disp([mouse(imouse).expt(iexp).date '-m' num2str(length(mouse(imouse).expt(iexp).cells(cellInd_auROC).ind))])
            end
            end
            
            tc_mvscr_all = cat(2,tc_mvscr_all,nanmean(resp_m,3)); 
            tc_crvsm_all = cat(2,tc_crvsm_all,nanmean(resp_cr,3)); 
            
            if any(nDirs_mvscr > 1)
                    resp_mvscr = squeeze(nanmean(nanmean(resp_m(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(resp_m(pre_win,:,:),3),1));
                resp_crvsm = squeeze(nanmean(nanmean(resp_cr(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(resp_cr(pre_win,:,:),3),1));
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
            
            if sum(isnan(val_temp)) > 0 
                disp('val_temp nan after mvscr')
            end
            
            if ~isempty(val_temp) | ~isempty(inv_temp)
            resp_val_all = cat(2,resp_val_all,squeeze(nanmean(nanmean(val_temp(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(val_temp(pre_win,:,:),3),1)));%cat(2,resp_val_all,resp_mvscr);
            resp_inv_all = cat(2,resp_inv_all,squeeze(nanmean(nanmean(inv_temp(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(inv_temp(pre_win,:,:),3),1)));%cat(2,resp_inv_all,resp_crvsm);
            resp_all_all = cat(2,resp_all_all,squeeze(nanmean(nanmean(all_temp(trans_win,:,:),3),1))-squeeze(nanmean(nanmean(all_temp(pre_win,:,:),3),1)));
            tc_val_all = cat(2,tc_val_all,nanmean(val_temp,3));%nanmean(resp_m,3));
            tc_inv_all = cat(2,tc_inv_all,nanmean(inv_temp,3));%nanmean(resp_cr,3));
            tc_all = cat(2,tc_all,nanmean(all_temp,3));
            end
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
            resp(fas).all = NaN(length(tt),length(cell_ind),sum(nDirs_favscr(nDirs_favscr > 1)),100);
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
            resp(fas).all(:,:,:,iboot) = resp_rTfa;
            resp_cr(:,:,:,iboot) = resp_rTcr;
%             if iboot == 1
%                tc_favscr = resp_rTfa;
%                tc_crvsfa = resp_rTcr;
%             end
            end
            resp(fas).all = mean(resp(fas).all,4);
            resp_cr = mean(resp_cr,4);
            else
            resp(fas).all = [];
            resp_cr = [];
%             tc_favscr = [];
%             tc_crvsfa = [];
            end
            
            tc_favscr_all = cat(2,tc_favscr_all,mean(resp(fas).all,3)); 
            tc_crvsfa_all = cat(2,tc_crvsfa_all,mean(resp_cr,3)); 
            
            
            if any(nDirs_favscr > 1)
                resp_favscr = squeeze(mean(mean(resp(fas).all(trans_win,:,:),3),1))-squeeze(mean(mean(resp(fas).all(pre_win,:,:),3),1));
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
            
%********* % press align tc
            cycFrames = mouse(imouse).expt(iexp).info.cyc_time;
            fr_temp = (cycFrames/mouse(imouse).expt(iexp).info.cyc_time_ms);
            minLongMs = 2800;
            minShortMs = 1400;
            minLongFrames = floor(minLongMs*fr_temp);
            minShortFrames = floor(minShortMs*fr_temp);
            cycsLong = ceil(minLongFrames/cycFrames);
            cycsShort = ceil(minShortFrames/cycFrames);
            
            if size(mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp,2) >= cycsLong
            v_long_temp = mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp{cycsLong};
            a_long_temp = mouse(imouse).expt(iexp).align(1).av(2).outcome(hits).cmlvCycResp{cycsLong};
            
            v_sht_temp = mouse(imouse).expt(iexp).align(1).av(1).outcome(hits).cmlvCycResp{cycsShort};
            a_sht_temp = mouse(imouse).expt(iexp).align(1).av(2).outcome(hits).cmlvCycResp{cycsShort};
            
            tc_vis_late = cat(2,tc_vis_late,mean(v_long_temp(1:minLongFrames+pre_event_frames,cell_ind,:),3));
            tc_aud_late = cat(2,tc_aud_late,mean(a_long_temp(1:minLongFrames+pre_event_frames,cell_ind,:),3));
            tc_vis_early = cat(2,tc_vis_early,mean(v_sht_temp(1:minShortFrames+pre_event_frames,cell_ind,:),3));
            tc_aud_early = cat(2,tc_aud_early,mean(a_sht_temp(1:minShortFrames+pre_event_frames,cell_ind,:),3));
            
            end
       

          
        end
        i = i+1;
    end
end

%scatter
figure(auc_fig);
[h,p] = ttest(resp_val_all,resp_inv_all,'alpha',0.05/(size(resp_inv_all,2)));
subplot(nbins,3,1+((iAuc-1)*3))
scatter(resp_val_all,resp_inv_all,50,'k.')
hold on
errorbarxy(nanmean(resp_val_all),nanmean(resp_inv_all),nanstd(resp_val_all(~isnan(resp_val_all)))/length(resp_val_all),nanstd(resp_inv_all)/length(resp_inv_all(~isnan(resp_inv_all))),{'ro','r','r'});
xlim([-0.05 0.1])
ylim([-0.05 0.1])
hold on
plot([-10:0.1:20],[-10:0.1:20],'k--')
axis square
xlabel('valid all')
ylabel('invalid all')
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells']})
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];['all val,all inv; p=' num2str(p) ', ' num2str(size(resp_inv_all,2)) '  cells']})

%cdf
figure(auc_fig);
subplot(nbins,3,2+((iAuc-1)*3))
hVvsInv = cdfplot(resp_val_all);
hVvsInv.Color = 'k';
hold on
hInvvsV = cdfplot(resp_inv_all);
hInvvsV.Color = 'c';
[h p] = kstest2(resp_val_all,resp_inv_all);
xlim([-0.05 0.05])
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
axis square

%time-courses
    %create tc across all trials, aligned to previous base stim
    nCells = size(tc_val_all,2);
    tc_val_lastBaseAlign(:,start:(start+nCells)-1) = tc_val_all;
    tc_inv_lastBaseAlign(:,start:(start+nCells)-1) = tc_inv_all;
    start = start+nCells;

figure(auc_fig)
subplot(nbins,3,3+((iAuc-1)*3))
tc1 = nanmean(tc_val_all,2)-nanmean(nanmean(tc_val_all(pre_win,:),2));
err1 = nanstd(tc_val_all,[],2)/sqrt(size(tc_val_all(~isnan(sum(tc_inv_all))),2));
shadedErrorBar(ttMs,tc1,err1,'k');
hold on
tc2 = nanmean(tc_inv_all,2)-nanmean(nanmean(tc_inv_all(pre_win,:),2));
err2 = nanstd(tc_inv_all,[],2)/sqrt(size(tc_inv_all(~isnan(sum(tc_inv_all))),2));
shadedErrorBar(ttMs,tc2,err2,'c');
hold on
xlim([-15 20]/(cycTime/cycTimeMs))
ylim([-0.01 0.05])
vline(baseStimFrames/(cycTime/cycTimeMs),':k')
hold on
vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs),'--r')
hold on
h = vline(([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs))-cycTimeMs,'--r');
h(1).Color = [0.5 0 0];
h(2).Color = [0.5 0 0];
hold on
vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
hold on
vline(([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs))-cycTimeMs, '--k')
% title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];[num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']});
title({['ROC bin ' num2str(iAuc) ', ' perc_ind_maxDir(iAuc) ':' perc_ind_maxDir(iAuc+1)];[num2str(sum(n_all)) ' all val, inv;' num2str(size(resp_val_all,2)) ' cells']});
axis square


figure(valOverlay)
subplot(1,2,1)
hold on
l = shadedErrorBar(ttMs,tc1,err1);
l.mainLine.Color = aucColors(iAuc,:);
l.mainLine.LineWidth = 3;
legData(iAuc) = l.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('valid trials')
legendInfo{iAuc} = ['Bin ' num2str(iAuc)];
subplot(1,2,2)
hold on
l = shadedErrorBar(ttMs,tc2,err2);
l.mainLine.Color = aucColors(iAuc,:);
l.mainLine.LineWidth = 3;
legData(iAuc) = l.mainLine;
hold on
vline(baseStimFrames/(cycTime/cycTimeMs),'k:')
xlim([-1000 1000])
ylim([-0.01 0.05])
xlabel('time (ms)')
ylabel('dF/F')
title('invalid trials')
legendInfo{iAuc} = ['Bin ' num2str(iAuc)];

figure(heatmap_fig);
subplot(nbins,2,1+((iAuc-1)*2))
temp_hm = bsxfun(@minus,tc_val_all,nanmean(tc_val_all(pre_win,:),1));
% hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*3000),:)');
% hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*2000)];
% hm.Parent.XTickLabel = [0 1 2]
hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000),:)');
hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*250) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*500)];
hm.Parent.XTickLabel = [0 0.25 0.5];
xlabel('time (s)')
ylabel('cells')
title('valid')
axis square
clim([-0.1 0.1])
colorbar
subplot(nbins,2,2+((iAuc-1)*2))
temp_hm = bsxfun(@minus,tc_inv_all,nanmean(tc_inv_all(pre_win,:),1));
% hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*3000),:)');
% hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*2000)];
% hm.Parent.XTickLabel = [0 1 2]
hm = imagesc(temp_hm(20:(pre_event_frames-20)+floor((cycTime/cycTimeMs)*1000),:)');
hm.Parent.XTick = [pre_event_frames-20 (pre_event_frames-20)+floor((cycTime/cycTimeMs)*250) (pre_event_frames-20)+floor((cycTime/cycTimeMs)*500)];
hm.Parent.XTickLabel = [0 0.25 0.5];
xlabel('time (s)')
ylabel('cells')
title('invalid')
axis square
clim([-0.1 0.1])
colorbar
% % 
% % %*********%press align figures
% % early_win = (11:33)+trans_win(1);
% % late_win = ((66:88)+pre_event_frames)-(trans_win(1)-pre_event_frames);
% % resp_v_late = nanmean(nanmean(tc_vis_late(late_win,:),1),3);
% % resp_a_late = nanmean(nanmean(tc_aud_late(late_win,:),1),3);
% % resp_v_early = nanmean(nanmean(tc_vis_early(early_win,:),1),3);
% % resp_a_early = nanmean(nanmean(tc_aud_early(early_win,:),1),3);
% % 
% % 
% % %scatter
% % %early analysis window
% % figure(auc_fig_press_early)
% % [h,p] = ttest(resp_v_early,resp_a_early,'alpha',0.05/(size(resp_a_early,2)));
% % subplot(nbins,3,1+((iAuc-1)*3))
% % scatter(resp_v_early,resp_a_early,50,'k.')
% % hold on
% % errorbarxy(nanmean(resp_v_early),nanmean(resp_a_early),nanstd(resp_v_early)/length(resp_v_early(~isnan(resp_v_early))),nanstd(resp_a_early)/length(resp_a_early(~isnan(resp_a_early)),{'ro','r','r'});
% % xlim([-0.05 0.1])
% % ylim([-0.05 0.1])
% % hold on
% % plot([-10:0.1:20],[-10:0.1:20],'k--')
% % axis square
% % xlabel('visual')
% % ylabel('auditory')
% % title({['ROC bin ' num2str(iAuc)];['early-press align; p=' num2str(p) ', ' num2str(size(resp_v_early,2)) '  cells']})
% % 
% % %late analyssi window
% % figure(auc_fig_press_late)
% % [h,p] = ttest(resp_v_late,resp_a_late,'alpha',0.05/(size(resp_a_late,2)));
% % subplot(nbins,3,1+((iAuc-1)*3))
% % scatter(resp_v_late,resp_a_late,50,'k.')
% % hold on
% % errorbarxy(nanmean(resp_v_late),nanmean(resp_a_late),nanstd(resp_v_late)/length(resp_v_late),nanstd(resp_a_late)/length(resp_a_late),{'ro','r','r'});
% % xlim([-0.05 0.1])
% % ylim([-0.05 0.1])
% % hold on
% % plot([-10:0.1:20],[-10:0.1:20],'k--')
% % axis square
% % xlabel('visual')
% % ylabel('auditory')
% % title({['ROC bin ' num2str(iAuc)];['late-press align; p=' num2str(p) ', ' num2str(size(resp_v_late,2)) '  cells']})
% % 
% % %cdf
% % %early
% % figure(auc_fig_press_early);
% % subplot(nbins,3,2+((iAuc-1)*3))
% % hVvsInv = cdfplot(resp_v_early);
% % hVvsInv.Color = 'g';
% % hold on
% % hInvvsV = cdfplot(resp_a_early);
% % hInvvsV.Color = 'k';
% % [h p] = kstest2(resp_v_early,resp_a_early);
% % xlim([-0.05 0.05])
% % % title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
% % title({['ROC bin ' num2str(iAuc)];['early; p= ' num2str(p)]})
% % axis square
% % 
% % %late
% % figure(auc_fig_press_late);
% % subplot(nbins,3,2+((iAuc-1)*3))
% % hVvsInv = cdfplot(resp_v_late);
% % hVvsInv.Color = 'g';
% % hold on
% % hInvvsV = cdfplot(resp_a_late);
% % hInvvsV.Color = 'k';
% % [h p] = kstest2(resp_v_late,resp_a_late);
% % xlim([-0.05 0.05])
% % % title({['ROC bin ' num2str(iAuc) ', ' percentile_ind(iAuc) ':' percentile_ind(iAuc+1)];['valid vs invalid ' mouse(1).expt(1).cells(cellsInd).name ' pref; p= ' num2str(p)]})
% % title({['ROC bin ' num2str(iAuc)];['late; p= ' num2str(p)]})
% % axis square
% % 
% % %time-courses
% % %early
% % baseStimFrames_press = 0:cycTime:size(tc_vis_early,1)-1;
% % ttMs_press = (-pre_event_frames:minShortFrames-1)/(cycTime/cycTimeMs);
% % figure(auc_fig_press_early)
% % subplot(nbins,3,3+((iAuc-1)*3))
% % tc1 = nanmean(tc_vis_early,2)-nanmean(nanmean(tc_vis_early(pre_win,:),2));
% % err1 = nanstd(tc_vis_early,[],2)/sqrt(size(tc_vis_early(~isnan(tc_vis_early)),2));
% % shadedErrorBar(ttMs_press,tc1,err1,'g');
% % hold on
% % tc2 = nanmean(tc_aud_early,2)-nanmean(nanmean(tc_aud_early(pre_win,:),2));
% % err2 = nanstd(tc_aud_early,[],2)/sqrt(size(tc_aud_early(~isnan(tc_aud_early)),2));
% % shadedErrorBar(ttMs_press,tc2,err2,'k');
% % hold on
% % xlim([-15 minShortFrames-1]/(cycTime/cycTimeMs))
% % ylim([-0.01 0.05])
% % vline(baseStimFrames_press/(cycTime/cycTimeMs),':k')
% % hold on
% % vline([early_win(1)-pre_event_frames early_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
% % hold on
% % title({['ROC bin ' num2str(iAuc)]});
% % % axis square
% % 
% % %late
% % baseStimFrames_press = 0:cycTime:size(tc_vis_late,1)-1;
% % ttMs_press = (-pre_event_frames:minLongFrames-1)/(cycTime/cycTimeMs);
% % figure(auc_fig_press_late)
% % subplot(nbins,3,3+((iAuc-1)*3))
% % tc1 = nanmean(tc_vis_late,2)-nanmean(nanmean(tc_vis_late(pre_win,:),2));
% % err1 = nanstd(tc_vis_late,[],2)/sqrt(size(tc_vis_late(~isnan(tc_vis_late)),2));
% % shadedErrorBar(ttMs_press,tc1,err1,'g');
% % hold on
% % tc2 = nanmean(tc_aud_late,2)-nanmean(nanmean(tc_aud_late(pre_win,:),2));
% % err2 = nanstd(tc_aud_late,[],2)/sqrt(size(tc_aud_late(~isnan(tc_aud_late)),2));
% % shadedErrorBar(ttMs_press,tc2,err2,'k');
% % hold on
% % xlim([-15 minLongFrames-1]/(cycTime/cycTimeMs))
% % ylim([-0.01 0.05])
% % vline(baseStimFrames_press/(cycTime/cycTimeMs),':k')
% % hold on
% % vline([late_win(1)-pre_event_frames late_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
% % hold on
% % title({['ROC bin ' num2str(iAuc)]});
% % % axis square
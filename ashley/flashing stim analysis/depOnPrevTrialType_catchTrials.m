tc_hvsfa_all_prevTr = cell(1,3);
tc_favsh_all_prevTr = cell(1,3);
tc_mvscr_all_prevTr = cell(1,3);
tc_crvsm_all_prevTr = cell(1,3);
tc_hvscr_all_prevTr = cell(1,3);
tc_crvsh_all_prevTr = cell(1,3);
tc_hvsm_all_prevTr = cell(1,3);
tc_mvsh_all_prevTr = cell(1,3);
tc_favscr_all_prevTr = cell(1,3);
tc_crvsfa_all_prevTr = cell(1,3);

tc_val_all_prevTr = cell(1,3);
tc_inv_all_prevTr = cell(1,3);
tc_all_prevTr = cell(1,3);

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

hits = 1;
misses = 2;
fas = 3;
crs = 4;
vis_mat_name = cell(1,3);
a = 1;
for avInd = [1 3 4]
    vis_mat_name{a} = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).name;
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
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp,'unif',false)),3,length(cdirs));
            nDirs_h = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_h = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(misses).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(misses).stimResp,'unif',false)),3,length(cdirs));
            nDirs_m = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(misses).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_m = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(fas).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(fas).stimResp,'unif',false)),3,length(cdirs));
            nDirs_fa = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(fas).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
                nT = [ones(1,length(ind)) sz(3,:)];
                nDirs_fa = nT(ind3);
            end
            ckL = cell2mat(cellfun(@length,cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp,'unif',false),'unif',false));
            if all(ckL == 3)
            sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp,'unif',false)),3,length(cdirs));
            nDirs_cr = sz(3,:);
            else
                ind = find(ckL < 3);
                ind2 = find(ckL == 3);
                [j ind3] = sort([ind ind2]);
                sz = reshape(cell2mat(cellfun(@size,mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp(ind2),'unif',false)),3,length(cdirs)-length(ind));
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
            resp_h_dir_temp = cell(length(cdirs),1);
            resp_fa_dir_temp = cell(length(cdirs),1);
            if any(nDirs_hvsfa > 1)          
            resp_h = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            resp_fa = NaN(length(tt),length(cell_ind),sum(nDirs_hvsfa(nDirs_hvsfa > 1)),100);
            for iboot = 1:100
                resp_rTh = [];
                resp_rTfa = [];
                rH_dir_temp = cell(length(cdirs),1);
                rFA_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_hvsfa(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rH_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_h_dir_temp{idir} = cat(4,resp_h_dir_temp{idir},rH_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir))));
                    rFA_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsfa(idir)));
                    resp_fa_dir_temp{idir} = cat(4,resp_fa_dir_temp{idir},rFA_dir_temp{idir});
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_fa(:,:,:,iboot) = resp_rTfa;
            end
            resp_h = nanmean(resp_h,4);
            resp_fa = nanmean(resp_fa,4);
            
            
            else
            resp_h = [];
            resp_fa = [];
            end
            
            
            tc_hvsfa_all_prevTr{a} = cat(2,tc_hvsfa_all_prevTr{a},mean(resp_h,3));
            tc_favsh_all_prevTr{a} = cat(2,tc_favsh_all_prevTr{a},mean(resp_fa,3));
            
            if any(nDirs_hvsfa > 1)
            val_temp = resp_h;
            inv_temp = resp_fa;
            all_temp = cat(3,resp_h,resp_fa);
            else
                val_temp = [];
                inv_temp = [];
                all_temp = [];
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
            for iboot = 1:100
                resp_rTm = [];
                resp_rTcr = [];
                rM_dir_temp = cell(length(cdirs),1);
                rCR_dir_temp = cell(length(cdirs),1);
            for idir = 1:length(cdirs)
                if nDirs_mvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rM_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_m_dir_temp{idir} = cat(4,resp_m_dir_temp{idir},rM_dir_temp{idir});
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir))));
                    rCR_dir_temp{idir} = rT(:,cell_ind,randperm(size(rT,3),nDirs_mvscr(idir)));
                    resp_cr_dir_temp{idir} = cat(4,resp_cr_dir_temp{idir},rCR_dir_temp{idir});
                end
            end
            end
            resp_m(:,:,:,iboot) = resp_rTm;
            resp_cr(:,:,:,iboot) = resp_rTcr;            
            
            resp_m = nanmean(resp_m,4);
            resp_cr = nanmean(resp_cr,4);
            else
                resp_m = [];
                resp_cr = [];
            end
            
            tc_mvscr_all_prevTr{a} = cat(2,tc_mvscr_all_prevTr{a},mean(resp_m,3)); 
            tc_crvsm_all_prevTr{a} = cat(2,tc_crvsm_all_prevTr{a},mean(resp_cr,3)); 
            
            if any(nDirs_mvscr > 1)
                val_temp = cat(3,val_temp,resp_m);
                inv_temp = cat(3,inv_temp,resp_cr);
                all_temp = cat(3,all_temp,resp_m,resp_cr);
            else
                val_temp = val_temp;
                inv_temp = inv_temp;
                all_temp = all_temp;
            end
            
            tc_val_all_prevTr{a} = cat(2,tc_val_all_prevTr{a},mean(val_temp,3));%mean(resp_m,3));
            tc_inv_all_prevTr{a} = cat(2,tc_inv_all_prevTr{a},mean(inv_temp,3));%mean(resp_cr,3));
            tc_all_prevTr{a} = cat(2,tc_all_prevTr{a},mean(all_temp,3));


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
            for iboot = 1:100
            resp_rTh = [];
            resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_hvscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvscr(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_cr(:,:,:,iboot) = resp_rTcr;
            end
            resp_h = mean(resp_h,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_h = [];
            resp_cr = [];
            end   
            
            tc_hvscr_all_prevTr{a} = cat(2,tc_hvscr_all_prevTr{a},mean(resp_h,3)); 
            tc_crvsh_all_prevTr{a} = cat(2,tc_crvsh_all_prevTr{a},mean(resp_cr,3));  

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
            for iboot = 1:100
                resp_rTh = [];
                resp_rTm = [];
            for idir = 1:length(cdirs)
                if nDirs_hvsm(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(hits).stimResp{idir};
                    resp_rTh = cat(3,resp_rTh,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(misses).stimResp{idir};
                    resp_rTm = cat(3,resp_rTm,rT(:,cell_ind,randperm(size(rT,3),nDirs_hvsm(idir))));
                end
            end
            resp_h(:,:,:,iboot) = resp_rTh;
            resp_m(:,:,:,iboot) = resp_rTm;
            end
            resp_h = mean(resp_h,4);
            resp_m = mean(resp_m,4);
            else
            resp_h = [];
            resp_m = [];
            end
            
            
            tc_hvsm_all_prevTr{a} = cat(2,tc_hvsm_all_prevTr{a},mean(resp_h,3)); 
            tc_mvsh_all_prevTr{a} = cat(2,tc_mvsh_all_prevTr{a},mean(resp_m,3)); 

          
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
            for iboot = 1:100
                resp_rTfa = [];
                resp_rTcr = [];
            for idir = 1:length(cdirs)
                if nDirs_favscr(idir) > 1
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(fas).stimResp{idir};
                    resp_rTfa = cat(3,resp_rTfa,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                    rT = mouse(imouse).expt(iexp).align(catchAlign).av(avInd).outcome(crs).stimResp{idir};
                    resp_rTcr = cat(3,resp_rTcr,rT(:,cell_ind,randperm(size(rT,3),nDirs_favscr(idir))));
                end
            end
            resp_fa(:,:,:,iboot) = resp_rTfa;
            resp_cr(:,:,:,iboot) = resp_rTcr;
            end
            resp_fa = mean(resp_fa,4);
            resp_cr = mean(resp_cr,4);
            else
            resp_fa = [];
            resp_cr = [];
            end
            
            tc_favscr_all_prevTr{a} = cat(2,tc_favscr_all_prevTr{a},mean(resp_fa,3)); 
            tc_crvsfa_all_prevTr{a} = cat(2,tc_crvsfa_all_prevTr{a},mean(resp_cr,3)); 
            
        end
    end
end
    a = a+1;
end
%% summary stats
tc_hvsfa_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_hvsfa_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_hvsfa_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_favsh_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_favsh_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_favsh_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_mvscr_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_mvscr_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_mvscr_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_crvsm_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_crvsm_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_crvsm_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_hvsm_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_hvsm_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_hvsm_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_mvsh_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_mvsh_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_mvsh_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_favscr_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_favscr_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_favscr_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_crvsfa_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_crvsfa_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_crvsfa_all_prevTr,'unif',false),'unif',false),'unif',false);

tc_val_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_val_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_val_all_prevTr,'unif',false),'unif',false),'unif',false);
tc_inv_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x,2),tc_inv_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),cellfun(@(x) mean(x,2),tc_inv_all_prevTr,'unif',false),'unif',false),'unif',false);

resp_hvsfa_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_hvsfa_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_hvsfa_all_prevTr,'unif',false),'unif',false);
resp_favsh_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_favsh_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_favsh_all_prevTr,'unif',false),'unif',false);
resp_mvscr_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_mvscr_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_mvscr_all_prevTr,'unif',false),'unif',false);
resp_crvsm_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_crvsm_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_crvsm_all_prevTr,'unif',false),'unif',false);
resp_hvsm_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_hvsm_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_hvsm_all_prevTr,'unif',false),'unif',false);
resp_mvsh_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_mvsh_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_mvsh_all_prevTr,'unif',false),'unif',false);
resp_favscr_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_favscr_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_favscr_all_prevTr,'unif',false),'unif',false);
resp_crvsfa_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_crvsfa_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_crvsfa_all_prevTr,'unif',false),'unif',false);

resp_val_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_val_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_val_all_prevTr,'unif',false),'unif',false);
resp_inv_prevTr_mean = cellfun(@minus, cellfun(@(x) mean(x(trans_win,:),1),tc_inv_all_prevTr,'unif',false) , cellfun(@(y) mean(y(pre_win,:),1),tc_inv_all_prevTr,'unif',false),'unif',false);


tc_hvsfa_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_hvsfa_all_prevTr,'unif',false);
tc_favsh_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_favsh_all_prevTr,'unif',false);
tc_mvscr_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_mvscr_all_prevTr,'unif',false);
tc_crvsm_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_crvsm_all_prevTr,'unif',false);
tc_hvsm_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_hvsm_all_prevTr,'unif',false);
tc_mvsh_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_mvsh_all_prevTr,'unif',false);
tc_favscr_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_favscr_all_prevTr,'unif',false);
tc_crvsfa_prevTr_ste =  cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_crvsfa_all_prevTr,'unif',false);

tc_val_prevTr_ste = cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_val_all_prevTr,'unif',false);
tc_inv_prevTr_ste = cellfun(@(x) std(x,[],2)./sqrt(size(x,2)),tc_inv_all_prevTr,'unif',false);
%% analysis strings

varStr = {'_val';'_inv';'_hvsfa';'_favsh';'_mvscr';'_crvsm';'_hvsm';'_mvsh';'_favscr';'_crvsfa'};
colStr = {'k';'c';'k';'c';'r';'b';'k';'r';'c';'b'};
trTypeStr = {'all val';'all inv';'hits val';'hits inv';'misses val';'misses inv';'hits val';'misses val';'hits inv';'misses inv'};

%% plotting - timecourses
% all valid vs all invalid trials
plotMat = cat(1,1:2:10,2:2:10);

for ifig = 1:5;
tr1 = plotMat(1,ifig);
tr2 = plotMat(2,ifig);

prevTrialAnalysisFig = figure;
% suptitle({[titleStr ' - resp matched across target types']; 'black-all val,cyan-all inv'})
suptitle({[titleStr ' - resp matched across target types']; [colStr{tr1} '-' trTypeStr{tr1} ',' colStr{tr2} '-' trTypeStr{tr2}]})

% all val vs inv
subplot(4,2,1)
l1_ind = 1;
l2_ind = 1;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l1_ind},'Location','northwest')

% vis-prev vs aud-prev
subplot(4,2,3)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l2_ind}']),colStr{tr1});
if strcmp(colStr{tr1},'k')
l2.mainLine.Color = l2.mainLine.Color+0.5;
l2.patch.FaceColor = l2.patch.FaceColor./2;
else
l2.mainLine.Color = l2.mainLine.Color./2;
l2.patch.FaceColor = l2.mainLine.Color+0.25;
end
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
% vis-prev vs aud-prev
subplot(4,2,4)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l1_ind}']),colStr{tr2});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
l2.mainLine.Color = l2.mainLine.Color./2;
l2.patch.FaceColor = l2.mainLine.Color+0.25;
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')

% vis-prev vs vis-prev
subplot(4,2,5)
l1_ind = 2;
l2_ind = 2;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
% aud-prev vs aud-prev
subplot(4,2,6)
l1_ind = 3;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')

% vis-prev vsaud-prev
subplot(4,2,7)
l1_ind = 2;
l2_ind = 3;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
% aud-prev vs vis-prev
subplot(4,2,8)
l1_ind = 3;
l2_ind = 2;
l1 = shadedErrorBar(ttMs,eval(['tc' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['tc' varStr{tr1} '_prevTr_ste{l1_ind}']),colStr{tr1});
hold on
l2 = shadedErrorBar(ttMs,eval(['tc' varStr{tr2} '_prevTr_mean{l2_ind}']),eval(['tc' varStr{tr2} '_prevTr_ste{l2_ind}']),colStr{tr2});
hold on
legend([l1.mainLine l2.mainLine],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')

figure(prevTrialAnalysisFig)
for iplot = 1:8;
    subplot(4,2,iplot)
    hold on
    xlim([-200 1000])
    ylim([-0.01 0.05])
    vline(baseStimFrames/(cycTime/cycTimeMs),':k')
    vline([trans_win(1)-pre_event_frames trans_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--r')
    hold on
    vline([pre_win(1)-pre_event_frames pre_win(end)-pre_event_frames]/(cycTime/cycTimeMs), '--k')
    xlabel('t(ms)')
    ylabel('dF/F')
end

print([fnout 'target_align_respByPrevTrial' varStr{tr1} varStr{tr2} datasetStr '.pdf'], '-dpdf');

end

%% transient response to target and kstest across population

for ifig = 1:5;
tr1 = plotMat(1,ifig);
tr2 = plotMat(2,ifig);

prevTrialAnalysisCDFFig = figure;
% suptitle({[titleStr ' - resp matched across target types']; 'black-all val,cyan-all inv'})
suptitle({[titleStr ' - resp matched, kstest']; [colStr{tr1} '-' trTypeStr{tr1} ',' colStr{tr2} '-' trTypeStr{tr2}]})

% all val vs inv
subplot(4,2,1)
l1_ind = 1;
l2_ind = 1;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l1_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])

% vis-prev vs aud-prev
subplot(4,2,3)
l1_ind = 2;
l2_ind = 3;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr1};
if strcmp(colStr{tr1},'k')
l2.Color = l2.Color+0.5;
else
l2.Color = l2.Color./2;
end
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr1} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])
% vis-prev vs aud-prev
subplot(4,2,4)
l1_ind = 2;
l2_ind = 3;
l1 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr2};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
l2.Color = l2.Color./2;
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr2} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])

% vis-prev vs vis-prev
subplot(4,2,5)
l1_ind = 2;
l2_ind = 2;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])
% aud-prev vs aud-prev
subplot(4,2,6)
l1_ind = 3;
l2_ind = 3;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])

% vis-prev vs aud-prev
subplot(4,2,7)
l1_ind = 2;
l2_ind = 3;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])
% aud-prev vs vis-prev
subplot(4,2,8)
l1_ind = 3;
l2_ind = 2;
l1 = cdfplot(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']));
l1.Color = colStr{tr1};
hold on
l2 = cdfplot(eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
l2.Color = colStr{tr2};
hold on
legend([l1 l2],vis_mat_name{l1_ind},vis_mat_name{l2_ind},'Location','northwest')
[h,p] = kstest2(eval(['resp' varStr{tr1} '_prevTr_mean{l1_ind}']),eval(['resp' varStr{tr2} '_prevTr_mean{l2_ind}']));
title(['p = ' num2str(p)])

figure(prevTrialAnalysisCDFFig)
for iplot = 1:8;
    subplot(4,2,iplot)
    hold on
%     xlim([-200 1000])
%     ylim([-0.01 0.05])
    xlabel('dF/F')
    ylabel('cmlv frctn')
end

print([fnout 'target_align_CDFByPrevTrial' varStr{tr1} varStr{tr2} datasetStr '.pdf'], '-dpdf');
end
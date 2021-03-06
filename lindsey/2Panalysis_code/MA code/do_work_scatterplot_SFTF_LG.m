function RegInfo = do_work_scatterplot_SFTF(RegInfo);
%goal is to start with RegInfo, and get SF and TF and Xi estimates for all
%cells and confidence intervals..

%first, for each cell, get timecourses, sort by stimtypes and by running
%(and later, by eyeposition and by whether the 3D correg was good)

if isfield(RegInfo,'Stats') == 0
    RegInfo.Stats = [];
end
if isfield(RegInfo.Stats,'Fit_struct') == 0
    RegInfo.Stats.Fit_struct = [];
end
if isfield(RegInfo.Stats.Fit_struct,'Shuf') == 0
    RegInfo.Stats.Fit_struct.Shuf = 1;
    RegInfo.Stats.Fit_struct.True = 1;
end


%plan: make a structure AllStat_IndivMA that includes all relevant info for
%each individual expti
EXPT = RegInfo.All.EXPT;
WIN_PRE_USE = RegInfo.Runs(1).WIN_PRE;
WIN_POST_USE = RegInfo.Runs(1).WIN_POST;

%WIN_PRE_USE = [-2 0];
%WIN_POST_USE = [0 5];
%RegInfo.All.WIN_PRE_USE = WIN_PRE_USE;
%RegInfo.All.WIN_POST_USE = WIN_POST_USE;

tc_USE_str = 'tc_VBnpsub';
Ntot1_USE =  RegInfo.Runs(1).Ntot1;
t_USE = RegInfo.Runs(1).t(1:Ntot1_USE);

AllStat_IndivMA = [];

%first, for each cell, get timecourses, sort by stimtypes and by running
%(and later, by eyeposition and by whether the 3D correg was good)
tc_USE = [];
Nruns = length(RegInfo.Runs);
for count = 1:Nruns
   %add neuropil timecourse following all cell timecourses..
    tmp =  RegInfo.Runs(count).nptc_VBnpsub;

    %    tc_USE = [tc_USE; eval(['RegInfo.Runs(count).',tc_USE_str])];
    tc_USE = [tc_USE; [eval(['RegInfo.Runs(count).',tc_USE_str]) tmp(:,1)]];
%    tc_USE = [tc_USE; eval(['RegInfo.Runs(count).',tc_USE_str])];
end


%%% Add additional computations for computing cell max Xsectional area and
%%% height, as well as cross-correlation
Ncells6 = RegInfo.Cell.nCells;
Area0 = RegInfo.Cell.area;
Stats_cellarea = zeros(Ncells6,4); 
tmp00 = RegInfo.Cell.labelimg;
scaleXY = RegInfo.All.FOV(1)./128;
XC_cells = zeros(Ncells6,Ncells6);
Is_adjacent = zeros(Ncells6,Ncells6);
%choose prestim times to do the correlation: 
t1 = [-4:5]'*ones(1,size(tc_USE,1)/10);
t = reshape(t1,size(t1,1)*size(t1,2),1);
ind_pre00 = find(t>=-2 & t<0);
for count6 = 1:Ncells6
   tmp0 = tmp00 == count6;
   tmp1 = squeeze(sum(sum(tmp0,1),2));
   XSection_max = max(tmp1);
   tmp2 = find(tmp1);
   Height = max(tmp2) - min(tmp2) + 1;
   Stats_cellarea(count6,:) = [XSection_max Height Area0(count6) round(10*sqrt(XSection_max*scaleXY/pi)*2)./10];
end

Mask_crossection = ones(Ncells6+1,1);
thresh_crosssection = 7; %um
Mask_crossection(find(Stats_cellarea(:,4) <  thresh_crosssection)) = 0;

 OLD = 1;
lowcut1 = 50;
 if OLD == 1
     %Now compute XC of each cell and each other cells
    tc1 = tc_USE(:,count6);
   for count6b = 1:(count6-1)
       tc2 = tc_USE(:,count6b);

       
       %       tmp3 = corrcoef(tc1(ind_pre00),tc2(ind_pre00));
       tmp3 = corrcoef(low_cutKO(tc1(ind_pre00),lowcut1),low_cutKO(tc2(ind_pre00),lowcut1));
       XC_cells(count6,count6b) = tmp3(2,1);       
   end
 end
 
for count6 = 1:Ncells6
    tmp5a = double(RegInfo.Cell.labelimg == count6);
   for count6b = 1:(count6-1)
       tmp5b = double(RegInfo.Cell.labelimg == count6b); 
       
       tmp5c = tmp5a + tmp5b;
       tmp6 = bwlabeln(tmp5c,6);
       if max(max(max(tmp6)))<2
           Is_adjacent(count6,count6b) = 1;
       end
   end
end

Is_adjacent = diag(Mask_crossection(1:Ncells6))*Is_adjacent;
Is_adjacent = Is_adjacent + Is_adjacent';


%remove chunks with smaller xsections
Mask_nosmallbits = ones(Ncells6+1,1);
Is_adjacent2  = Is_adjacent;
for count6 = 1:Ncells6
    ind = [count6; find(Is_adjacent2(:,count6))];
    ind2 = find( Stats_cellarea(ind,1) < max( Stats_cellarea(ind,1) ));
    if length(ind2)>0
        Mask_nosmallbits(ind(ind2)) = 0;
        Is_adjacent2(ind(ind2),:) = 0;
        Is_adjacent2(:,ind(ind2)) = 0;
    end
end
Stats_cellarea = [Stats_cellarea Mask_nosmallbits(1:Ncells6) Mask_crossection(1:Ncells6)];

RegInfo.Stats.Stats_cellarea = Stats_cellarea;
RegInfo.Stats.thresh_crosssection = thresh_crosssection;
RegInfo.Stats.XC_cells = XC_cells;
RegInfo.Stats.lowcut1 = lowcut1;
RegInfo.Stats.labelimg = RegInfo.Cell.labelimg;
RegInfo.Stats.avg_norm = RegInfo.Cell.avg_norm;
%%% done with additional calcs of cell size and correlation, etc. 



Npts = size(tc_USE,1);
Ncells = size(tc_USE,2);
Ntrials = Npts/Ntot1_USE;

%tc_USE2 = reshape(tc_USE,Ntot1_USE,Ntrials,Ncells);
%make tc_USE2 have pre and post
%add an extra 10s onto tc_USE for wrap around: 
tc_USEB = [tc_USE; tc_USE(end-9:end,:)];
ind_TC = -5:10;
Trigs_TC = 6:10:(length(tc_USE));
Ntrials_TC = length(Trigs_TC);
tc_USE2 = zeros(length(ind_TC),Ntrials,Ncells); %%LG look here
for count_TC = 1:Ntrials
    tc_USE2(:,count_TC,:) = tc_USEB(Trigs_TC(count_TC) + ind_TC,:);
end

%tc_USE2_mean = squeeze(mean(tc_USE2,2));
ind_pre = find(t_USE > WIN_PRE_USE(1) & t_USE<=WIN_PRE_USE(2));
ind_post = find(t_USE > WIN_POST_USE(1) & t_USE<=WIN_POST_USE(2));
tc_resp_post = squeeze(mean(tc_USE2(ind_post,:,:),1));
tc_resp_pre = squeeze(mean(tc_USE2(ind_pre,:,:),1));
tc_resp = (tc_resp_post - tc_resp_pre)./tc_resp_pre;

tc_resp_pre2 = zeros(size(tc_USE2));
for count1 = 1:size(tc_USE2,1)
    tc_resp_pre2(count1,:,:) = tc_resp_pre;
end
tc_resp2 = (tc_USE2 - tc_resp_pre2)./tc_resp_pre2;


%now get Wheel info:
Wheel_mat_USE = RegInfo.State.Wheel_mat > 0; %for now, any motion is considered motion..
Wheel_mat_USE2 = reshape(Wheel_mat_USE,Ntot1_USE,Ntrials)';
%Wheel_mat_USE3 = sum(Wheel_mat_USE2(:,[ind_pre ind_post]),2)>0; %for now, any running during response window or prestim window counts..
Wheel_mat_USE3a = sum(Wheel_mat_USE2(:,[ind_post]),2); %for now, any running during response wwindow counts..
Wheel_mat_USE3 = sum(Wheel_mat_USE2(:,[ind_post]),2)>0; %for now, any running during response wwindow counts..

Mask_trials_USE = Wheel_mat_USE3; %can add eyetracking here, EEG, etc..


%%%%% figure out stim randomization

RandStim_info = RegInfo.Rand.RandStim_info;
count_protocol = 1; %count_prot_USE;
%plane = PLANE;
%str_prot = ['prot',num2str(count_protocol),'_plane',num2str(plane)];

userun = RandStim_info.prot(count_protocol).userun0
TFvec = RandStim_info.prot(count_protocol).TFvec;
SFvec = RandStim_info.prot(count_protocol).SFvec;
pos = RandStim_info.prot(count_protocol).pos;
StimXpos = RandStim_info.prot(count_protocol).StimXpos;
StimYpos = RandStim_info.prot(count_protocol).StimYpos;
Nstimtypes0000 = RandStim_info.prot(count_protocol).Nstimtypes0000;
oris = RandStim_info.prot(count_protocol).oris;
contrastvec = RandStim_info.prot(count_protocol).contrastvec;
userun0 = RandStim_info.prot(count_protocol).userun0;
userunDR = RandStim_info.prot(count_protocol).userunDR;
Fs1 = RandStim_info.prot(count_protocol).Fs1;
Fs1_orig = RandStim_info.prot(count_protocol).Fs1_orig;
Toff1 = RandStim_info.prot(count_protocol).Toff1;
Ton1 = RandStim_info.prot(count_protocol).Ton1;
Noff1 = RandStim_info.prot(count_protocol).Noff1;
Non1 = RandStim_info.prot(count_protocol).Non1;
Ntot1 = RandStim_info.prot(count_protocol).Ntot1;
DIR_rand0 = RandStim_info.general.DIR_rand0;
file_mat = RandStim_info.general.file_mat;
randvec = RandStim_info.general.randvec;
ind_start = Non1;
ind_end = Non1; %%note: using different end times, so choosing the longest one for now..
ind_start_vec = ones(1,Nstimtypes0000)*ind_start;
ind_end_vec = ones(1,Nstimtypes0000)*ind_end;

rate = 1;


%    load_exptinfo_Aug2010

Seq2 = [];
for count_run10 = 1:length(userun)
    OLD = 0;
    if OLD == 1
        tmp = randvec(userun(count_run10),2);
        tmp = rem(tmp,8);
        tmp(find(tmp==0)) = 8;
        tmp_DIR = deblank(DIR_rand0(count_prot_USE,:));
        DIR_rand = dir([tmp_DIR,['*_',num2str(tmp),'__','*.mat']]);
        if length(DIR_rand) == 1
            file_rand = [tmp_DIR,DIR_rand(1).name];
        else
            error('error: wrong # of matches');
            return
        end
        load([file_rand]);
    else
        %load Seq from RegInfo
        Seq = RegInfo.Runs(count_run10).Seq;
    end

    Npts = length(RegInfo.Runs(count_run10).tc_VBnpsub);
    Nframesperstim0 = Ntot1_USE;
%    Ntrials_USE = min(Npts./Nframesperstim0, size(Seq,1));
    Ntrials_USE = Npts./Nframesperstim0;
    if size(Seq,1) < Ntrials_USE; %for 110516 and 110517, each run had 10 extra blank trials at end
        Npad = Ntrials_USE - size(Seq,1);
        Seq = [Seq; -1*ones(Npad,size(Seq,2))];
    end
    

    Seq = Seq(1:Ntrials_USE,:);
    if size(Seq,2) == 7 %a quick fix for 110205 ONLY
        Seq = [Seq, ones(size(Seq,1),1)*36 ones(size(Seq,1),1)*-18];
    end
    %Seq2 = [Seq2; Seq];
    Seq2 = [Seq2; Seq(:,:)];
end


Seqposition = [];
Seqposition2 = [];

%        Nstim1 = Nframes1./(Ntot1*Fs1);
%for now, only use if each run is the same ********************* change
Nstim1 = min(size(Seq2,1),Ntrials);

Nstimtypes1 = length(TFvec);
Nstimtypes = Nstimtypes1;
NstimtypesMA = Nstimtypes1;

Nindvec = [];
tmp_mask_vec = zeros(Nstim1,Nstimtypes1);
Stimvalue = zeros(Nstim1,1);
Stimvalue2 = zeros(Nstim1,Nstimtypes);
for count_seq = 1:Nstimtypes1
    tmp_mask = ...
        Seq2(1:Nstim1,4) == TFvec(count_seq) & ...
        Seq2(1:Nstim1,5) == SFvec(count_seq) & ...
        Seq2(1:Nstim1,6) == oris(count_seq) & ...
        Seq2(1:Nstim1,7) == contrastvec(count_seq) & ...
        Seq2(1:Nstim1,8) == StimXpos(count_seq) & ...
        Seq2(1:Nstim1,9) == StimYpos(count_seq);
    ind = find(tmp_mask);
    Stimvalue(ind) = count_seq;
    Stimvalue2(ind,count_seq) = count_seq;
    tmp_mask_vec(:,count_seq) = tmp_mask;
    Seqposition(count_seq).ind = ind;
    Nindvec = [Nindvec; length(ind)];
    Stimvec =  [TFvec(count_seq) SFvec(count_seq) oris(count_seq) contrastvec(count_seq) ...
        StimXpos(count_seq) StimYpos(count_seq)];
    Seqposition(count_seq).TFSFetc = Stimvec;
    Seqposition2 = [Seqposition2; [ones(length(ind),1)*[count_seq Stimvec] ind]];
end


%%%%%
%Organize SF, TF etc:

tmp = sort(TFvec(1:(Nstimtypes-1)));
TFvec000 = tmp([1 find(diff(tmp))+1]);
tmp = sort(SFvec(1:(Nstimtypes-1)));
SFvec000 = tmp([1 find(diff(tmp))+1]);

TFuse = [.5 1 2 4 8 15 24];
SFuse = [.02 .04 .08 .16 .32];
ind_TFuse = [];
for count = 1:length(TFuse)
    ind2 = find(TFuse(count) == TFvec000);
    if length(ind2)>0
        ind_TFuse = [ind_TFuse; [ind2(1) TFuse(count)]];
    end
end
NTF = length(ind_TFuse);

ind_SFuse = [];
for count = 1:length(SFuse)
    ind2 = find(SFuse(count) == SFvec000);
    if length(ind2)>0
        ind_SFuse = [ind_SFuse; [ind2(1) SFuse(count)]];
    end
end
ind_SFuse = flipud(ind_SFuse);
NSF = length(ind_SFuse);
%%%end stim randomization work

PLOTRUN_NORUN_RUN_ALLRUN_vec = [1:3] %run = 1, norun= 2, allrun = 3
RUN_EO_vec = [3]; %evens vs odds, set to 3 to take all trials

for RUN_EO = RUN_EO_vec;
    Mask_stimtype = zeros(length(Wheel_mat_USE3),Nstimtypes);
    for count_stim = 1:NstimtypesMA
        Mask_stim = Stimvalue2(:,count_stim) == count_stim;
        ind_Mask_stim = find(Mask_stim);
        if RUN_EO == 1;
            Mask_stimtype(ind_Mask_stim(1:2:end),count_stim) = count_stim;
        elseif RUN_EO == 2;
            Mask_stimtype(ind_Mask_stim(2:2:end),count_stim) = count_stim;
        elseif RUN_EO == 3;
            Mask_stimtype(ind_Mask_stim(1:1:end),count_stim) = count_stim;
        end
    end

    REJECT_LARGERESP = 0;
    if REJECT_LARGERESP == 1
        %choose an ROI and reject huge negative responses
        tmp = squeeze(mean(tc_resp,2));
        tmp_Mean = mean(tmp);
        tmp_std = std(tmp);
        Nstd = 2;

        %    NoRej_mat = ((tmp > (tmp_Mean - Nstd*tmp_std)) & (tmp < (tmp_Mean + Nstd*tmp_std)));
        NoRej_mat = ( (tmp > (tmp_Mean - Nstd*tmp_std)));
    else
        NoRej_mat = ones(Ntrials,1);
    end

    Ind_struct = [];

    for count = 1:NstimtypesMA
        Ind_struct(count).run_trials = find(Mask_stimtype(:,count) == count & Wheel_mat_USE3 == 1 & NoRej_mat);
        Ind_struct(count).norun_trials = find(Mask_stimtype(:,count) == count & Wheel_mat_USE3 == 0 & NoRej_mat);
        Ind_struct(count).allrun_trials = find(Mask_stimtype(:,count) == count & NoRej_mat);
    end


    %%first, for allrun, find cells where at least 1 stimtype is
    %%robustly driven (~p<.001)
    str_run = 'allrun';


    %now compute means:
    eval(['Im_N_',str_run,' = zeros(NstimtypesMA,1);']);
    eval(['Im_mat_',str_run,' = zeros(Ncells,NstimtypesMA);']);
    eval(['Im_mat_std_',str_run,' = zeros(Ncells,NstimtypesMA);']);
    for count = 1:NstimtypesMA
        ind_run_1 = eval(['Ind_struct(count).',str_run,'_trials;']);
        eval(['Im_N_',str_run,'(count) = length(ind_run_1);']);
        if length(ind_run_1)>0
            eval(['Im_mat_',str_run,'(:,count) = mean(tc_resp(ind_run_1,:),1);']);
            eval(['Im_mat_std_',str_run,'(:,count) = std(tc_resp(ind_run_1,:),[],1);']);
        else
            %make all stimtypes where there is no data equal to NaN
            eval(['Im_mat_',str_run,'(:,count) = NaN;']);
            eval(['Im_mat_std_',str_run,'(:,count) = NaN;']);        
        end
    end

    
    eval(['Im_mat_USE = Im_mat_',str_run,';']);
    eval(['Im_mat_std_USE = Im_mat_std_',str_run,';']);
    eval(['Im_N_USE = Im_N_',str_run,';']);
    Im_mat_USE_true = Im_mat_USE;
    Im_mat_std_USE_true = Im_mat_std_USE;
    Im_N_USE_true = Im_N_USE;


    %now generate the mean for fitting:

    %first, root out crappy cells as ones without at least one
    %significant, at significance level alpha2
   
   TF_vec0 = ind_TFuse(:,2);
    SF_vec0 = ind_SFuse(:,2);
    [sfsf,tftf]=meshgrid(SF_vec0,TF_vec0);
    grid2.sfsf = sfsf;
    grid2.tftf = tftf;
%for estimating highcutoff sample more precisely, and estimate 50% and 10%
%high cutoff.. 
dSF = median(diff(log2(SF_vec0)));
dTF = median(diff(log2(TF_vec0)));
SF_vec00 = log2(SF_vec0(1)):(dSF/10):log2(SF_vec0(end));
TF_vec00 = log2(TF_vec0(1)):(dTF/10):log2(TF_vec0(end));
[sfsf00,tftf00]=meshgrid(SF_vec00,TF_vec00);
grid2.sfsf00 = sfsf00;
grid2.tftf00 = tftf00;

    
    Fit_struct = [];

    %repeat this for true data and for Nshuf shuffles
    Nshuf = 500; %100; %=0 means only run for true data
    Mask_cells_USE_ALL2 = zeros(Ncells,3);
    for count_shuf = 0:Nshuf
        for count_run = 1:length(PLOTRUN_NORUN_RUN_ALLRUN_vec)
            PLOTRUN_NORUN_RUN_ALLRUN = PLOTRUN_NORUN_RUN_ALLRUN_vec(count_run);
            if PLOTRUN_NORUN_RUN_ALLRUN ==1
                str_run = 'run';
            elseif PLOTRUN_NORUN_RUN_ALLRUN ==2
                str_run = 'norun';
            else
                str_run = 'allrun';
            end

            display([EXPT,'; ',str_run,': count_shuf = ',num2str(count_shuf)]);

            %now compute means:
            eval(['Im_N_',str_run,' = zeros(Nstimtypes,1);']);
            eval(['Im_mat_',str_run,' = zeros(Ncells,NstimtypesMA);']);
            eval(['Im_mat_std_',str_run,' = zeros(Ncells,NstimtypesMA);']);
            for count = 1:NstimtypesMA
                ind_run_1 = eval(['Ind_struct(count).',str_run,'_trials;']);
                if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                    ind_run_1 = ind_run_1(randsample(length(ind_run_1),length(ind_run_1),1));
                end
                    
                eval(['Im_N_',str_run,'(count) = length(ind_run_1);']);
                if length(ind_run_1)>0
                    eval(['Im_mat_',str_run,'(:,count) = mean(tc_resp(ind_run_1,:),1);']);
                    eval(['Im_mat_std_',str_run,'(:,count) = std(tc_resp(ind_run_1,:),[],1);']);
                else
                    %make all stimtypes where there is no data equal to NaN
                    eval(['Im_mat_',str_run,'(:,count) = NaN;']);
                    eval(['Im_mat_std_',str_run,'(:,count) = NaN;']);
                end
            end

            eval(['Im_mat_USE = Im_mat_',str_run,';']);
            eval(['Im_mat_std_USE = Im_mat_std_',str_run,';']);
            eval(['Im_N_USE = Im_N_',str_run,';']);

            
            if count_shuf == 0
                %first, root out crappy cells as ones without at least one
                %significant, at significance level alpha2
%%            
               Info_ttest_mat = ones(Ncells,NstimtypesMA,3);
                for count = 1:NstimtypesMA
                    ind_run_1 = eval(['Ind_struct(count).',str_run,'_trials;']);
                    if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                        ind_run_1 = ind_run_1(randsample(length(ind_run_1),length(ind_run_1),1));
                    end

                    %eval(['Im_N_',str_run,'(count) = length(ind_run_1);']);
                    if length(ind_run_1)>0
                        Rpost = tc_resp_post(ind_run_1,:);
                        Rpre = tc_resp_pre(ind_run_1,:);
                        
                        tmp2 = (Rpost-Rpre)./Rpre;
                        tmp2mean = mean(tmp2,1);
                        tmp2std = std(tmp2,[],1);
                        N = size(tmp2,1);
                        tstat1 = tmp2mean./(tmp2std./sqrt(N));
                        p_ttestA = 1 - tcdf(tstat1,N-1);
                        %                H_ttest = min(p_ttest,[],2) < alpha2;
                        %[h,p] =
                        
                        p_ttestB = zeros(1,Ncells);
                        alphaB = .05./(Nstimtypes-1);
                        for cellnum = 1:Ncells;
%                            [h_ttestB1,p_ttestB1] = ttest2(Rpre(:,cellnum),Rpost(:,cellnum),alphaB,'left','equal');
                            [h_ttestB1,p_ttestB1] = ttest2(tc_resp_pre(:,cellnum),Rpost(:,cellnum),alphaB,'left','unequal');
                            p_ttestB(cellnum) = p_ttestB1;
                        end
                        Info_ttest_mat(:,count,1) = p_ttestA;
                        Info_ttest_mat(:,count,2) = p_ttestB;
                        %eval(['Im_mat_',str_run,'(:,count) = mean(tc_resp(ind_run_1,:),1);']);
                        %eval(['Im_mat_std_',str_run,'(:,count) = std(tc_resp(ind_run_1,:),[],1);']);
                        %hi = 1
                        %pause(10)
                        
                    end
                end                
                
                alpha2 = .05./(NstimtypesMA-1); %Bonferroni correction, doesn't require indep; don't include blank in deg freedom
                
                p_ttest = squeeze(Info_ttest_mat(:,:,2));
                H_ttest = min(p_ttest(:,1:(NstimtypesMA-1)),[],2) < alpha2;
                %Im_N_USE2 = (ones(Ncells,1)*Im_N_USE');
                %tstat = Im_mat_USE./(Im_mat_std_USE./sqrt(Im_N_USE2));
                %p_ttest = 1 - tcdf(tstat,Im_N_USE2-1);
                %H_ttest = min(p_ttest,[],2) < alpha2;

%                alpha2 = .01; %.05./(NstimtypesMA-1); %Bonferroni correction, doesn't require indep; don't include blank in deg freedom

%                Im_N_USE2 = (ones(Ncells,1)*Im_N_USE');
%                tstat = Im_mat_USE./(Im_mat_std_USE./sqrt(Im_N_USE2));
%                p_ttest = 1 - tcdf(tstat,Im_N_USE2-1);
%                H_ttest = min(p_ttest,[],2) < alpha2;

                Mask_cells_USE = H_ttest .* Mask_nosmallbits .* Mask_crossection; % 
                %force last value = 1 (neuropil)
                Mask_cells_USE(end) = 1;

                Mask_cells_USE = H_ttest;

                if ~isempty(find(isnan(Im_mat_USE(:,1:NstimtypesMA))))
                    Mask_cells_USE(:) = 0;
                end
                
                %                eval(['RegInfo.Stats.tstat_',str_run,' = tstat;'])
                eval(['RegInfo.Stats.p_ttest_',str_run,' = p_ttest;'])
                eval(['RegInfo.Stats.H_ttest_',str_run,' = H_ttest;'])
                eval(['RegInfo.Stats.Mask_cells_USE_',str_run,' = Mask_cells_USE;'])

                Mask_cells_USE_ALL2(:,PLOTRUN_NORUN_RUN_ALLRUN) = Mask_cells_USE;
            else
                Mask_cells_USE = Mask_cells_USE_ALL2(:,PLOTRUN_NORUN_RUN_ALLRUN);
            end
            
            for count_cell = 1:Ncells;
            
                    if Mask_cells_USE(count_cell) == 1
                        a = Im_mat_USE(count_cell,:);
                        b = reshape(a(1:NSF*NTF),NSF,NTF);
                        b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
                        data = b2';


                    %for now, clean up analysis by shifting curves up so
                    %all is positive.
                    %all <0 set to 0
                    ind0 = find(data<0);
                    if ~isempty(ind0)
                        %data(ind0) = 0;
                        data = data - min(min(data));
                    end

                        if count_shuf == 0
                            PLOTIT_FIT = 0;
                            SAVEALLDATA = 1;
                            Fit_2Dellipse_MA2
                            eval(['Fit_struct(count_cell).True.s_',str_run,' = s;']);


                            %                        imagesc(data'); colormap('gray')
                            %                        hi = 1
                            %                        pause(100)
                        else
                            SAVEALLDATA = 0;
                            PLOTIT_FIT = 0;
                            Fit_2Dellipse_MA2
                            %                        hi = 1
                            eval(['Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,' = s;']);
%                            if count_shuf == 1
%                                hi = 1
%                                pause(100)
%                            end
                        end
                    title([str_run,': cell #',num2str(count_cell)]);
                end
            end
        end

    end
end

RegInfo.Stats.Fit_struct = Fit_struct;
RegInfo.Stats.ind_SFuse = ind_SFuse;
RegInfo.Stats.ind_TFuse = ind_TFuse;
%%
if isfield(RegInfo,'Stats') == 0
    RegInfo.Stats = [];
end
if isfield(RegInfo.Stats,'Fit_struct') == 0
    RegInfo.Stats.Fit_struct = [];
end
if isfield(RegInfo.Stats.Fit_struct,'Shuf') == 0
    RegInfo.Stats.Fit_struct.Shuf = 1;
    RegInfo.Stats.Fit_struct.True = 1;
end

%now collect all the fits into on vector: 
for count_run = 1:length(PLOTRUN_NORUN_RUN_ALLRUN_vec)
    PLOTRUN_NORUN_RUN_ALLRUN = PLOTRUN_NORUN_RUN_ALLRUN_vec(count_run);
    if PLOTRUN_NORUN_RUN_ALLRUN ==1
        str_run = 'run';
    elseif PLOTRUN_NORUN_RUN_ALLRUN ==2
        str_run = 'norun';
    else
        str_run = 'allrun';
    end
    
    eval(['Im_N_',str_run,' = zeros(NstimtypesMA,1);']);
    eval(['Im_mat_',str_run,' = zeros(Ncells,NstimtypesMA);']);
    eval(['Im_mat_std_',str_run,' = zeros(Ncells,NstimtypesMA);']);


    
    Mask_cells_USE = Mask_cells_USE_ALL2(:,PLOTRUN_NORUN_RUN_ALLRUN);
    ind = find(Mask_cells_USE);
    %    ind = find(Mask_cells_USE_ALL2(:,PLOTRUN_NORUN_RUN_ALLRUN));
    if length(ind) > 0
        %add a new one for avg timecourse and std of avg timecourse:
        %NEWADD MA 110531
%        Ntc = Ntot1_USE; %size(tc_resp2,2);
        Ntc =length(ind_TC); %size(tc_resp2,2);
        
        eval(['Im_mat_tc_',str_run,' = zeros(Ncells,Ntc,NstimtypesMA);']);
        eval(['Im_mat_tc_std_',str_run,' = zeros(Ncells,Ntc,NstimtypesMA);']);

        eval(['tmp = Fit_struct(ind(1)).Shuf(1).s_',str_run,'.x;']);
%    eval(['tmp = RegInfo.Stats.Fit_struct(1).Shuf(1).s_',str_run,'.x;']);
    sz_fit = size(tmp);
    fit_shuf_vec = zeros(Ncells,sz_fit(2)+4,Nshuf);
    fit_true_vec = zeros(Ncells,sz_fit(2)+4);
    for count = 1:NstimtypesMA
        ind_run_1 = eval(['Ind_struct(count).',str_run,'_trials;']);
        %if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
        %    ind_run_1 = ind_run_1(randsample(length(ind_run_1),length(ind_run_1),1));
        %end

        eval(['Im_N_',str_run,'(count) = length(ind_run_1);']);
        if length(ind_run_1)>0
            eval(['Im_mat_',str_run,'(:,count) = mean(tc_resp(ind_run_1,:),1);']);
            eval(['Im_mat_std_',str_run,'(:,count) = std(tc_resp(ind_run_1,:),[],1);']);
  
           %add a new one for avg timecourse and std of avg timecourse: 
           %NEWADD MA 110531
           eval(['Im_mat_tc_',str_run,'(:,:,count) = squeeze(mean(tc_resp2(:,ind_run_1,:),2))'';']);
           eval(['Im_mat_tc_std_',str_run,'(:,:,count) = squeeze(std(tc_resp2(:,ind_run_1,:),[],2))'';']);

        else
            %make all stimtypes where there is no data equal to NaN
            eval(['Im_mat_',str_run,'(:,count) = NaN;']);
            eval(['Im_mat_std_',str_run,'(:,count) = NaN;']);

            %add a new one for avg timecourse and std of avg timecourse: 
           %NEWADD MA 110531
           eval(['Im_mat_tc_',str_run,'(:,:,count) = NaN;']);
           eval(['Im_mat_tc_std_',str_run,'(:,:,count) = NaN;']);
           
        end
    end

    
    eval(['RegInfo.Stats.Im_mat_',str_run,' = Im_mat_',str_run,';'])
    eval(['RegInfo.Stats.Im_mat_std_',str_run,' = Im_mat_std_',str_run,';'])
  
            %add a new one for avg timecourse and std of avg timecourse: 
           %NEWADD MA 110531
    eval(['RegInfo.Stats.Im_mat_tc_',str_run,' = Im_mat_tc_',str_run,';'])
    eval(['RegInfo.Stats.Im_mat_tc_std_',str_run,' = Im_mat_tc_std_',str_run,';'])  
    
    eval(['RegInfo.Stats.Im_N_std_',str_run,' = Im_N_',str_run,';'])
    RegInfo.Stats.Ind_struct = Ind_struct;
    %RegInfo.Stats.Im_mat_std_USE_true = Im_mat_std_USE_true;
    %RegInfo.Stats.Im_N_USE_true = Im_N_USE_true;

    %now extract fits
    for count_cell = 1:Ncells
        if Mask_cells_USE(count_cell) == 1
            if ~isempty(RegInfo.Stats.Fit_struct(count_cell).True)                
                eval(['tmp = RegInfo.Stats.Fit_struct(count_cell).True.s_',str_run,'.x;']);
                eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).True.s_',str_run,'.SFhicut_50];']);
                eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).True.s_',str_run,'.TFhicut_50];']);
                eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).True.s_',str_run,'.SFhicut_10];']);
                eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).True.s_',str_run,'.TFhicut_10];']);
                fit_true_vec(count_cell,:) = tmp;
            end
        end
    end

    for count_shuf = 1:Nshuf
        for count_cell = 1:Ncells
            if Mask_cells_USE(count_cell) == 1
                if ~isempty(RegInfo.Stats.Fit_struct(count_cell).Shuf)
                    eval(['tmp = RegInfo.Stats.Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,'.x;']);
                    eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,'.SFhicut_50];']);
                    eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,'.TFhicut_50];']);
                    eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,'.SFhicut_10];']);
                    eval(['tmp = [tmp RegInfo.Stats.Fit_struct(count_cell).Shuf(count_shuf).s_',str_run,'.TFhicut_10];']);
                    %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                    fit_shuf_vec(count_cell,:,count_shuf) = tmp;
                end
            end
        end
    end
    
    Npars = size(fit_shuf_vec,2);
    lbub_fits = zeros(Ncells,Npars,5);
    alpha_bound = .025;
    for count_cell = 1:Ncells
        if Mask_cells_USE(count_cell) == 1
            for count2 = 1:Npars
                tmp = squeeze(fit_shuf_vec(count_cell,count2,:));
                [i,j] = sort(tmp);
                ind_shuf_lb = ceil(Nshuf*alpha_bound);
                ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
                lbub_fits(count_cell,count2,1) = i(ind_shuf_lb);
                lbub_fits(count_cell,count2,2) = i(ind_shuf_ub);
                lbub_fits(count_cell,count2,3) = mean(i); 
                lbub_fits(count_cell,count2,5) = std(i);

            end
            %now take means from truedata fit:
            lbub_fits(count_cell,:,4) = fit_true_vec(count_cell,:);
        end 
    end
    
    eval(['RegInfo.Stats.lbub_fits_',str_run,' = lbub_fits;'])
    eval(['RegInfo.Stats.alpha_bound = alpha_bound;']);

    
  
    
    
    %now
    
    
    eval(['RegInfo.Stats.fit_shuf_vec_',str_run,' = fit_shuf_vec;'])
    eval(['RegInfo.Stats.fit_true_vec_',str_run,' = fit_true_vec;'])

    end

%    eval(['RegInfo.Stats.Fit_struct_',str_run,' = Fit_struct;'])
%    eval(['RegInfo.Stats.ORI_vec0_',str_run,' = ORI_vec0; '])
%    eval(['RegInfo.Stats.SF_vec0_',str_run,' = SF_vec0;'])
%    eval(['RegInfo.Stats.alpha2_',str_run,' = alpha2;'])
end


RegInfo.Stats.Fit_struct = Fit_struct;
RegInfo.Stats.TF_vec0 = TF_vec0; 
RegInfo.Stats.SF_vec0 = SF_vec0; 
RegInfo.Stats.alpha2 = alpha2;
%RegInfo.Stats.tstat = tstat; 
RegInfo.Stats.p_ttest = p_ttest;
RegInfo.Stats.H_ttest = H_ttest;
RegInfo.Stats.Mask_cells_USE = Mask_cells_USE;


%RegInfo.Stats.Im_mat_USE_true = Im_mat_USE_true;
%RegInfo.Stats.Im_mat_std_USE_true = Im_mat_std_USE_true;
%RegInfo.Stats.Im_N_USE_true = Im_N_USE_true;
    

alloutfname=RegInfo.All.alloutfname;
save([alloutfname,'RegInfo.mat'],'RegInfo');

%RegInfo.Stats.Fit_struct(1).Shuf(100)


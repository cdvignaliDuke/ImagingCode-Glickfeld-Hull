
'0209_CFJ_F'; %name of Mfile with all parameters..
%eval([Mfile]);
%assemble a directory of mat files in order, the write to a tif file (or 2
%tiff files)

%Mfile = 'PARAMS_110205_DR1_F'; %name of Mfile with all parameters..
%Mfile = 'PARAMS_110205_DR3_F'; %name of Mfile with all parameters..
%Mfile = 'PARAMS_11

%PWD = 'F:\users\demetris\2pdata\110209_CFJ\';
%can specify PWD or get PWD from running the PARAMS file first

plane = 1; %plane number
count_protocol =1; %=1 => 1st protocol ran that day on that plane

RUNREG = 0; %run analysis on registered version the files
%WIN_PRE = [8 10];
%WIN_POST = [1 5];


% 5 sec on window
WIN_POST = [1 5];
% 1 sec on window
%WIN_POST = [0 2];

WIN_PRE = [-2 0];


%window in seconds

PWD_analysis = [PWD,'analysis\'];
if isempty(dir(PWD_analysis))
    eval(['!mkdir ',PWD_analysis]);
end


filename_RandStim = [PWD,'RandStim_info_plane',num2str(plane),'.mat'];
load(filename_RandStim);
str_prot = ['prot',num2str(count_protocol),'_plane',num2str(plane)];

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
Toff1 = RandStim_info.prot(count_protocol).Toff1;
Ton1 = RandStim_info.prot(count_protocol).Ton1;
Noff1 = RandStim_info.prot(count_protocol).Noff1;
Non1 = RandStim_info.prot(count_protocol).Non1;
Ntot1 = RandStim_info.prot(count_protocol).Ntot1;
DIR_rand0 = RandStim_info.general.DIR_rand0;
file_mat = RandStim_info.general.file_mat;
randvec = RandStim_info.general.randvec;

t = [1/Fs1:1/Fs1:2*Ntot1/Fs1] - Noff1/Fs1;
Nstimtypes = length(TFvec);
Nruns = length(userun);

for count_run = 1:length(userun)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% load imaging data
    
    file_USE = deblank(file_mat(userun(count_run),:));
    if RUNREG == 1
        tmp = findstr(file_USE,'.tif');
        tmp2 = file_USE(tmp(end)-3:tmp(end)-1);
        if strcmp(tmp2,'reg') == 0
            file_USE = [file_USE(1:(tmp(end)-1)),'_reg.tif'];
        end
    end
    stack = readtiff([PWD,file_USE]);
    sz1 = size(stack);
    Nframes1 = size(stack,3);
    %add an extra stimulus presentation to the end, for analysis purposes..
    stack2 = zeros(sz1(1),sz1(2),sz1(3) + Ntot1);
    stack2(:,:,1:Nframes1) = stack;
    stack2(:,:,Nframes1 + [1:Ntot1]) = stack(:,:,1:Ntot1);
    stack = stack2;
    clear('stack2');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load stimulus randomization:
    tmp = randvec(userun(count_run),2);
    tmp = rem(tmp,8);
    tmp(find(tmp==0)) = 8;
    tmp_DIR = deblank(DIR_rand0(count_protocol,:));
    DIR_rand = dir([tmp_DIR,['*_',num2str(tmp),'__','*.mat']]);
    
    %DIR_rand = dir([tmp_DIR]);
    if length(DIR_rand) >= 1
        file_rand = [tmp_DIR,DIR_rand(1).name];
    else
        error('error: wrong # of matches');
        return
    end
    load([file_rand]);
    if size(Seq,2) == 7 %a quick fix for 110205 ONLY
        Seq = [Seq, ones(size(Seq,1),1)*36 ones(size(Seq,1),1)*-18];
    end

    %        Nstim1 = Nframes1./(Ntot1*Fs1);
   % Nstim1 = min(size(Seq,1),Nframes1./(Ntot1)); %%NOTE : CHECK THIS
   
   % Nstim1 = min(size(Seq,1),Nframes1./(Ntot1)); %%NOTE : CHECK THIS
 
       Nstim1 = floor(Nframes1./(Ntot1)); %%NOTE : CHECK THIS
 
    Seqposition = [];
    Seqposition.ind = [];
    Seqposition2 = [];
    New_Seqposition= [];
    New_Seqposition.ind= [];
    Nstimtypes1 = length(TFvec);
    
    for count_seq = 1:Nstimtypes1
        ind = find( ...
            Seq(1:Nstim1,4) == TFvec(count_seq) & ...
            Seq(1:Nstim1,5) == SFvec(count_seq) & ...
            Seq(1:Nstim1,6) == oris(count_seq) & ...
            Seq(1:Nstim1,7) == contrastvec(count_seq) & ...
            Seq(1:Nstim1,8) == StimXpos(count_seq) & ...
            Seq(1:Nstim1,9) == StimYpos(count_seq));
        %ind = reshape(ind,length(ind),1));
        
        Seqposition(count_seq).ind = ind;
        Stimvec =  [TFvec(count_seq) SFvec(count_seq) oris(count_seq) contrastvec(count_seq) ...
            StimXpos(pos(count_seq)) StimYpos(pos(count_seq))];
        Seqposition(count_seq).TFSFetc = Stimvec;
        Seqposition2 = [Seqposition2; [ones(length(ind),1)*[count_seq Stimvec] ind]];
    end
    if count_run == 1
        Big_Seqposition = [];
        Big_Seqposition.ind = [];
        Nstim_runningtotal = 0;
    end

    for count_seq = 1:Nstimtypes1
        if count_run == 1
            Big_Seqposition(count_seq).TFSFetc =  Seqposition(count_seq).TFSFetc;
            Big_Seqposition(count_seq).ind = Seqposition(count_seq).ind;
        end
        if count_run > 1
            New_Seqposition(count_seq).ind = Seqposition(count_seq).ind + Nstim_runningtotal;
            Big_Seqposition(count_seq).ind = [Big_Seqposition(count_seq).ind; New_Seqposition(count_seq).ind];
        end    
    end
    
    Nstim_runningtotal = Nstim_runningtotal + Nstim1;
    
    %werite a new, truncated
     
    %stack2 = stack(:,:,1:(Nstim1*Ntot1));
    %write this to disk, adding something to end of filename
    %ind2 = findstr(file_USE,'.tif');
    
    
    %file_USE2 = [file_USE(1:ind2(1)-1),'_trunc.tif']; 
    %writetif(stack2,[PWD,file_USE2]);
    
    
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %reorder the stimuli in the order in
    %stack2 = zeros(Nstim1, sz1(1), sz1(2),Ntot1);
    %for count = 1:Nstim1
    %    ind2 = Seqposition2(count,8);
    %    ind_use = (ind2-1)*Ntot1 + [1:Ntot1];
    %    stack2(count,:,:,:) = stack(:,:,ind_use);
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%now average all the stimuli of a certain type together (dF/Fmaps
    %%%%and raw TCs
    stack3 = zeros(sz1(1), sz1(2),Nstimtypes1*Ntot1);
    stack4 = zeros(sz1(1), sz1(2),Nstimtypes1);

    ind_pre = find(t>WIN_PRE(1) & t<=WIN_PRE(2));
    ind_post = find(t>WIN_POST(1) & t<=WIN_POST(2));

    Npertype = zeros(Nstimtypes1,1);
    for count_seq = 1:Nstimtypes1
        tmp = zeros(sz1(1),sz1(2),Ntot1);
        tmpB = zeros(sz1(1),sz1(2));
        ind = Seqposition(count_seq).ind;
        Npertype(count_seq) = length(ind);
        for count2 = 1:length(ind)
            ind2 = ind(count2);
            ind_use = (ind2-1)*Ntot1 + [1:Ntot1];
            tmp0 =  double(squeeze(stack(:,:,ind_use)));
            tmp = tmp + tmp0;

            ind_use_pre = (ind2-1)*Ntot1 + [ind_pre];
            ind_use_post = (ind2-1)*Ntot1 + [ind_post];
            tmpB_pre = squeeze(mean(stack(:,:,ind_use_pre),3));
            tmpB_post = squeeze(mean(stack(:,:,ind_use_post),3));
            tmpB_dF = (tmpB_post - tmpB_pre)./tmpB_pre;
            tmpB = tmpB + tmpB_dF;            
        end
        tmp = tmp./length(ind);
        tmpB = tmpB./length(ind);
        stack3(:,:,(count_seq-1)*Ntot1 + [1:Ntot1]) = tmp;
        stack4(:,:,count_seq) = tmpB;

        if count_run == 1
            Big_stack4 = zeros(length(userun),sz1(1),sz1(2),Nstimtypes);
            Big_stack3 = zeros(length(userun),sz1(1),sz1(2),Nstimtypes*Ntot1);
        end
        Big_stack4(count_run,:,:,:) = stack4;
        Big_stack3(count_run,:,:,:) = stack3;
    end
    
    Npertype(count_seq) = length(ind);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now output tif files fr looking at average dF/F
    ind = findstr(file_USE,'.tif');
    file_USE2 = [PWD_analysis,file_USE(1:(ind(1)-1)),'_tc_Nstim',num2str(Nstimtypes1),'.tif'];
    writetiff(uint16(stack3),file_USE2,'uint16');
    
    ind = findstr(file_USE,'.tif');
    if strfind(file_USE, '_plane') > 1
        amin = 12;
        bmin = 13;
    elseif strfind(file_USE, '_reg') > 1
        amin = 9;
        bmin = 10; 
    else
        amin = 1;
        bmin= 2;
    end
    
    fn_out = [PWD_analysis, file_USE(1:(ind(1)-amin)), '_Seqposition.mat'];
    save(fn_out, 'Seqposition');

    if isnumber(file_USE((ind(1)-bmin))) %MA ADD
       bmin = bmin + 1;
    end
    
    fn_out = [PWD_analysis, file_USE(1:(ind(1)-bmin)), num2str(userun) '_Big_Seqposition.mat'];
    save(fn_out, 'Big_Seqposition');
    
    ind = findstr(file_USE,'.tif');
    file_USE2 = [PWD_analysis,file_USE(1:(ind(1)-amin)),'_dF_F_Nstim',num2str(Nstimtypes1), ...
        '_POST_',num2str(WIN_POST(1)),'_',num2str(WIN_POST(2)),'.tif'];
    writetiff((stack4),file_USE2,'double');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now output tif files fr looking at average dF/F, averaged across
%%% all runs in a protocol
stack4_mean = squeeze(mean(Big_stack4,1));
stack3_mean = squeeze(mean(Big_stack3,1));
if length(userun)>1
ind = findstr(file_USE,'.tif');
file_USE2 = [PWD_analysis,file_USE(1:(ind(1)-1)),'_dF_F_Nstim',num2str(Nstimtypes1), ...
    '_POST_',num2str(WIN_POST(1)),'_',num2str(WIN_POST(2)),'_MEAN.tif'];
%   file_USE2 = [PWD,file_USE(1:(ind(1)-1)),'_dF_F_Nstim',num2str(Nstimtypes1),'_MEAN.tif'];
writetiff((stack4_mean),file_USE2,'double');

ind = findstr(file_USE,'.tif');
file_USE2 = [PWD_analysis,file_USE(1:(ind(1)-1)),'_tc_Nstim',num2str(Nstimtypes1), ...
    '_POST_',num2str(WIN_POST(1)),'_',num2str(WIN_POST(2)),'_TCMEAN.tif'];
%file_USE2 = [PWD,file_USE(1:(ind(1)-1)),'_dF_F_Nstim',num2str(Nstimtypes1),'_MEAN.tif'];
%writetiff((stack3_mean),file_USE2,'double');
end
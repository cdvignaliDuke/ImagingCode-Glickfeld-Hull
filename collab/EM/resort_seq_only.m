plane = 1; %plane number
count_protocol = 1; %=1 => 1st protocol ran that day on that plane

RUNREG = 0; %run analysis on registered version the files
WIN_PRE = [8 10];
%WIN_POST = [1 4]
WIN_POST = [1 5]; %window in seconds
%WIN_POST = [5 8]; %window in seconds
%WIN_POST = [2 5]; %window in seconds
%WIN_POST = [10 18]; %window in seconds

PWD_analysis = fullfile(PWD,'analysis\');
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

    Nframes1 = sizetiff([PWD,file_USE])

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
            Big_Seqposition(count_seq).ind = Seqposition(count_seq).ind;
            Big_Seqposition(count_seq).TFSFetc =  Seqposition(count_seq).TFSFetc;
        end
        if count_run > 1
            New_Seqposition(count_seq).ind = Seqposition(count_seq).ind + Nstim_runningtotal;
            Big_Seqposition(count_seq).ind = [Big_Seqposition(count_seq).ind; New_Seqposition(count_seq).ind];
            Big_Seqposition(count_seq).TFSFetc =  Seqposition(count_seq).TFSFetc;
        end    
    end
    
    Nstim_runningtotal = Nstim_runningtotal + Nstim1;
  
    ind = findstr(file_USE,'.tif');
    if ind>20
        min = 9;
        bmin = 10;
    else
        min = 1;
        bmin= 2;
    end
    
    fn_out = [PWD_analysis, file_USE(1:(ind(1)-min)), '_Seqposition.mat'];
    save(fn_out, 'Seqposition');

    fn_out = [PWD_analysis, file_USE(1:(ind(1)-bmin)), num2str(userun) '_Big_Seqposition.mat'];
    save(fn_out, 'Big_Seqposition');
end

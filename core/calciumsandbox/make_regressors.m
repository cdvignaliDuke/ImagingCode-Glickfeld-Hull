function regressors = make_regressors (s, lowcut_off_frame, Ca_ker)
%  s.Non, s.Noff, s.Nstim_per_run, s.Nrep

nframes_per_stim =(s.Non+s.Noff);
nframes_per_run = nframes_per_stim * s.Nstim_per_run;
Nframes = nframes_per_run * s.Nrep;

Nlowcut=floor(log2(Nframes/lowcut_off_frame))+2;
Nregressors=s.Nstim_per_run + Nlowcut*2 + 1;

regressors=zeros(Nframes,Nregressors);

stim_frames = getepochs (s.Noff, s.Non, s.Nstim_per_run, s.Nrep);

for st=1:s.Nstim_per_run
    for tr=1:s.Nrep
        regressors(stim_frames{tr,st},st)=1;
    end
end


for st=1:s.Nstim_per_run
    temp=conv(Ca_ker,regressors(:,st));
    regressors(:,st)=temp(1:Nframes);
end


for i=1:Nlowcut
    period=Nframes/(2^(i-2));
    regressors(:,s.Nstim_per_run+i*2-1)=cos([1:Nframes]*2*pi/period);
    regressors(:,s.Nstim_per_run+i*2)=sin([1:Nframes]*2*pi/period);
end

regressors(:,Nregressors)=1;
end




function [ProcTimeCourses,data_tables]=tcProcess (timeCourses, s, baseline_percentile, leica)
%TCPROCESS Kenichi's standard processing stream
% [PROCTIMECOURSES,DATA_TABLES]=TCPROCESS (TIMECOURSES, S, BASELINE_PERCENTILE, LEICA)
% timeCourses Nt * Ncells
% structure of exp parameters: s.Non, s.Noff, s.Nstim_per_run, s.Nrep,
% lowcut_sec: low cut filter (sec) default: trial (run) duration *2
% baseline_percentile: set the baseline at this value of percentile from the lowest. default: 30 (%)
%
%   10.31.2008 Kenichi Ohki


if nargin < 3 
    baseline_percentile = 30 ;
end
if nargin < 4
    leica = 0 ;
end

pre_stim=1;

%lowcut_frame=lowcut_sec*s.framerate;
    
% setup parameters
[nframes, ncells]=size(timeCourses);
expt=[];
expt.dummy=[];
expt = setupstim(expt,s.Non, s.Noff, s.Nstim_per_run, nframes);
expt.stims = getepochs (s.Noff, s.Non, s.Nstim_per_run, s.Nrep,pre_stim,leica);

% preprocessing of time courses

ProcTimeCourses.lowcut = tcLowCut (timeCourses, expt.trialdur*2, 'gaussian', 1);
ProcTimeCourses.lowcut(expt.dur+1:end,:)=[];
ProcTimeCourses.lowcut_avg=squeeze(sum(reshape(ProcTimeCourses.lowcut,expt.trialdur,expt.ntrials,ncells),2))./expt.ntrials;


temp=sort(ProcTimeCourses.lowcut);
base = temp(round(nframes*baseline_percentile/100),:);
ProcTimeCourses.norm = (ProcTimeCourses.lowcut-repmat(base,expt.dur,1))./repmat(base,expt.dur,1)*100;   % percent signal change
ProcTimeCourses.norm_avg=squeeze(sum(reshape(ProcTimeCourses.norm,expt.trialdur,expt.ntrials,ncells),2))./expt.ntrials;
data_tables = tcEpochAverage(ProcTimeCourses.lowcut,expt.stims);


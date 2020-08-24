% 
% spikes.spiketime_S = [spikes_9.spiketime_S;spikes_11.spiketime_S;spikes_13.spiketime_S];
% spikes.isi_ms = [spikes_9.isi_ms;spikes_11.isi_ms;spikes_13.isi_ms];
% spikes.isi_vilation = [spikes_9.isi_vilation, spikes_11.isi_vilation,spikes_13.isi_vilation];
% spikes.cluster = [spikes_9.cluster;spikes_11.cluster;spikes_13.cluster];
% spikes.channelposlist = spikes_9.channelposlist;
% spikes.template = cat(3,cat(3,spikes_9.template, spikes_11.template),spikes_13.template);
% spikes.template_pos = cat(1,cat(1,spikes_9.template_pos, spikes_11.template_pos),spikes_13.template_pos);
% spikes.template_big = cat(2,cat(2,spikes_9.template_big, spikes_11.template_big),spikes_13.template_big);
% spikes.chID = cat(1,cat(1,spikes_9.chID, spikes_11.chID),spikes_13.chID);
% 
% spikes.waveform.mean = [spikes_9.waveform.mean, spikes_11.waveform.mean, spikes_13.waveform.mean];
% spikes.waveform.all = [spikes_9.waveform.all, spikes_11.waveform.all, spikes_13.waveform.all];
% spikes.waveform.pos = cat(1,cat(1,spikes_9.waveform.pos, spikes_11.waveform.pos),spikes_13.waveform.pos);
% spikes.waveform.big = cat(2,cat(2,spikes_9.waveform.big, spikes_11.waveform.big),spikes_13.waveform.big);
% spikes.waveform.chID = [spikes_9.waveform.chID; spikes_11.waveform.chID; spikes_13.waveform.chID];
% spikes.waveform.bigall = [spikes_9.waveform.bigall, spikes_11.waveform.bigall,spikes_13.waveform.bigall]; 
% 
% spikes.waveform.stats.pt_ms = [spikes_9.waveform.stats.pt_ms;spikes_11.waveform.stats.pt_ms;spikes_13.waveform.stats.pt_ms];
% spikes.waveform.stats.pt_ratio = [spikes_9.waveform.stats.pt_ratio;spikes_11.waveform.stats.pt_ratio;spikes_13.waveform.stats.pt_ratio];
% spikes.waveform.stats.pt_distance = [spikes_9.waveform.stats.pt_distance;spikes_11.waveform.stats.pt_distance;spikes_13.waveform.stats.pt_distance];
% spikes.waveform.stats.asym = [spikes_9.waveform.stats.asym;spikes_11.waveform.stats.asym;spikes_13.waveform.stats.asym];
% spikes.waveform.stats.hfw = [spikes_9.waveform.stats.hfw;spikes_11.waveform.stats.hfw;spikes_13.waveform.stats.hfw];
% spikes.waveform.unit = spikes_9.waveform.unit; 
% 
% 
% spikes.wave_stats.pt_ms = [spikes_9.wave_stats.pt_ms;spikes_11.wave_stats.pt_ms;spikes_13.wave_stats.pt_ms];
% spikes.wave_stats.pt_ratio = [spikes_9.wave_stats.pt_ratio;spikes_11.wave_stats.pt_ratio;spikes_13.wave_stats.pt_ratio];
% spikes.wave_stats.pt_distance = [spikes_9.wave_stats.pt_distance;spikes_11.wave_stats.pt_distance;spikes_13.wave_stats.pt_distance];
% spikes.wave_stats.asym = [spikes_9.wave_stats.asym;spikes_11.wave_stats.asym;spikes_13.wave_stats.asym];
% spikes.wave_stats.hfw = [spikes_9.wave_stats.hfw;spikes_11.wave_stats.hfw;spikes_13.wave_stats.hfw];
% 
% clearvars -except spikes

spikes.spiketime_S = [spikes_9.spiketime_S;spikes_11.spiketime_S];
spikes.isi_ms = [spikes_9.isi_ms;spikes_11.isi_ms];
spikes.isi_vilation = [spikes_9.isi_vilation, spikes_11.isi_vilation];
spikes.cluster = [spikes_9.cluster;spikes_11.cluster];
spikes.channelposlist = spikes_9.channelposlist;
spikes.template = cat(3,spikes_9.template, spikes_11.template);
spikes.template_pos = cat(1,spikes_9.template_pos, spikes_11.template_pos);
spikes.template_big = cat(2,spikes_9.template_big, spikes_11.template_big);
spikes.chID = cat(1,spikes_9.chID, spikes_11.chID);

spikes.waveform.mean = [spikes_9.waveform.mean spikes_11.waveform.mean];
spikes.waveform.all = [spikes_9.waveform.all spikes_11.waveform.all];
spikes.waveform.pos = cat(1,spikes_9.waveform.pos, spikes_11.waveform.pos);
spikes.waveform.big = cat(2,spikes_9.waveform.big, spikes_11.waveform.big);
spikes.waveform.chID = [spikes_9.waveform.chID; spikes_11.waveform.chID];
spikes.waveform.bigall = [spikes_9.waveform.bigall spikes_11.waveform.bigall]; 

spikes.waveform.stats.pt_ms = [spikes_9.waveform.stats.pt_ms;spikes_11.waveform.stats.pt_ms];
spikes.waveform.stats.pt_ratio = [spikes_9.waveform.stats.pt_ratio;spikes_11.waveform.stats.pt_ratio];
spikes.waveform.stats.pt_distance = [spikes_9.waveform.stats.pt_distance;spikes_11.waveform.stats.pt_distance];
spikes.waveform.stats.asym = [spikes_9.waveform.stats.asym;spikes_11.waveform.stats.asym];
spikes.waveform.stats.hfw = [spikes_9.waveform.stats.hfw;spikes_11.waveform.stats.hfw];
spikes.waveform.unit = spikes_9.waveform.unit; 



spikes.wave_stats.pt_ms = [spikes_9.wave_stats.pt_ms;spikes_11.wave_stats.pt_ms];
spikes.wave_stats.pt_ratio = [spikes_9.wave_stats.pt_ratio;spikes_11.wave_stats.pt_ratio];
spikes.wave_stats.pt_distance = [spikes_9.wave_stats.pt_distance;spikes_11.wave_stats.pt_distance];
spikes.wave_stats.asym = [spikes_9.wave_stats.asym;spikes_11.wave_stats.asym];
spikes.wave_stats.hfw = [spikes_9.wave_stats.hfw;spikes_11.wave_stats.hfw];
% 
clearvars -except spikes


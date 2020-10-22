function [mu,sigma,CV,rate] = calc_spike_var(t,spikes,type,N_min)
%Comments added by JH 191106
%adapted from a function found on stack Exchange
%t = vector of timepoints within the sweeps you wish to calculate variables
%spikes = array of spikes times for each sweep. dim1=time dim2=sweep#
%N_min   # of spikes must be > N_min for a given sweep to count that sweep

%   MU = VCN_REGULARITY(T,SPIKES) returns the mean Inter-Spike Interval
%   (ISI) at the times contained in the vector T for the spike trains 
%   contained in the columns of the array SPIKES. Each column of SPIKES 
%   should consist of an ordered set of spike times (ordering is not 
%   checked) padded with NaNs. MU will contain NaNs for times at which it 
%   is not defined.
%
%   [MU,SIGMA,CV,RATE] = VCN_REGULARITY(...) also returns the standard 
%   deviation SIGMA, the coefficient of variation CV and the mean
%   inter-spike rate RATE.
%
%   VCN_REGULARITY(T,SPIKES,'hybrid') uses a hybrid algorithm to mimic the
%   results of the regularity analysis used by Young et al (1988). ISIs
%   only contribute to their statistics if they are within a certain
%   distance of their commencement, and that distance is controlled by the
%   internal parameter pc_hybrid. VCN_REGULARITY(T,SPIKES,'exact') gives the 
%   default exact calculation.
%
%   VCN_REGULARITY(T,SPIKES,TYPE,N) replaces the outputs with NaNs at times
%   for which fewer than N intervals are in progress.  If N is not
%   specified it is set to one fifth of the number of columns in SPIKES.


N_sweeps = size(spikes,2);
R = nan(size(spikes,1),N_sweeps);

%for each spike within each sweep, find the interval between the current
%spike and the next spike
for j = 1:N_sweeps
    for i = 1:length(find(~isnan(spikes(:,j)))) - 1 %for a given spike time...
        if spikes(i+1,j)<t(end) & spikes(i,j) > t(1) %first check if both this spike and the following spike are within analysis window
            R(i,j) = spikes(i+1,j) - spikes(i,j); %then calculate interval between this spike and the next spike. 
        end
    end
end

R = reshape(R,1,size(R,1)*size(R,2));
R = R(find(~isnan(R)));

if length(R) > 5
    mu = mean(R);
    sigma = std(R);
    rate = NaN;
    CV = sigma./mu;
else
    mu = [];
    sigma = [];
    rate = [];
    CV = [];
end
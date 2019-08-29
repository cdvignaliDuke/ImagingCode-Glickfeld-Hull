% Alex's code: identify spikes of raw fluorescence in 2P. 
function [spklogic] = detectspikes(data, threshold, refractoryperiod, minspkwidth)

% INPUTS
% data              = frames x cell
% threshold         = number of std dev above mean to set threshold
% refractoryperiod  = number of frames between acceptable spikes
% minspkwidth       = minimum number of frames for each spike

% OUTPUTS
% spkmat            = array containing spike indices (matches data size)

% calculate trial specific thresholds
data     = [diff(data,1); data(end,:)];
std_val  = threshold*std(data,0,1);
mean_val = mean(data,1);
thresh   = mean_val + std_val;

% threshold the data array based on trial specific thresholds
spklogic = bsxfun(@gt, data, thresh);

% loop through values to remove indices that are within same spiking event
for cells = 1:size(data,2)

   groups = bwconncomp(spklogic(:,cells));
   inds   = groups.PixelIdxList;

   for spk = 1:length(inds)

       if numel(inds{spk}) < minspkwidth
           spklogic(inds{spk}, cells) = false;

           continue
       end

       if numel(inds{spk}>1)
           [~,ind2keep] = max(data(inds{spk}, cells));

           if length(ind2keep) > 1
               keyboard
           end

           spklogic(inds{spk}(inds{spk} ~= inds{spk}(ind2keep)), cells) =  false;

       end

   end
end

% remove spiking events that are within refractory period of each other
for cells = 1:size(data,2)
   groups = bwconncomp(spklogic(:,cells));
   inds   = groups.PixelIdxList;

   if numel(inds)<=1 || isempty(inds)
       continue
   end

   ind2remove = false(1,length(inds));
   for i = 1:length(inds)-1
       inds_i  = inds{i};

       for j = (i + 1) :length(inds)
           inds_j  = inds{j};

           if (inds_j - inds_i) > refractoryperiod
               continue
           end

           spkrt_i = data(inds_i, cells);
           spkrt_j = data(inds_j, cells);

           if spkrt_i > spkrt_j
              ind2remove(j) = true;
           else
              ind2remove(i) = true; 
           end

       end
   end

   spklogic([inds{ind2remove}], cells) = false;

end 



end


%function to be used in WF_CRP_peak_mag_means to calculate peak magnitudes
function cell_out = min_trial_finder(cell_in);
%inputs cell arrays which contain the pre- or post- reward magnitude values
%for each trial for each animal
%outputs  another cell array but each animal has the same number of trials

%make a vector of all the non zero lengths of items in a cell 
cell_len = cellfun(@length, cell_in);
cell_len(find(cell_len<1)) = NaN;

%take the smallest value from that vector
min_cell_len = nanmin(cell_len);

%use linspace to select an evenly spaced sampl of that number of elements from each element of the cell array
cell_out = []; 
for session_num = 1:length(cell_in)
    %skip empty cells
    if isempty(cell_in{session_num})
        cell_out{session_num} = [];
    else
        lin_samp_ind = round(linspace(1, length(cell_in{session_num}), min_cell_len));
        lin_samp = cell_in{session_num}(lin_samp_ind);
        cell_out{session_num} = lin_samp;
    end
end



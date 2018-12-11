function output_cells = cell_cat_selector(cell_cats, all_cell_cats, ITI_lick_resp_cells, ITI_bout_resp_cells);
    if strcmp(cell_cats, 'allresp_cells')
        output_cells = all_cell_cats.allresp_cells;
    elseif strcmp(cell_cats, 'allresp_cells_pos')
        output_cells = all_cell_cats.allresp_cells_pos;
    elseif strcmp(cell_cats, 'cue_cells_pos')
        output_cells = all_cell_cats.cue_cells_pos;
    elseif strcmp(cell_cats, 'rew_cells_pos')
        output_cells = all_cell_cats.rew_cells_pos;
    elseif strcmp(cell_cats, 'allresp_cells_neg')
        output_cells = all_cell_cats.allresp_cells_neg;
    elseif strcmp(cell_cats, 'rew_cells_neg')
        output_cells = all_cell_cats.rew_cells_neg;
    elseif strcmp(cell_cats, 'cue_cells_neg')
        output_cells = all_cell_cats.cue_cells_neg;
    elseif strcmp(cell_cats, 'rew_cells')
        output_cells = all_cell_cats.rew_cells;
        
    elseif strcmp(cell_cats, 'OR_Cue_resp_cells_pos')
         output_cells = zeros(1,length(all_cell_cats.allresp_cells));
         output_cells([all_cell_cats.OR_Cue_resp_cells_pos]) = 1; 
    elseif strcmp(cell_cats, 'OR_Rew_resp_cells_pos')   
        output_cells = zeros(1,length(all_cell_cats.allresp_cells));
        output_cells([all_cell_cats.OR_Rew_resp_cells_pos]) = 1; 
    elseif strcmp(cell_cats, 'NR_Cue_resp_cells_pos')
        output_cells = zeros(1,length(all_cell_cats.allresp_cells));
        output_cells([all_cell_cats.NR_Cue_resp_cells_pos]) = 1;
    elseif strcmp(cell_cats, 'NR_Rew_resp_cells_pos')
        output_cells = zeros(1,length(all_cell_cats.allresp_cells));
        output_cells([all_cell_cats.NR_Rew_resp_cells_pos]) = 1;
    elseif strcmp(cell_cats, 'NR_Rew_resp_cells_neg')
        output_cells = zeros(1,length(all_cell_cats.allresp_cells));
        output_cells([all_cell_cats.NR_Rew_resp_cells_neg]) = 1;
        
    elseif strcmp(cell_cats, 'ITI_lick_resp_cells') & ~isnan(ITI_lick_resp_cells)
        output_cells = ITI_lick_resp_cells;
    elseif strcmp(cell_cats, 'ITI_lick_resp_cells') & isnan(ITI_lick_resp_cells)
        output_cells = all_cell_cats.allresp_cells;
    elseif strcmp(cell_cats, 'ITI_bout_resp_cells') & ~isnan(ITI_bout_resp_cells)
        output_cells = ITI_bout_resp_cells;
    elseif strcmp(cell_cats, 'ITI_bout_resp_cells') & isnan(ITI_bout_resp_cells)
        output_cells = all_cell_cats.allresp_cells;
    end
    
    return
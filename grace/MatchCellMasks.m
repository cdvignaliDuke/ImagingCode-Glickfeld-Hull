% variables
date = '191114';
ImgFolder = strvcat('003');
run = strvcat('000');
mouse = 'i1306';
doFromRef = 1;
ref_date = '191113';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);

% load reference and current cell masks
doFromRef;
    ref_str = ['runs-' ref_run];
if size(ref_run,1)>1
    ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
end

refMaskData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_mask_cell.mat']));
ref_mask_cell =  refMaskData.mask_cell;
ncells_ref = [length(unique(ref_mask_cell(:)))]-1;

currentMaskData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']));
current_mask_cell =  currentMaskData.mask_cell;
ncells_current = [length(unique(current_mask_cell(:)))]-1;

% create reference and current cell masks
isCellInRefMask = false(1,ncells_ref);
cellMatchInRefMask = nan(1,ncells_ref);
nCellMatchInRefMask = nan(1,ncells_ref);
fractCellMatchInRefMask = nan(1,ncells_ref);
for icell = 1:ncells_ref
    idx = current_mask_cell == icell;
    isCellInRefMask(icell) = any(ref_mask_cell(idx)>0);
    if isCellInRefMask(icell) == 1
        temp = unique(ref_mask_cell(idx));
        nCellMatchInRefMask(icell) = length(find(temp>0));
        temp_ref_pix = length(find(ref_mask_cell(idx)>0));
        fractCellMatchInRefMask(icell) = temp_ref_pix./length(idx);
        %cellMatchInRefMask(icell) = max(temp);
    end
end


isCellInCurrentMask = false(1,ncells_current);
for icell = 1:ncells_current
    idx = ref_mask_cell == icell;
    isCellInCurrentMask(icell) = any(current_mask_cell(idx)>0);
end

% load and create reference TCs
refTcData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_TCs.mat']));
ref_npSub_tc = refTcData.npSub_tc;
ref_tc_data = ref_npSub_tc(:,isCellInRefMask);
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_refTCs.mat']), 'ref_tc_data');

% load and create current TCs
currentTcData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']));
current_npSub_tc = currentTcData.npSub_tc;
current_tc_data = current_npSub_tc(:,isCellInCurrentMask);
save(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_currentTCs.mat']), 'current_tc_data');


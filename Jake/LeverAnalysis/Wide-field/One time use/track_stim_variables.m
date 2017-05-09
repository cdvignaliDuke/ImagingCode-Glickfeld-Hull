clear;
days = {'170215_img78', '170216_img78', '170219_img78', '170220_img78', '170221_img78', '170222_img78', '170223_img78', '170224_img78', '170301_img78', '170302_img78', '170303_img78', '170305_img78', '170306_img78', '170307_img78', '170308_img78'};
bx_source      = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
bx_outputs_dir = ['Z:\Analysis\Cue_reward_pairing_analysis\BxAndAnalysisOutputs\BxOutputs\'];
old_cd = cd; %save old cd so I can restore it later

rewardSoundCell = {};
targetSoundCell = {};
for ii = 1:length(days)
    bx_out_dir  = [bx_outputs_dir days{ii} '_bx_outputs'];
    b_data = get_bx_data(bx_source, days{ii});  %find the correct behavior file and loads it.
    if isfield(b_data, 'doRewardSound');
            rewardSoundCell = setfield(rewardSoundCell, ['i', days{ii}(1:6)], b_data.doRewardSound);
    else
         rewardSoundCell = setfield(rewardSoundCell, ['i', days{ii}(1:6)], 1);
    end
    targetSoundCell = setfield(targetSoundCell, ['i', days{ii}(1:6)], b_data.soundTargetAmplitude);
end
rewardSoundCell
targetSoundCell
%to be run within summary_PSTHs_CRP_JH and other scripts which need to load CRP variables

function [dest_sub, dest_sub_spikes, iti_lick_dir, data_dir, session, session_date, mouse_num, mouse_ID, rID] = get_sess_info(session, runID, crp_dir, ITI_dir);
    session_date = session(1:6);
    if session(end-2) == 'g'
        mouse_num = ['9', session(end-1:end)];
        mouse_ID = ['img', session(end-1:end)];
    elseif session(end-2) =='0'
        mouse_num = session(end-2:end);
        mouse_ID = ['img', session(end-2:end)];
    end
    
    %set pathnames
    for rID = 1:2
        if exist([crp_dir, session_date, '_', runID{rID}, '_', mouse_ID], 'file') == 7
            break
        end
    end
    if strcmp(session, '170426_img91') | strcmp(session, '170522_img89')
        rID=2;
    end
    dest_sub  = fullfile(crp_dir, [session_date, '_', runID{rID}, '_', mouse_ID],'\');
    dest_sub_spikes  = fullfile(dest_sub, 'spike_outputs', '\');
    if ~isempty(ITI_dir)
        iti_lick_dir = fullfile(ITI_dir, session, '\');
    else
         iti_lick_dir = [];
    end
    data_dir = fullfile('Z:\Data\2P_imaging', session, mouse_ID, '\');
end
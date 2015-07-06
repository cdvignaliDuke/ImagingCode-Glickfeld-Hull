function use_ev =remove_events_by_lever_state(ev, lever_state, t_beg, t_end, enable_state)

% use only events w/o lever press after
ev =   round(ev);
ev = ev(ev > -t_beg & ev < length(lever_state) - t_end);
use_ev = [];
for i=1:length(ev)
    if(all(lever_state(ev(i)+t_beg:ev(i)+t_end) == enable_state))
        use_ev(end+1) = ev(i);
    end
end


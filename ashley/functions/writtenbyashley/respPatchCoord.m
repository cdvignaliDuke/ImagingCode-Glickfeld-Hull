function [x,y] = respPatchCoord(resp_win,y_lim)
    x = [resp_win(1) resp_win(end) resp_win(end) resp_win(1)];
    y = [y_lim(1) y_lim(1) y_lim(2) y_lim(2)];
end
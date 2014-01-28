function offsetlabel(sel,delta,str)

if nargin < 3
    str = '';
end

xl = xlim;

labels = cellstr(cat(2,ones(length(sel),1)*str, num2str(sel')));
text(xl(2)*ones(1,length(sel)),delta*[0:length(sel)-1], labels );

return;
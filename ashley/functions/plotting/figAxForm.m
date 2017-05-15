function figAxForm(fig)
if isempty(fig)
    fig = gca;
end
fig.TickDir = 'out';
fig.Box = 'off';
axis square

end
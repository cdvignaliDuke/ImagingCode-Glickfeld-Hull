function [output] = handle0statemodel4SBML(model,tspan)
% handle0statemodel4SBML: just a really really stupid function. 
% But the SBML people want to be able to simulate models that do not
% contain any states. And the toolbox should be able to pass the SBML
% test suite cases!

[pn_xyz,pv_xyz] = SBparameters(model);
[vn_xyz,vf_xyz] = SBvariables(model);
[rn_xyz,rf_xyz] = SBreactions(model);

% assign parameters
for kloop=1:length(pn_xyz),
    eval(sprintf('%s = %g;',pn_xyz{kloop},pv_xyz(kloop)));
end
% evaluate variables
for kloop=1:length(vn_xyz),
    eval(sprintf('%s = %s;',vn_xyz{kloop},vf_xyz{kloop}));
end
% evaluate reactions
for kloop=1:length(rn_xyz),
    eval(sprintf('%s = %s;',rn_xyz{kloop},rf_xyz{kloop}));
end

% build variablevalues
variablevalues = [];
for kloop=1:length(vn_xyz),
    eval(sprintf('variablevalues(1:length(tspan),%d) = %s;',kloop,vn_xyz{kloop}));
end
% build reactionvalues
reactionvalues = [];
for kloop=1:length(rn_xyz),
    eval(sprintf('reactionvalues(1:length(tspan),%d) = %s;',kloop,rn_xyz{kloop}));
end

% build output
output.time = tspan;
output.states = {};
output.statevalues = [];
output.algebraic = {};
output.algebraicvalues = [];
output.variables = vn_xyz;
output.variablevalues = variablevalues;
output.reactions = rn_xyz;
output.reactionvalues = reactionvalues;
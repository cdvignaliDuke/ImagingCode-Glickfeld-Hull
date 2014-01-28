% This script builds the functions TranslateSBML and 
% OutputSBML for the Windows system.

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% Copyright 2008 by Henning Schmidt, henning@sbtoolbox2.org

disp('Building SBML import/export functions.');

mexcSB('TranslateSBML.c libsbml.lib -I.')
mexcSB('OutputSBML.c libsbml.lib -I.')
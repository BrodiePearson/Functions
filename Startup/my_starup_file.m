%-------------------------------------------------------------------------
% Paths
%-------------------------------------------------------------------------

% add function path
addpath(genpath('../Functions'));

% add private function path
addpath(genpath('../Private_Functions'));

%-------------------------------------------------------------------------
% MATLAB Theme
%-------------------------------------------------------------------------
schemer_import(darksteel)

%-------------------------------------------------------------------------
% Plotting Defaults
%-------------------------------------------------------------------------

% set interpreter to LATEX
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

% line properties
set(groot,'defaultLineLineWidth',3)

% colors
set(groot,'defaultFigureColor','w')

% fonts
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultFontName','New Century Schoolbook')
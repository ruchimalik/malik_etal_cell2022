function opt = mktunecurve2ddrawopt(varargin)
% MKTUNECURVE2DDRAWOPT options for plotting tuning curves
  
% take in tunecurveopt
% polar plots
opt = struct('whitebg','1',... % or 'bars', or 'stacked', or '2d'
             'colormap', 1-gray(64),... % for non-RGB maps
             'xlim', [],...
             'ylim',[],...
             'autolim', [],... % set x/ylims at a fixed distance from the track
             'datalim',[],...
             'shareddatalim', true,... % share datalims across multiple cls
             'drawcolorscale', true,... % draw colorbar or colorbox
             'drawtrackmap', true,...
             'drawspikes', true,...
             'drawarrows', true,...
             'drawtitle', true,...
             'trackmapcol', [0.5 0.5 0.5]);

opt = parseArgs(varargin, opt);
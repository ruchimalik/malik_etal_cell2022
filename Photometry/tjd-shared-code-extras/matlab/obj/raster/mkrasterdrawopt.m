function dopt = mkrasterdrawopt(varargin)
% MKRASTERDRAWOPT display options for drawraster
  
% make a sensible cmap for a rasterplot
% raster is very finicky, since it uses indexes into the cmap to draw all
% those markers pixel by pixel. The first 4 entries in its map must be as
% below, followed by an hsv colormap, so that 'parmshades' will work
%
% first four colors are reserved: black, white, 2 raster bgs (lighter and
% darker gray), i.e. the gray bands present to indicate that there's a
% cell there, even when no spiking)
  
  % color of rasters when no spikes present
  rasterplotbg_black = [0 0 0]; 
  rasterbgcolor_black = [0.15 0.15 0.15];
  rasterplotbg_white = 1-rasterplotbg_black;
  rasterbgcolor_white = 1-rasterbgcolor_black;

  % colored spike rasters
  %
  cmaprast = [rasterplotbg_black; rasterbgcolor_black; ...
              rasterplotbg_white; rasterbgcolor_white; ...
              hsv(128)];
  
  dopt = struct(...
      'cmap', cmaprast,... 
      'cmapwin', [],...
      'clim', [],...
      ...
      'drawlegend', true,...
      'drawxaxis', false,...
      'drawyaxis', false,...      ...
      'dispname', []);
      
  dopt = parseArgsLite(varargin, dopt);

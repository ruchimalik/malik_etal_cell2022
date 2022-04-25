function dopt = mksegdrawopt(varargin)
% MKSEGDRAWOPT make display structure for seg object
%
% $$$  -segdrawopt
% $$$    -colororder
% $$     -drawxaxis: show ticklabels on x-axis  
% $$$    -drawyaxis: show segment names on y-axis

  dopt = struct(...
      'colororder', [],... % (default hsv(8)/adjusted if whitebg)
      'drawxaxis', true,...
      'drawyaxis', true,...
      'segheight', 0.9,... % fractional height of each segment
      'fontsize', [],...
      'tickdir', 'out',... % inward ticks are confusing
      'ydir', 'reverse',... % so seg lists go top-to-bottom
      'drawedges', false,... % draw edges on segment rectangles
      'alpha', 1 ...
      );

  dopt = parseArgsLite(varargin,dopt);

function[occ] = occup(varargin)
% OCCUP - calculate occupany histogram
% $Id$

  a = struct(...
      'ep',[],... % e.pos struct
      'posi',[],... % pos records to use
      'params',{{}},... % pos parameters to calculate occupancy on
      'histedges',{{}}); % hist bins
  
  a = parseArgs(varargin,a);

  if ~isempty(a.posi),
    pdat = getepdi(a.ep, a.posi, a.params);
  else
    pdat = getepd(a.ep, a.params);
  end
  
  % also works for empty data arrays
  occ = histn(pdat, a.histedges, 'noedgebin');

      
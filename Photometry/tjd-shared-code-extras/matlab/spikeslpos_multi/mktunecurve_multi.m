function tc = mktunecurve_multi(varargin)
% MKTUNECURVE make tuning curves (caching of precomputed values)
%
% $$$    -tunecurve
% $$$      -tc (struct to use/fill in options from)
% $$$      -e (not in return struct)
% $$$      -edesc (not in input struct)
% $$$      -clnos : clnos to generate tuning curves for
% $$$      -tunecurveopt : tunecurve options cf. mktunecurveopt 
% $$$
% $$$      [-dispname - name to be displayed]
% $$$
% $$$      {outputs/cache}
% $$$      -tcdats : tuning curves
% $$$      -histedges : histedges for parameter histogram 
% $$$      -occ : rat occupancy in tuneparam bins
% $$$
% $$$      -cache: previous tunecurve/s, also for caching
% 
%  to force a recalc, don't include cache
% 
%  todo:
%   -test that we have all the inputs
  
  tc = struct(...
      'type', 'tunecurve',...
      ...
      'name', [],... % 'inputargs'
      'e', [],... 
      'edesc',[],...
      'clnos', [],... 
      'tunecurveopt', [],...
      ...
      'tcdats', [],... % 'outputargs' 
      'occ', [],...
      ...
      'template',[],...
      'cache', [],...
      'cache_hit', false);
  

  % parse args, including 'template' for old args
  tc = obj_reparse(tc, varargin);

  if islogical(tc.clnos),
    tc.clnos = find(tc.clnos);
  end
  
  % avoid dups; also sorts clnos
  tc.clnos = unique(tc.clnos);

  % reuse tunecurves from cache if possible
  tc = obj_cachesearch(tc);

  % return cache hit if found

  if tc.cache_hit,
    tc = obj_cleanup(tc);
    return
  end
  
  % we didn't get a cache hit, recalc
  
  % can't recalc without the 'e' struct
  if isempty(tc.e),
    error(['No matching object found in cache; not enough data to recalculate,' ...
           'try passing in ''e''.']);
  end

  % default: make tc for all clusters
  if isempty(tc.clnos),
    tc.clnos = 1:length(tc.e.epoch(tc.tunecurveopt.runepochi).cl);
  end
  
  % we have to recalc the tuning curves
  [tc.tcdats tc.occ] = dotunecurve_multi(...
      'e', tc.e,...
      'clnos', tc.clnos,...
      'opt', tc.tunecurveopt);
  
  tc = obj_cleanup(tc);
  return;
  

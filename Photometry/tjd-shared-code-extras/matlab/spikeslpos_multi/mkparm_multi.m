function parm = mkparm_multi(varargin)
% MKPARM mk a parm struct
  
% use the first param in a tuning curve for data. if there is a second
% tuning param, use it for marker color. 
  
  parm =  struct(...
      'type', 'parm',...
      ...
      'e', [],... %inputs
      'epochi', [],...
      'timewin', [],...
      'params', [],...
      'tunecurveopt', [],... % we can extract tuneparms from a tunecurve
      ...
      'tpfdat', [],... %outputs
      ...
      'template', [],...
      'cache', [],...
      'cache_hit', false);
      
  parm = obj_reparse(parm, varargin);

  parm = obj_cachesearch(parm);
  
  if parm.cache_hit
    parm = obj_cleanup(parm);
    return;
  end

  % use params from tunecurveopt, if provided
  if ~isempty(parm.tunecurveopt),
    params = parm.tunecurveopt.tuneparams;
  else
    params = parm.params;
  end
  
  % Get animal position data (all params used in tuning curve)
  parm.tpfdat=getepd(parm.e.epoch(parm.epochi).pos, ...
                     'time', ...
                     params{:});

  % select time window of interest
  if ~isempty(parm.timewin), % default use all
    % get pos records in xstart-xend window
    % (so we can plot pos outside of where there are spikes)
    parm.tpfdat = parm.tpfdat(parm.timewin(1) < parm.tpfdat(:,1) & ...
                              parm.tpfdat(:,1) < parm.timewin(2),:);
  end
  
  parm = obj_cleanup(parm);
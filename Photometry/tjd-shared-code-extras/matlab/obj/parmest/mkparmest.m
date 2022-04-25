function [pe tc clh] = mkparmest(varargin)
% MKPARMEST make a reconstruction (PARaMeter ESTimation) based on spiking data
%
% $$$     -inputs:
% $$$      -e : experimental structure, with clusters
% $$$      -edesc : 'desc' field from 'e'
% $$$      -clnos : cluster numbers to use in estimate
% $$$      -clperm : permutation of cluster->tunecurve mapping
% $$$      -timewin : timewindow to make estimate for
% $$$      -timebinsize : size of timebin to use for estimation
% $$$      -overlapsperbin: sliding window estimate
% $$$      -timebincenter: fit timebins inside timewin, centered
% $$$      -parmestopt : options structure for parmest
% $$$      -tunecurveopt : options structure for tunecurves
% $$$      -tunecurve : pre-computed tunecurve object
% $$$
% $$$      {outputs/cache}
% $$$      -map : posest vs time map 
% $$$      -timehistedgecenters: centers of first/last time bins
% $$$      -parmhistedgecenters centers of first/last param bins
% $$$      -modei: index of peak bin in each column
% $$$      -mode: value of map at each modei
% $$$      -spikecount: # of spikes contributing to each column
% $$$      -cellcount: # of unique cells contributing to each column
% $$$      -parm: actual value of reconstructed parm, if available
% $$$      -label: label for plotting (from tunecurve)
% $$$
% $$$     optional data to use for caching:
% $$$      -cache : previous parmest/s, to use for caching
% 
%  to force a recalc, don't include cache.
% 
%  todo:
%   -test that we have all the inputs
  
  % initialize outputs
  [pe tc clh] = deal([]); %#ok
  
  pe = struct(...
      'type', 'parmest',...
      ...
      'e', [],...
      'edesc',[],...
      'clnos', [],... 
      'clperm', [],... 
      'timewin',[],...
      'timebinsize', [],... 
      'overlapsperbin', [],...
      'timebincenter', false,...
      'histbins', [],...
      'firstspiketriggers', [],...
      'parmestopt', [],... 
      'tunecurveopt', [],...
      'tunecurve', [],... % if you are sure you know which tunecurve to use
      ...
      'map', [],... 
      'timehistedgecenters', [],...
      'parmhistedgecenters', [],...
      'mode',[],...
      'modei',[],...
      'mode_parm', [],...
      'mode_marg', [],...
      'modei_marg', [],...
      'mode_parm_marg', [],...
      'spikecount', [],...
      'cellcount', [],...
      'posi', [],...
      'parm', [],...
      'parmi', [],...
      'label', [],...
      ...
      'template',[],...
      'cache', [],...
      'cache_hit', false);
      
      %      'parmwini', [],... % for efficiency in, e.g. thetamap.m
      
  % handle 'template' object arg
  pe = obj_reparse(pe,varargin);    
  
  if islogical(pe.clnos)
    pe.clnos = find(pe.clnos);
  end
  
  if ~isempty(pe.tunecurve),
    
    if ~isempty(pe.tunecurveopt)
      error('can''t provide a tunecurve and a tunecurveopt');
    else
      pe.tunecurveopt = pe.tunecurve.tunecurveopt;
    end
    
    if ~isempty(pe.clnos)
      error('can''t provide a tunecurve and ''clnos''');
    else
      pe.clnos = pe.tunecurve.clnos;
    end
    
    tc = pe.tunecurve;
    
    pe.tunecurve = []; % (it's available in the caller; save some memory)
    
  end
  
  % clnos = [] means use all clusters (do before cachesearch);
  if isempty(pe.clnos),
    pe.clnos = 1:length(pe.e.cl);
  end
  
  % search cache for appropriate object
  pe = obj_cachesearch(pe);
  
  if pe.cache_hit
    pe = obj_cleanup(pe);
    return;
  end
   
  
  %%% no cache hit, make new parmest maps

  if isempty(pe.e),
    error(['No matching object found in cache; not enough data to recalculate,' ...
           'try passing in ''e''.']);
  end
  
  % clnos = [] means use all clusters
  if isempty(pe.clnos),
    clnos = 1:length(pe.e.cl);
  else
    clnos = pe.clnos;
  end
  
  if isempty(tc)
    % note: tc returned may be superset of requested clnos
    tc = mktunecurve('e',pe.e,...
                     'clnos',clnos,...
                     'tunecurveopt', pe.tunecurveopt,...
                     'cache', pe.cache);

    % select appropriate tuning curves for e.cl(clnos)
    [dummy tcdatsi] = ismember(clnos,tc.clnos); %#ok
    colons = repmat({':'},1,ndims(tc.tcdats)-1);
    tcdats = tc.tcdats(tcdatsi,colons{:});
  else
    tcdats = tc.tcdats;
  end
  
  ntcdims = ndims(tc.tcdats)-1;
  
  if ~isfield(pe.parmestopt, 'decim_type')
    pe.parmestopt.decim_type = [];
  end
  
  % get spiking histograms for e.cl(clnos)
  clh = mkclhist(...
      'e', pe.e,...
      'clnos', clnos,...
      'histbins', pe.histbins,...
      'timewin', pe.timewin,...
      'timebinsize', pe.timebinsize,...
      'overlapsperbin', pe.overlapsperbin,...
      'timebincenter', pe.timebincenter,...
      'firstspiketriggers', pe.firstspiketriggers,...
      'binarize', pe.parmestopt.binarizect,...
      'decim_frac', pe.parmestopt.decim_frac,...
      'decim_type', pe.parmestopt.decim_type,...
      'rate', strcmp(pe.parmestopt.method,'basis'),...
      'cache', pe.cache);

  pe.spikecount = clh.spikecount;
  pe.cellcount = clh.cellcount;
  
% $$$   % select particular rows of the tuning curve to use
% $$$   if ~isempty(pe.parmwini),
% $$$     tcdats = tcdats(:,pe.parmwini(1):pe.parmwini(2));
% $$$   end
  
  % create parmest map
  pe.map = doparmest(...
      'clh', clh,...
      'tcdats', tcdats,...
      'occ', tc.occ,...
      'clperm', pe.clperm,...
      'opt', pe.parmestopt);
  

  if isempty(pe.map),
    pe.label = tc.tunecurveopt.name;
    pe = obj_cleanup(pe);
    return
  end
  
  % get coordinates of map edges (for imagesc)

  pe.timehistedgecenters = clh.histedgecenters;
  
  for dim = 1:ntcdims
    parmhistedges = tc.tunecurveopt.histedges{dim};
    parmctrs{dim} = ctrs(parmhistedges); %#ok
    if ~isempty(tc.tunecurveopt.keepbins) && ...
          ~isempty(tc.tunecurveopt.keepbins{dim}),
      parmctrs{dim} = parmctrs{dim}(tc.tunecurveopt.keepbins{dim});
    end
    % save ends for export (saves some memory, I guess?)
    pe.parmhistedgecenters{dim} = parmctrs{dim}([1 end]);
  end
  
  %%% get n-dimensional 'mode' (max of pdf for each estimate)
  [pe.mode pe.modei] = max(pe.map(:,:),[],2);
  mapsz = size(pe.map);
  mapndims = ndims(pe.map);
  mapsubs = cell(1,mapndims-1); %#ok --needs to be defined so ind2sub
                                %knows how many outputs to calc
  [mapsubs{:}] = ind2sub(mapsz(2:end),pe.modei);
  pe.modei = cell2mat(mapsubs);
  
  %%% get modes of marginals

  if ntcdims == 1,
    pe.mode_marg = pe.mode;
    pe.modei_marg = pe.modei;
  
  elseif ntcdims == 2,
    for dim = 1:ntcdims
      mapsum = reshape(sum(pe.map,4-dim), size(pe.map,1),[]);
      [pe.mode_marg(:,dim) pe.modei_marg(:,dim)] = max(mapsum,[],2);
    end
  end
      
  %%% get point estimates in parameter space

  for dim = 1:ntcdims

    % the n-D joint point estimate
    pe.mode_parm(:,dim) = parmctrs{dim}(pe.modei(:,dim));

    % the independent point estimates based on the mode of the marginals
    pe.mode_parm_marg(:,dim) = parmctrs{dim}(pe.modei_marg(:,dim));
  
  end
    
  %%% get 'parm'--mean actual value of parameter, if available
  
  % times of all pos records
  pos_t = getepd(pe.e.pos,'time');

  % time of middle of each parmest bin
  if ~isempty(pe.timebinsize),
    if isempty(pe.overlapsperbin)
      pe_t = pe.timehistedgecenters(1):pe.timebinsize:pe.timehistedgecenters(2);
    else
      pe_t = pe.timehistedgecenters(1):...
             pe.timebinsize/pe.overlapsperbin:...
             pe.timehistedgecenters(2);
    end
  else
    pe_t = mean(pe.histbins,2)';
  end

  
  % get tuneparam at closest pos record to center
  pe.posi = findnearest(pos_t, pe_t', 'sorted');
  pe.parm = getepdi(pe.e.pos, pe.posi, pe.tunecurveopt.tuneparams{:});
  
  % get corresponding bin indexes for actual parms (useful for confusion
  % matrix, e.g.)
  
  for dim = 1:ntcdims
    [dum dimparmi] = histc(pe.parm(:,dim), pe.tunecurveopt.histedges{dim});
    if ~isempty(pe.tunecurveopt.keepbins) &&...
                ~isempty(pe.tunecurveopt.keepbins{dim})
      keepbins = pe.tunecurveopt.keepbins{dim};
      % we assign values corresponding to non-kept bins to 0 (just like out of
      % range values), and adjust index of later bins so that they
      % index properly
      newbins = (1:length(keepbins)) + -cumsum(~keepbins);
      newbins(~keepbins) = 0;
      dimparmi = newbins(dimparmi);
    end
    pe.parmi(:,dim) = dimparmi;
  end
  
  %%% create label for plotting from tunecurve name
  pe.label = tc.tunecurveopt.name;
  
  pe = obj_cleanup(pe);
  

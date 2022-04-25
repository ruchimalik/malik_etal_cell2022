function clh = mkclhist(varargin)
% MKCLHIST get spike firing histograms for clusters (with caching)

  clh = struct(...
      'type', 'clhist',...
       ...
      'e',[],... % inputs
      'edesc',[],...
      'clnos', [],... 
      'timewin',[],...
      'timebinsize',[],...
      'overlapsperbin', [],...
      'timebincenter', false,... % do not extend past timewin, center
          ...                    % remaining bins
      'histbins', [],...
      'binarize', false,...
      'decim_frac', [],... % resample # of spikes
      'decim_type', [],... % resample # of spikes
      'rate', false,...
      'firstspiketriggers', [],... % only report first spike after each time in list
      ...
      'clhists', [],... % outputs
      'histbinsize', [],...
      'nbins', [],...
      'histedgecenters', [],...
      'spikecount', [],...
      'cellcount', [],...
      ...
      'template',[],...
      'cache',[],...
      'cache_hit',false);
  
      
  % handle 'template' object arg
  clh = obj_reparse(clh,varargin);
  
  % search cache for appropriate object
  clh = obj_cachesearch(clh);
  
  if clh.cache_hit
    clh = obj_cleanup(clh);
    return;
  end

  
  %%%% no cache hit, make a new clhist
  cls = clh.e.cl(clh.clnos);
  ncls = length(cls);
  
  if ~isempty(clh.decim_frac) && ...
        strcmp(clh.decim_type, 'dropspikes') && ...
        (clh.decim_frac>1 || clh.decim_frac<0); 
    error('''decim_frac'' must be on [0,1] for ''dropspikes'' decimation'); 
  end
      
  
  if ~isempty(clh.overlapsperbin) && clh.overlapsperbin ~= 1;
    overlap = true;
    if clh.overlapsperbin - fix(clh.overlapsperbin) ~= 0
      error('''overlapsperbin'' must be an integer');
    end
  else
    overlap = false;
  end
    
  if ~isempty(clh.histbins)
    % histbins provided
    
    if any([clh.timewin clh.timebinsize clh.overlapsperbin]),
      error('''can''t provide ''histbins'' and ''timewin''/''timebinsize''/''overlapsperbin''');
    end
    
    if size(clh.histbins,2) ~= 2,
      error ('''histbins'' must be mx2 array of times');
    end
    
    % HACK ALERT
    % we treat inter-bin times as extra histogram bins, then throw them
    % away before returning.
    
    % the *eventual* number of bins
    clh.nbins = size(clh.histbins,1); % the # of *pre-reshape* rows
    histbins = reshape(clh.histbins',1,[]);
    
    % bins can be of different sizes, use to divide later to get rates. 
    clh.histbinsize = diff(clh.histbins,[],2)';
    
  else
    % calculate histbins
    if overlap, 
      % Fudge with multiple of bins, then combine neighboring bins later!
      timebinsize = clh.timebinsize/clh.overlapsperbin;
    else
      timebinsize = clh.timebinsize;
    end
      
    % handle special case of timebinsize = tend-tstart 
    if clh.timewin(2) - clh.timewin(1) == timebinsize,
      histbins = clh.timewin;
    else
      if ~clh.timebincenter
        % make sure we have at least one bin, starting at tstart, and that all
        % bins are the same size (note this means we often analyze data past
        % time tend).
        histbins = clh.timewin(1):timebinsize:clh.timewin(2)+timebinsize;
      else
        % do not extend past ends of timebin, center remaining bins
        nbins = floor(diff(clh.timewin)./timebinsize);
        binsctr = mean(clh.timewin);
        binsdur = nbins*timebinsize;
        histbins = (binsctr-binsdur/2):timebinsize:(binsctr+binsdur/2);
      end 
    end
    
    clh.nbins = length(histbins)-1;
    clh.histbinsize = clh.timebinsize; % report requested binsize, even for overlaps!
  end

  if clh.nbins == 0
    clh.histedgecenters = [];
  elseif ~isempty(clh.timebinsize)
    clh.histedgecenters = [histbins(1) + clh.timebinsize./2 ...
                        histbins(end) - clh.timebinsize./2];
  else
    clh.histedgecenters = [mean(histbins(1:2)) mean(histbins(end-1:end))];
  end

  % build rate histograms and scale as requested
  tstart = min(histbins(:));
  tend = max(histbins(:));

  % initialize histograms, counters
  clh.clhists = zeros(ncls,clh.nbins);

  clh.spikecount = zeros(1,clh.nbins);
  clh.cellcount = clh.spikecount;

  if ncls == 0 || clh.nbins == 0,
    return
  end
  
  % same across all cls
  ti = getepindex(cls(1),'time');
  
  for j = 1:ncls, 
    
% $$$     % get all spike times 
% $$$     clts = getepd(cls(j), 'time');
    
    % get spikes in tstart-tend window 
    clts = cls(j).dat(tstart < cls(j).dat(:,ti) &...
                      cls(j).dat(:,ti) < tend);
    
    %%% filter for first spike after each time in vector clh.firstspiketriggers
    if ~isempty(clh.firstspiketriggers),
      % get all trigger times in tstart-tend window
      triggers = clh.firstspiketriggers(tstart < clh.firstspiketriggers & ...
                                      clh.firstspiketriggers < tend);

      % if no valid trigger, can be no first spike after trigger in time window
      if isempty(triggers),
        clts = [];
      elseif ~isempty(clts),
        % get times of nearest spikes to triggers
        nearesti = findnearest(clts, triggers, 'sorted');
        nearestt = clts(nearesti);
        % findnearest can find spike before trigger, in this case get next spike
        addonei = nearestt < triggers;
        nearesti(addonei) = nearesti(addonei)+1;
        % if a spike is the first after several triggers, only count it once
        nearesti = unique(nearesti);
        % there may only be one spike in the list
        nearesti = nearesti(nearesti <=length(clts));
        clts = clts(nearesti);
      end
    end
    
    % get firing histogram over time
    if any(clts) % only increment if spikes

      % if overlapping bins (possible), then use a big for-loop hack,
      % really slow--use 'overlapsperbin' if possible.
      if ~any(diff(histbins)<0)
        % normal case, no overlap
        clhist = histc(clts',histbins);
        clhist = clhist(1:end-1); %last bin meaningless
      else
        clhist = [];
        % iterate over and concatenate the non-overlapping parts
        ends = find(diff(histbins)<0);
        starts = [1 ends+1];
        ends = [ends length(histbins)];
        for k = 1:length(starts)
          clhist = [clhist histc(clts', histbins(starts(k):ends(k)))];
        end
      end
      
      % HACK ALERT - throw away the nonsense 'spacer' clhists, see above
      if ~isempty(clh.histbins),
        clhist = clhist(1:2:end);
      end

      % add to output array
      clh.clhists(j,:) = clhist;
      
    end
  end
  
  if ~isempty(clh.decim_frac)
    
    switch(clh.decim_type)
      case('dropspikes')
       
       % randomly delete spikes according to the fraction (0,1) requested

       % can't just make random draws from the pdf for spikes to delete, since we
       % could ask to remove more spikes than exist for a particular
       % cluster. Need to draw randomly from existing spikes.

       % how many spikes total do we need to get rid of, per bin?
       % (jitter is to avoid rounding bias with integers)
       del_nspikes = round(jitter(sum(clh.clhists)*(1-clh.decim_frac),1e-6));

       for k = 1:clh.nbins
         if del_nspikes(k) > 0;
           
           % create list of spikes (identified by cell #) 
           spikelist = []; %#ok
           for j = 1:ncls
             spikelist = [spikelist repmat(j,1,clh.clhists(j,k))];
           end
           
           % randomly delete requested # of spikes
           spikelist = spikelist(randperm(numel(spikelist)));
           spikelist(1:del_nspikes(k)) = [];

           % recreate hist array from remaining spikes
           if isempty(spikelist)
             clh.clhists(:,k) = zeros(ncls,1);
           else
             clh.clhists(:,k) = accumarray(spikelist',...
                                           ones(size(spikelist)), ...
                                           [ncls 1]);
           end
         end
         
       end
    
     case 'scalerate'
      % just scale # of spikes by requested fraction
      % results in non-integer spike #s
      clh.clhists = clh.clhists .* clh.decim_frac;
    
     otherwise
      error(['decim_frac requsted, but not decim_type specified in ' ...
             'parmestopt']);
      
    end
  
  end

  % # of spikes per time bin
  clh.spikecount = sum(clh.clhists);
  % # of cells contributing
  clh.cellcount = sum(clh.clhists>0);

  if overlap
    % we run a box filter over the histogram made up of several smaller
    % time bins
    box = ones(1,clh.overlapsperbin);
    clh.clhists = conv2(clh.clhists, box, 'valid');
    clh.spikecount = conv2(clh.spikecount, box, 'valid');
    clh.cellcount = conv2(clh.cellcount, box, 'valid');
    clh.nbins = size(clh.clhists,2);
  end
  
  if clh.rate
    % Use firing rate rather than spike count for basis method, helps
    % keep estimate stable over large range of time bins
    clh.clhists = clh.clhists ./ clh.histbinsize;
    clh.spikecount = clh.spikecount ./ clh.histbinsize;
  end
  
  if clh.binarize,
    % only count *whether* a given cell fired (used for rasters)
    clh.clhists = clh.clhists ~= 0;
  end
        
  clh = obj_cleanup(clh);
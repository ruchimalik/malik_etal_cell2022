function hs = mkhistedges(varargin)
% MKHISTEDGES several ways to specify histogram edges
  
  edges  = struct(...
      'edges',[],...
      'nbins',[],...
      'binsize',[],...
      'range',[]);
  
  hs = parseArgs(varargin,hs);
  
  
  % histogram edge specification 
  %
  % 'edges' will get passed directly to histn calls (see 'doc histc'),
  % allows for non-uniform binning, greatest control. No other options
  % may be provided.
  %
  % 'range' specifies the range of data to be binned ([m x 2])
  % 
  % only one of 'nbins'/'binsize' can be provided. If neither is
  % provided, default to histnbins = 50
  %
  % 'binsize': size of bins. If not evenly divisible by length of
  % requested range, we add a bin at the end to capture the full
  % range. (use nbins instead!)
  %
  % 'nbins' number of bins to use. Usually preferred as we don't have the
  % problem of the last bin going off beyond the requested range.
  %
  % 'keepbins', when using 'edges' or 'nbins', specify which bins to keep
  % and in which order when making tunecurve (useful for discarding the
  % low-velocity slice of position reconstructions, e.g.). logical index
  
  nbins = [];
  
  if ~isnumeric(hs.edges) ||...
        ~isnumeric(hs.range) ||...
        ~isnumeric(hs.nbins) ||...
        ~isnumeric(hs.binsize),
    error('args must be numeric');
  end
  
  if ~isempty(hs.edges),
    if(~isempty(hs.range) || ...
       ~isempty(hs.nbins) || ...
       ~isempty(hs.binsize)),
      error('if ''edges'' is provided, no other args may be provided');
    end
    nbins = length(hs.edges)-1;
  else
    
    % no edges provided, make sure we have 1 of nbins/binsize
    switch sum([~isempty(hs.nbins) ~isempty(hs.binsize)]),
     case 0,
      hs.nbins = 50; % default of 50 bins, linearly spaced
     case 1,
      % OK!
     case 2,
      error('only one of ''nbins'' and ''binsize'' may be provided');
    end
    if ~isempty(hs.nbins),
      nbins = hs.nbins;
    end
  end
    
  
  if ~isempty(a.opt.histspecs.edges), % done like dinner!
    histedges = a.opt.histspecs.edges; %#ok
    
  else % we need to create nbins or binsize bins in given range
    if isempty(a.opt.histspecs.range),
        histrange = a.opt.histspecs.range; %#ok
      else
        
      end
      
      if ~isempty(a.opt.histspecs.nbins),
        histedges = linspace(histrange(1), histrange(2), a.opt.histspecs.nbins+1);
      else % use binsize
        histedges = histrange(1):a.opt.histspecs.binsize:histrange(2)+a.opt.histspecs.binsize;
        if histedges(end-1) >= histrange(2),
          histedges(end) = [];
        end
      end
      
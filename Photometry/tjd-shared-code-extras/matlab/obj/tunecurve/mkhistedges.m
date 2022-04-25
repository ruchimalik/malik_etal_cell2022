function ed = mkhistedges(varargin)
% MKHISTEDGES several ways to specify histogram edges
  
  a  = struct(...
      'edges',[],...
      'nbins',[],...
      'binsize',[],...
      'range',[],...
      'log10range',[]);
  
  a = parseArgsLite(varargin,a);
  
  
  % histogram edge specification 
  %
  % 'edges' will get passed directly to histn calls (see 'doc histc'),
  % allows for non-uniform binning, greatest control. No other options
  % may be provided.
  %
  % 'range' specifies the range of data to be binned ([m x 2])
  %
  % 'log10range' specifies the range to be binned in powers of 10 --
  % e.g. logrange = [-3 0] creates bins from 0.001 to 1. not compatible
  % with 'binsize' argument.
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
  
  if ~isnumeric(a.edges) ||...
        ~isnumeric(a.range) ||...
        ~isnumeric(a.log10range) ||...
        ~isnumeric(a.nbins) ||...
        ~isnumeric(a.binsize),
    error('args must be numeric');
  end
  
  if ~isempty(a.edges),
    if(~isempty(a.range) || ...
       ~isempty(a.log10range) || ...       
       ~isempty(a.nbins) || ...
       ~isempty(a.binsize)),
      error('if ''edges'' is provided, no other args may be provided');
    else
      ed = a.edges;
      return;
    end
  end

  % no edges provided
  if isempty(a.range) && isempty(a.log10range),
    error('''range'' or ''log10range'' must be provided if edges are not specified');
  end
  
  switch sum([~isempty(a.nbins) ~isempty(a.binsize)]),
   case 0,
    a.nbins = 50; % default of 50 bins, linearly spaced
   case 1,
    % OK!
   case 2,
    error('only one of ''nbins'' and ''binsize'' may be provided');
  end

  if ~isempty(a.nbins),
    if ~isempty(a.log10range),
      ed = logspace(a.log10range(1), a.log10range(2), a.nbins+1);
    else
      ed = linspace(a.range(1), a.range(2), a.nbins+1);
    end
  else 
    % use binsize: last bin *includes* end of range (i.e. may go past
    % end of range)
    if ~isempty(a.log10range),
      error('Can''t provide ''log10range'' and ''binsize'' -- use ''nbins''');
    end
    ed = a.range(1):a.binsize:a.range(2)+a.binsize;
    if ed (end-1) >= a.range(2),
      ed(end) = [];
    end
  end
  
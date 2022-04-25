function seg = mkseg(varargin)
% MKSEG make a seg object to be plotted by drawseg
  
  seg = struct(...
      'type', 'seg',...
      ...
      'seglists', {{}},... % list of segments to plot
      'segstruct', [],...
      'segfields', {{}},...
      'segnames', {{}},... % names of segments
      'timewin', [],... % time window (use all)
      ...
      'template', [],...
      'cache', [],...
      'cache_hit', false);

  seg = obj_reparse(seg, varargin);

  if ~isempty(seg.segstruct),
    if ~isempty(seg.seglists)
      error('can''t provide segslists and segstruct');
    end
    if isempty(seg.segfields), 
      usefields = fieldnames(seg.segstruct);
    else
      usefields = seg.segfields;
    end
    for fn = 1:numel(usefields)
      seg.seglists{fn} = seg.segstruct.(usefields{fn});
    end
  end
  
  % Don't bother doing searching cache, just pass it back out if it was passed
  % in

  if ~isempty(seg.segnames),
    if numel(seg.segnames) ~= numel(seg.seglists),
      error(['If segnames are provided, must be same number of names as ' ...
             'seglists']);
    end
  end
  

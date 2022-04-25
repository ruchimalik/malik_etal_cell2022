function hs = drawseg(varargin)
% DRAWSEG plot a seg object, using prefs from segopt, segdrawopt

  hs = [];
  
  a = struct('seg', [],...
             'opt', [],...
             'plottypeno', [],...
             'whitebg', [],...
             'ax', [],...
             'dpi', [],...
             'xlim', []);
  
  
  a = parseArgsLite(varargin,a);
  
  %%% setup
  
  if isempty(a.opt),
    a.opt = mksegdrawopt();
  end
  
  if isempty(a.ax),
    a.ax = gca;
  end

  % patches drawn by seg_plot require zbuffer
  figh = ancestor(a.ax, 'Figure');
  if strcmpi(get(figh, 'renderer'), 'zbuffer');
    warning('drawsegs requires zbuffer renderer');
    set(figh, 'Renderer', 'Zbuffer');
  end
  
  % get axis size
  oldunits = get(a.ax,'units');
  set(a.ax,'units','pixels');
  axpospix = get(a.ax,'position');
  set(a.ax,'units', oldunits);
  
  % set 'linear' (non-'log') scaling
  set(a.ax, 'yscale', 'linear');

  % choose a default colororder; rotate it depending on plot number
  if isempty(a.opt.colororder),
    a.opt.colororder = hsv(8);
    if ~isempty(a.plottypeno),
      a.opt.colororder= circshift(a.opt.colororder,-a.plottypeno+1);
    end
    if a.whitebg,
      a.opt.colororder= rgb2hsv(a.opt.colororder);
      a.opt.colororder(:,3) = 0.7;
      a.opt.colororder= hsv2rgb(a.opt.colororder);
    end
  end
  
  
  %%% draw segs
  
  %seg_plot, from fkSegments, takes care of finding segs within range, etc
  
  nlists = numel(a.seg.seglists);
  ht = a.opt.segheight;
  yoff = (1:nlists)- (ht/2); % center segs on integers

  % get enough colors for seg_plot, repeat if necessary
  color = repmat(a.opt.colororder, ceil(nlists/size(a.opt.colororder,1)),1);
  color = color(1:nlists,:);

  if a.opt.drawedges,
    edgecolor = color;
  else
    edgecolor = 'none';
  end
  
  seg_plot(...
      a.seg.seglists,...
      'FaceColor', color,...
      'EdgeColor', edgecolor,...
      'Axis', a.ax,...
      'Xlim', a.seg.timewin,...
      'YOffset', yoff,...
      'Alpha', a.opt.alpha,...
      'Height', ht);
  
  set(a.ax, 'YTick', 1:nlists);
  if ~isempty(a.seg.segnames)
    set(a.ax, 'YTickLabel', a.seg.segnames);
    if ~isempty(a.opt.fontsize)
      set(a.ax, 'FontSize', a.opt.fontsize);
    end
  end
  ylim(a.ax, [0 nlists] + 1/2);

  set(a.ax, 'tickdir', a.opt.tickdir);
  
  % this way segs are listed top to bottom
  set(a.ax, 'ydir', a.opt.ydir)
  
  %%% wrap-up
  
  %%% set xlim (usu according to requested timewin_plot)
  if ~isempty(a.xlim),
    xlim(a.ax, a.xlim),
  elseif ~isempty(a.seg.timewin)
    xlim(a.seg.timewin);
  end
  
  % remove numbers from x/y axes if requested. Leave ticks
  if ~a.opt.drawxaxis,
    set(a.ax,'xticklabel', []);
  end

  if ~a.opt.drawyaxis,
    set(a.ax,'yticklabel', []);
  end
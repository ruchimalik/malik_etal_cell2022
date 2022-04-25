function hs = drawparmest(varargin)
% DRAWPARMEST plot a parmest object, using parmestdrawopt

  hs = [];

  a = struct(...
      'parmest', [],...
      'opt', [],...
      'parm', [],...
      'parmdrawopt', [],...
      'whitebg',[],...
      'ax', [],...
      'xlim', []);
  
  a = parseArgsLite(varargin,a);
  
  if isempty(a.ax),
    a.ax = gca;
  end

  if isempty(a.opt),
    a.opt = mkparmestdrawopt();
  end

  % set/get the colormap
  if ~isempty(a.opt.cmap),
    cmap = colormap(a.ax, a.opt.cmap);
  else
    cmap = colormap(a.ax);
  end
  
  if isempty(a.parmdrawopt),
    a.parmdrawopt = mkparmdrawopt();
  end

  % get map from parmest
  map = a.parmest.map;

  if ndims(map) > 3,
    warning('can''t draw parmests > 3 dimensions--ignoring higher dims');
    ones_idx = repmat({1},1,ndims(map) - 3);
    map = map(:,:,:,ones_idx{:});
  end
  
  % the output from parmest has time in first dim, then the n dimensions
  % of the tuning curve in subsequent dimensions. When drawing, we will
  % put time on the x-axis, the first tc parameter on the y-axis, and map
  % the second tc parameter to colors (using draw2d).

  %map = reshape(shiftdim(map,1),size(map,2),size(map,1),size(map,3));
  map = permute(map,[2 1 3]);
  
  % if requested, flatten map to marginal of first estimated parm
  if a.opt.marginalize,
    map = sum(map,3);
  end
  
  [nrows ncols nslices] = size(map);

  % update dispname
  if isempty(a.opt.dispname)
    if ~isempty(a.parm),
      a.opt.dispname = [a.parm.tunecurveopt.name '_est'];
    else
      a.opt.dispname = 'parmest';
    end      
    if a.opt.modemap,
      a.opt.dispname = [a.opt.dispname '_mode'];
      if a.opt.modethresh > 0;
        a.opt.dispname = [a.opt.dispname '_th'];
      end
    end
  end
  
  if a.opt.modemap,
    
    th = a.parmest.mode > a.opt.modethresh;
    modei = a.parmest.modei(th,:);

    % build array of indexes for accumarray
    ind = [modei(:,1) find(th)];
    arrsz = [nrows ncols];
    
    % add in 3d dimension, if needed
    if size(modei,2) > 1,
      ind = [ind modei(:,2)];
      arrsz = [nrows ncols 2];
    end
    
    % use accumarray to make a (logical) image to plot
    map = logical(accumarray(ind, 1, ...
                             arrsz,@sum,0));
    % map = sparse(a.parmest.modei(:,1), 1:ncols, a.parmest.modei(:,2), nrows, ncols);
    %    map = sparse(a.parmest.modei(th), find(th), a.parmest.mode(th), nrows, ncols);
    % th = a.parmest.mode > a.opt.modethresh;
    % map = sparse(a.parmest.modei(th), find(th), 1, nrows, ncols);
  end
  
  if a.opt.binarize,
    map = map>0;
  end 
  
  % sharpen/scale map, if requested
  if ~isempty(a.opt.autoscale),
    
    % normalize maps by requested percentile *of non-zero elements*
    normf = prctile(map(map ~= 0), ...
                    a.opt.autoscale*100);
    if normf == 0 || isnan(normf);
      % makes each non-zero value a very large number when we divide
      % by it; 'binarizes' the map
      normf = realmin; 
    end
    map = map./normf;
  end
  
  if ~isempty(a.opt.sharpen)
    map(:) = map(:).^a.opt.sharpen;
  end
  
  if ~isempty(a.opt.scale)
    map(:) = map(:).*a.opt.scale;
  end

  %%% draw map data
  
  if ~islogical(map), 
    map = real(map); 
  end
  
  h = draw2d('data', map,...
             'ax', a.ax,...
             'RGBchan', a.opt.RGBchan,...
             'cmap', cmap,...
             'cmapwin', a.opt.cmapwin,...
             'cmapwin_whitebg', a.opt.cmapwin_whitebg,...
             'datalims', [0 1],...
             'zcolnorm', false,... % renormalize columns to sum to 1 after scale/crop
             'xhistendctrs', a.parmest.timehistedgecenters,...
             'yhistendctrs', a.parmest.parmhistedgecenters{1},...
             'whitebg', a.whitebg,...
             'shading', 'flat');
  %             'zlims', [0.1 0.5; 0.1 0.5; 0.1 0.5]);

  hs = [hs; h];     
  
  %%% plot behavior parameters
  if ~isempty(a.parm)
    a.parmdrawopt.legend = []; % we'll draw our own here
    h =  drawparm('parm', a.parm,...
                  'opt', a.parmdrawopt,...
                  'whitebg', a.whitebg,...
                  'ax', a.ax);
    hs = [hs; h];
  end
  
  %%% draw legend
  if a.opt.drawlegend,
    h = contdrawlegend('chanlabels', a.parmest.label,...
                       'dispname', a.opt.dispname,...
                       'plottype', 'image',...
                       'cmaptop', a.parmdrawopt.ptcolor(1,:),...
                       'whitebg', a.whitebg,...
                     'ax', a.ax);
    hs = [hs; h];
  end
  
  
  %%% draw scale bar
  if a.opt.drawscalebar,
    h = scalebar('ax', a.ax,...
                 'barlen', a.parmest.timebinsize,...
                 'bartxt', ['binsize = ' timestringfmt(a.parmest.timebinsize)],...
                 'xl', a.xlim,...
                 'whitebg', a.whitebg);
% $$$                'color', [0.7 0.7 0.7]);
    hs = [hs; h];
  end
    
  %%% set xlim (usu according to requested timewin_plot)
  if ~isempty(a.xlim),
    xlim(a.ax, a.xlim),
  end
  
  lims = objbounds(hs);

  if ~isempty(a.parmdrawopt.datalim) && numel(a.parmdrawopt.datalim)>1
    yl = [min([a.parmdrawopt.datalim(1) lims(3)]),...
          max([a.parmdrawopt.datalim(2) lims(4)])];
  else 
    yl = lims(3:4);
  end
  ylim(a.ax, yl);
  
  % remove numbers from x/y axes if requested. Leave ticks
  if ~a.opt.drawxaxis,
    set(a.ax,'xticklabel', []);
  end
  
  if ~a.opt.drawyaxis,
    set(a.ax,'yticklabel', []);
  end
  
  
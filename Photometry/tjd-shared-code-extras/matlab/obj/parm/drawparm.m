function hs = drawparm(varargin)
% DRAWPARM draw points from an e.pos parameter
  
  hs = [];
  
  a =  struct(...
      'parm', [],...
      'opt', [],...
      'whitebg', [],...
      'ax', [],...
      'xlim', []);
  
  a = parseArgsLite(varargin, a);
  
  if a.opt.draw,
    
    % get axis size
    oldunits = get(a.ax,'units');
    set(a.ax,'units','pixels');
    axpospix = get(a.ax,'position');
    set(a.ax,'units', oldunits);
    
    if ~a.whitebg,
      set(a.ax, 'color', [0 0 0]);
    else
      set(a.ax, 'color', [1 1 1]);
    end
    
    % subsample to make plotting quicker, but still allow for multiple
    % points/x-pixel value for steep lines
    axwidpix = axpospix(3);
    tpfdatl = size(a.parm.tpfdat, 1);
    if ~isempty(a.opt.subsampnpix) && (tpfdatl > axwidpix*a.opt.subsampnpix); % more samples than x * pixels
      subsampf = tpfdatl/(axwidpix*a.opt.subsampnpix);
      
      % low-pass filter (i.e. take moving average) of parms before subsampling
      smoothf = rectwin(round(subsampf./2));
      if smoothf > 1,
        a.parm.tpfdat = [a.parm.tpfdat(:,1) ... % don't bother smoothing time
                         filtfilt(smoothf, sum(smoothf), a.parm.tpfdat(:,(2:end)))];
      end
      subsampi = round(1:subsampf:tpfdatl);
      a.parm.tpfdat = a.parm.tpfdat(subsampi,:);
    end

    % if second tunecurve, bin data according to that histogram, assign
    % color based on bin.

    if ~isempty(a.parm.tpfdat)
      if numel(a.parm.tunecurveopt.tuneparams) > 1

        % get bin number of each value
        [dummy marker_bin] = histc(a.parm.tpfdat(:,3), ... 
                                   a.parm.tunecurveopt.histedges{2});

        % assign values from last ('edges') bin to bin 0
        marker_bin(marker_bin == length(dummy)) = 0;
        
        switch length(a.parm.tunecurveopt.histedges{2})
        
         case 4, % 3 bins: make middle bin bgcolor
          bin1 = marker_bin == 3;         
          bin2 = marker_bin == 1;
          bin0 = ~bin1 & ~bin2; % includes middle bin (#2)
          
         case 2 % 1 bin: make everything not in that bin bgcolor
          bin1 = marker_bin == 1;
          bin2 = [];
          bin0 = ~bin1;
          
         otherwise
          error('only 1 or 3 bins supported for color mapping of parm markers');
          
        end
        
        % plot bgcolor first, so it's in the background (we're using
        % 'painters' renderer, where these things matter)
        if ~isempty(bin0)
          mcol = a.opt.bgcolor;
          h = plot(a.ax, a.parm.tpfdat(bin0,1), a.parm.tpfdat(bin0,2),'.',...
                   'linestyle', 'none', 'marker', a.opt.marker, ...
                   'markeredgecolor', mcol,...
                   'markerfacecolor', mcol,...
                   'markersize', a.opt.markersize);
          hs = [hs; h];
        end
        
        if ~isempty(bin2)
          mcol = a.opt.ptcolor(2,:);
          h = plot(a.ax, a.parm.tpfdat(bin2,1), a.parm.tpfdat(bin2,2),'.',...
                   'linestyle', 'none', 'marker', a.opt.marker, ...
                   'markeredgecolor', mcol,...
                   'markerfacecolor', mcol,...
                   'markersize', a.opt.markersize);
          hs = [hs; h];
        end
        
        if ~isempty(bin1)
          mcol = a.opt.ptcolor(1,:);
          h = plot(a.ax, a.parm.tpfdat(bin1,1), a.parm.tpfdat(bin1,2),'.',...
                   'linestyle', 'none', 'marker', a.opt.marker, ...
                   'markeredgecolor', mcol,...
                   'markerfacecolor', mcol,...
                   'markersize', a.opt.markersize);
          hs = [hs; h];
        end
        
      else
        % no different marker colors
      
        h = plot(a.ax, a.parm.tpfdat(:,1), a.parm.tpfdat(:,2),'.',...
                 'color', a.opt.ptcolor(1,:));
        hs = [hs; h];
        
      end
    end
      
% $$$     % scatter has problems! If you pass in RGB values when using the
% $$$     % painters renderer, it barfs
% $$$     h = scatter (a.parm.tpfdat(:,1), ...
% $$$                  a.parm.tpfdat(:,2),20,mcol,'+');

  end

  
  % draw legend even if we don't draw the parm
  
  if a.opt.legend,
    
    if ~isempty(a.opt.dispname),
      name = a.opt.dispname;
    else
      if ~isempty(a.parm.tunecurveopt.name)
        name = a.parm.tunecurveopt.name;
      else
        name = a.parm.tunecurveopt.tuneparams{1};
      end
    end
             
    h = contdrawlegend('chanlabels', name,... 
                       'dispname', a.opt.dispname,...
                       'plottype', 'image',...
                       'cmaptop', a.opt.ptcolor(1,:),...
                       'whitebg', a.whitebg,...
                       'ax', a.ax);
    hs = [hs; h];
  end

  
  %%% set xlim (usu according to requested timewin_plot)
  if ~isempty(a.xlim),
    xlim(a.ax, a.xlim),
  end
  
  if ~isempty(a.opt.datalim),
    ylim(a.ax, a.opt.datalim),
  end

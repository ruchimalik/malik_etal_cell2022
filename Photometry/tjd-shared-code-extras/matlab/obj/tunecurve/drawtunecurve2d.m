function [figh axh cbaxh arh] = drawtunecurve2d(varargin)
% DRAWTUNECURVE2D plot a 2-D tuning curve from a tunecurves object
%
% todo:


  [figh axh cbaxh arh] = deal([]);
  
  a = struct('tunecurve', [],...
             'clnos', [],...  
             'ax', [],...
             'fig', [],...
             'e', [],...
             'opt', []);
  
  a = parseArgs(varargin, a);

  %%% parse arguments / error checking
  % shorter names
  tc = a.tunecurve;
  e = a.e;
 
  if isempty(a.opt),
    a.opt = mktunecurve2ddrawopt();
  end
  
  if isempty(a.clnos)
    error('must specify clnos');
  end

  if ~isempty(a.ax) && ~isempty(a.fig)
    error('Can''t provide both ax and fig');
  end
  
  if ~isempty(a.ax) || ~isempty(a.fig)
    if length(a.clnos) > 1 || length(a.fig) > 1,
      error('if ax or fig is provided, can only plot one cluster');
    end
  end

  if ~isempty(a.e) && ~strcmp(e.desc, tc.edesc),
    warning('''desc'' field of ''e'' doesn''t match ''edesc'' field of ''tunecurve''');
  end

  if a.opt.whitebg,
    textcol = 'k';
  else
    textcol = 'w';
  end
  
  % how many dimensions for tuning curve
  tcsize = size(tc.tcdats);
  tcndims = length(tcsize);

  if tcndims ~= 3 ... % clno * pos_x * pox_y, e.g.
        && ~(tcndims == 4 && tcsize(4) <= 3), % clno * pos_x * pox_y, * 2

    error(['Only 3-D tuning curves (or 4-D with up to 3 bins in 4th dim) ' ...
           'can be plotted (e.g. [clno x pos_x x pos_y x vel]']);
  end

  
  %%% draw the plots already!
  % get indexes of a.clnos to plot
  if islogical(a.clnos),
    a.clnos = find(a.clnos);
  end
  [dummy cli] = ismember(a.clnos,tc.clnos); %#ok

  % range of data in tcdats
  imxlim = tc.tunecurveopt.histedges{2}([1 end]);
  imylim = tc.tunecurveopt.histedges{1}([1 end]);
  
  % tunecurve param limits (shared across figs)
  if ~isempty(a.opt.autolim) && ~isempty(e),
    xlim = [min(e.track.splinevals(:,1)) - a.opt.autolim,...
            max(e.track.splinevals(:,1)) + a.opt.autolim];

    ylim = [min(e.track.splinevals(:,2)) - a.opt.autolim,...
            max(e.track.splinevals(:,2)) + a.opt.autolim];
  else
    if ~isempty(a.opt.xlim),
      xlim = a.opt.xlim;
    else
      xlim = imxlim;
    end
    if ~isempty(a.opt.ylim),
      ylim = a.opt.ylim;
    else
      ylim = imylim;
    end
  end
  
  % shared datalim/zlims across cls
  if ~isempty(a.opt.datalim),
    datalim = a.opt.datalim;
    zlimscale = false;
  else
    if a.opt.shareddatalim, 
      if tcndims == 4, % we want max of sum across channels
                       % (cf. draw2d->zlimscale)
        sumall = sum(tc.tcdats(cli,:,:,:),4);
      else % ndim = 3, just take max of data
        sumall = tc.tcdats(cli,:,:);
      end
      datalim = [0 max(sumall(:))];
      zlimscale = false;
    else
      datalim = [];
      zlimscale = true;
    end
  end
    
  for j = 1:length(cli);
    
    clij = cli(j);
    
    if ~isempty(a.fig),
      figh(j) = a.fig; %#ok
      clf(figh(j), 'reset');
      axh(j) = axes('parent', figh(j)); %#ok
    else 
      if ~isempty(a.ax),
        figh(j) = ancestor(a.ax,{'figure' 'uipanel'});
        axh(j) = a.ax;
      else
        figh(j) = figure;
        axh = axes('parent', figh(j));
      end
    end
    
    % get data to plot
    if ~a.opt.drawspikes
      tcdats = squeeze(tc.tcdats(clij,:,:,:)); % extra ':' are just ignored

      [hs zlims] = draw2d('ax', axh(j),...
                          'cmap',a.opt.colormap,...
                          'data', tcdats,...
                          'xhistendedges', imxlim,...
                          'yhistendedges', imylim,...
                          'datalims', datalim,...
                          'whitebg', a.opt.whitebg,...
                          'zlimscale', zlimscale); %#ok

      if a.opt.drawcolorscale,
        if tcndims == 3,
          colorbar('peer', axh(j), 'WestOutside')
        else % tcndims = 4
          
          set(axh(j), 'units', 'normalized');
          axpos = get(axh(j), 'position');
          cbaxh(j) = axes('units', 'normalized',...
                          'outerposition', [axpos(1:2) 0.25 0.25],...
                          'parent', ancestor(axh(j), {'figure', 'uipanel'})); %#ok
          
          colorbox('ax', cbaxh(j),...
                   'whitebg', a.opt.whitebg,...
                   'range1', zlims(1,:),...
                   'range2', zlims(2,:),...
                   'label1', 'A\rightarrowB (Hz)',...
                   'label2', 'B\rightarrowA (Hz)');
        end
      end
    end
    
    if ~isempty(e) && a.opt.drawtrackmap,
      trackmap('ax', axh(j),...
               'e', e,...
               'plots', {'track' 'endlabels'},...
               'linewidth', 1,...
               'fontsize', 10,...
               'whitebg', a.opt.whitebg,...
               'color', a.opt.trackmapcol);
    end

    if a.opt.drawspikes
      if isempty(e)
        error('Must provide ''e'' to draw spikes');
      end
      
      % draw mix of AB/BA spikes in color
      spikei = getepd(e.cl(a.clnos), 'posi');
      switch tc.tunecurveopt.tuneparams{3}
       case 'lvelf'
        spikes = getepdi(e.pos,spikei,'centx', 'centy', 'lvelf');
        th = tc.tunecurveopt.histedges{3}([2 end-1]);
        vel = spikes(:,3);
        spikes(vel<th(1),3) = 0;
        spikes(vel>th(1) & vel<th(2),3) = 1;
        spikes(vel>th(2),3) = 2;
       case 'lspeedf'
        spikes = getepdi(e.pos,spikei,'centx', 'centy', 'lspeedf');
        th = tc.tunecurveopt.histedges{3}(1);
        vel = spikes(:,3);
        spikes(vel<th(1),3) = 1;
        spikes(vel>th(1),3) = 3;
       otherwise
        error(['only ''lvelf'' and ''lspeedf'' currently supported for ' ...
               'drawspikes']);
      end
      sh = scatter(axh(j), spikes(:,1), spikes(:,2), 5, spikes(:,3), 'filled');
      % make a red gray blue magenta colormap
      colormap(axh(j),[1 0 0;... 
                        0.6 0.6 0.6; ...
                        0 0 1; ...
                        0.8 0 1]);
      set(axh(j), 'clim', [0 3]);
    end
    
    
    if a.opt.drawarrows && ~isempty(e) && tcndims == 4,
      % draw arrows to indicate direction of colors
      % (use arrow.m from Mathworks file exchange)
      arendoff_cm = 10; % location of start of arrow relative to end
      arlen_cm = 30; % length of arrow in cm
      artracki = findnearest(e.track.l, ...
                             [arendoff_cm ...
                          (arendoff_cm+arlen_cm)  ...
                          (e.track.l(end)-arendoff_cm) ...
                          (e.track.l(end)-arendoff_cm - arlen_cm)]');
                          
      aroff_cm = -10; % distance to offset arrow (normal to track)
      aroffx_cm = aroff_cm * -sin(deg2rad(e.track.dir(artracki([1 3]))) + [0; pi]);
      aroffy_cm = aroff_cm * cos(deg2rad(e.track.dir(artracki([1 3]))) + [0; pi]);

      arx = e.track.splinevals(artracki,1) +...
            aroffx_cm([1 1 2 2]);
      ary = e.track.splinevals(artracki,2) +...
            aroffy_cm([1 1 2 2]);
      
      % stupid arrow.m will only draw into gca
      oldgca = gca;
      set(gcf, 'currentaxes', axh(j));
      
      arh(1) = arrow([arx(1) ary(1)], [arx(2) ary(2)], ...
                     'edgecolor', 'b',...
                     'facecolor', 'b',...
                     'width', 2);
      arh(2) = arrow([arx(3) ary(3)], [arx(4) ary(4)], ...
                     'edgecolor', 'r',...
                     'facecolor', 'r',...
                     'width', 2);
      set(gcf, 'currentaxes', oldgca);
    end
    
    if a.opt.drawtitle,
      title_str = [tc.edesc ' cl-' num2str(clij)];
      if ~isempty(e)
        title_str = [title_str ' (' e.cl(clij).name ').'];
      end
      title_str = [title_str 'Max FR: ' num2str(max(tcdats(:))) ' Hz. '];
      
      title(axh(j), title_str,...
            'interpreter', 'none', ...
            'color', textcol);
    end

    axis(axh(j), 'equal');
    axis(axh(j), 'tight');
    axis(axh(j), 'ij');
    set(axh(j), 'xlim', xlim);
    set(axh(j), 'ylim', ylim);

  end
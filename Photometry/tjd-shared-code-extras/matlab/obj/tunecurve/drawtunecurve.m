function hs = drawtunecurve(varargin)
% DRAWTUNECURVE plot tuning curves from a tunecurve object
% many subplots like linpfall
%
% todo:
% stacked bars like linpfmult? 
% 2D plot like linpf?
% polar plot for headdir data? (just to linear for now)

  hs = [];
  
  a = struct('tunecurve', [],...
             'clnos', [],...  % [] = use all
             'ax', [],...
             'fig', [],...
             'e', [],...
             'opt', []);
  
  a = parseArgs(varargin, a);
  
  % shorter names
  tc = a.tunecurve;
  e = a.e;
 
  if isempty(a.opt),
    a.opt = mktunecurvedrawopt();
  end
  
  if isempty(a.clnos)
    clnos = tc.clnos;
  else
    clnos = a.clnos;
  end

  if islogical(clnos),
    clnos = find(clnos);
  end
  
  if ~isempty(a.ax)
    if length(a.clnos) > 1,
      error('if ax is provided, can only plot one cluster');
    end
    oneax = true;
  else
    oneax = false;
  end

  % how many dimensions for tuning curve
  tcndims = ndims(tc.tcdats);

  if tcndims >3, % clno * pos * vel, e.g.
    error('Only 1D and 2D tuning curves can be plotted');
  elseif size(tc.tcdats,3)>2 % clno * pos * 2vel bins, e.g.
    error(['Only 2D tuning curves with 2 bins in last dimension ' ...
           'supported']);
  end

  % get indexes of clnos to plot
  [dummy cli] = ismember(clnos,tc.clnos); %#ok

  % get data to plot
  tcdats = tc.tcdats(cli,:,:);
  
  % how many curves to plot
  ncl = length(clnos);
  
  % get biggest value of all the data we're going to plot
  ymax = max(reshape(tcdats,[],1));
  
  if ~isempty(a.opt.ylims);
    yl = a.opt.ylims;
  elseif a.opt.ylimsequal,
    if tcndims == 2
      yl = [0 1.05*ymax];
    else
      yl = 1.05 * [-ymax ymax];
    end
  else
    yl = [];
  end
    
  if oneax,
    subaxh(1) = a.ax;
    ax = a.ax;
    h = subf_tcplot(ax, e, tcdats(1,:,:), tc.tunecurveopt.histedges{1}, ...
                    clnos(1), a.opt, yl);

    if ~a.opt.botlabel
      set(ax, 'xticklabel', {})
      xlabel(ax, []);
    else
      
      set(ax, 'xtick', 0:200:1200);
      % convert cm to m
      if strcmp(a.opt.units,'m');
        set(ax, 'xticklabel', {'0' '2' '4' '6' '8' '10' '12'}); 
      end
      xlabel(ax, ['Distance (' a.opt.units ')']); 
      hs = [hs;h];
    end
  else
    for fig = 1:ceil(ncl./a.opt.plotsperfig);
      if fig == 1 && ~isempty(a.fig),
        % use user figure as first figure, if provided
        set(0,'currentfigure', a.fig);
        clf;
      else
        figure;
      end
      %% set up figure
      set(gcf,'color',[1 1 1]);
      %print our backgrounds  
      set(gcf, 'InvertHardCopy', 'off');

      subp = 0;
      for k = ((fig-1)*a.opt.plotsperfig)+1 : min([ncl fig*a.opt.plotsperfig]);
        subp = subp + 1;

        subaxh(k) = subplot(a.opt.plotsperfig,1,subp, 'align');
        subax = subaxh(k);
        hs = [hs;subax];      

        h = subf_tcplot(subax, e, tcdats(k,:,:), tc.tunecurveopt.histedges{1}, ...
                        clnos(k), a.opt, yl);
        hs = [hs;h];      
        
        if ~a.opt.botlabel || (subp ~= a.opt.plotsperfig && subp ~= min([ncl fig*a.opt.plotsperfig]));
          %remove x labels for all but bottom plot
          set(subax, 'xticklabel', {})
        else
          set(subax, 'xtick', 0:200:1200);
          if strcmp(a.opt.units,'m');
            % convert cm to m
            set(subax, 'xticklabel', {'0' '2' '4' '6' '8' '10' '12'}); 
          end  
          xlabel(subax, ['Distance (' a.opt.units ')']); 
        end
      
      end
    end
  end
  
      
  if a.opt.drawarrows && tcndims == 3,
    % draw arrows to indicate direction of colors (only for 2D tcs)
    % (use arrow.m from Mathworks file exchange)

    % draw only into first axes
    axh = subaxh(1);
    
    arendoff_cm = 25; % location of start of arrow relative to end
    arlen_cm = 100; % length of arrow in cm
    arposy_frac = 0.6; % height of arrow on y-axis
    
    xl = xlim(axh);
    arposy_data = arposy_frac * ylim(axh);
    
    if ~isempty(e),
      endx = e.track.length;
    else
      endx = xl(2);
    end
    
    arx = [cumsum([xl(1)+arendoff_cm arlen_cm])...
           cumsum([endx-arendoff_cm -arlen_cm])];
    
    ary = arposy_data([2 2 1 1]);
    
    % arrow.m will only draw into gca
    oldgca = gca;
    oldgcf = gcf;
    set(0, 'currentfigure', ancestor(subaxh(1), 'figure'));
    set(gcf, 'currentaxes', subaxh(1));
    
    arh(1) = arrow([arx(1) ary(1)], [arx(2) ary(2)], ...
                   'edgecolor', 'b',...
                   'facecolor', 'b',...
                   'width', 2);
    arh(2) = arrow([arx(3) ary(3)], [arx(4) ary(4)], ...
                   'edgecolor', 'r',...
                   'facecolor', 'r',...
                   'width', 2);
    set(0, 'currentfigure', oldgcf);
    set(gcf, 'currentaxes', oldgca);
  end
  
  
  
function hs = subf_tcplot(ax, e, tcdat, histedges1, clno, dopt, yl)
  
  hs = [];
  
  if ~isempty(e) && clno > 0,
    clname = ['(' e.cl(clno).name ')'];
  else
    clname = '';
  end

  % if only one 'slice', plot it; if two slices, plot first below line,
  % second above
  if size(tcdat,3) > 1,
    tcdat(1,:,1) = -tcdat(1,:,1);
  end

  if max(tcdat(:)) % test for non-zero values
    for k = 1:size(tcdat,3); % either 1 or 2
      h = area(ax, ctrs(histedges1), tcdat(1,:,k));
      set(h, 'facecolor', dopt.col(k,:));
      set(h, 'edgecolor', 'none');
      hs = [hs;h];
      hold (ax, 'on');
    end      
    % plot line at center
    plot(ax, [histedges1(1) histedges1(end)], [0 0], 'k');
  end
  
  set(ax, 'color', dopt.bgcol);
  
  % remove ticks, add xgrid
  set(ax, 'TickLength', [0 0]);
  set(ax, 'XGrid', 'on');
  set(ax, 'YGrid', 'on');
  set(ax, 'YMinorGrid', 'on');
  
  % label the y-axis
  if dopt.leftlabel
    set(get(ax, 'ylabel'), 'string', sprintf('cl-%d,\n %s', clno, clname));
  else
    set(ax, 'ytick', []);
  end
  
  % use track length as xlim, if available
  if ~isempty(e)
    xlim(ax, [0 e.track.length]);
  else
    lims=objbounds(ax);
    if ~isempty(lims),
      xlim(ax, lims(1:2));
    end
  end
  
  if dopt.ylogscale,
    set(ax, 'yscale', 'log');
  end
  
  % set ylims
  if clno == 0,
    ylim(ax, [-0.1 0.1]);
  elseif ~isempty(yl) && ~all(yl == 0)
    ylim(ax, yl);
  end

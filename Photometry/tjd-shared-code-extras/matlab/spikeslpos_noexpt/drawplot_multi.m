function [cache hs out] = drawplot_multi(varargin)
% DRAWPLOT draw a plot object created by slmkplots (spikeslpos)

  [out hs] = deal([]);
  
  a = struct(...
      'plot', [],...
      'cache', [],...
      'ax', [],...
      'mintimebin',[],...
      'dpi', []);
  
  a = parseArgsLite(varargin, a);
  
  if isempty(a.plot),
    return;
  end
  
  cache = a.cache;
  
  % initialize struct to store objects in
  objs = [];

  if isempty(a.ax),
    a.ax = gca;
  end

  hold(a.ax, 'on');
  
  % fix for Matlab R2015a--allows seg rectangles to show through
  set(a.ax, 'color', 'none');
  
  switch a.plot.type
   case 'parmest'
    
    e = subf_gete(a.plot.dat);
    
    [objs.parmest objs.tunecurve objs.clhist] = ...
        mkparmest_multi('e',e,...
                        'epochi', a.plot.dat.epochi,...
                        'clnos',a.plot.dat.clnos,...
                        'clperm',a.plot.dat.clperm,...
                        'timewin', a.plot.dat.timewin,...
                        'timebinsize', max([a.mintimebin a.plot.dat.timebinsize]),...
                        'overlapsperbin', a.plot.dat.overlapsperbin,...
                        'parmestopt', a.plot.parmestopt,...
                        'tunecurveopt', a.plot.tunecurveopt,...
                        'cache',cache);
    
    % get behavior parameters object corresponding to the tunecurve used
    try
      objs.parm = ...
          mkparm_multi('e',e,...
                       'epochi', a.plot.dat.epochi,...
                       'timewin', a.plot.dat.timewin,...
                       'tunecurveopt', a.plot.tunecurveopt,...
                       'cache', cache);
    catch
      warning(['can''t get parm corresponding to tunecurve used to build ' ...
               'ordered rasters']);
      objs.parm = [];
    end
    
    hs = drawparmest('parmest', objs.parmest, ...
                     'opt', a.plot.parmestdrawopt, ...
                     'parm', objs.parm,...
                     'parmdrawopt', a.plot.parmdrawopt,...
                     'whitebg', a.plot.whitebg,...
                     'ax', a.ax,...
                     'xlim', a.plot.dat.timewin_plot);
    
   case 'raster'
    
    e = subf_gete(a.plot.dat);
    
    set(a.ax, 'units', 'pixels');
    axpos = get(a.ax,'position'); % for tick widths,etc
    
    if ~isempty(a.dpi)
      set(a.ax, 'units', 'inches');
      axpos_in = get(a.ax, 'position');
      imwidpix = axpos_in(3) * a.dpi;
      imhtpix = 4 * axpos_in(4) * a.dpi; %over-sample, as lost rows are
                                         %no biggie
      set(a.ax, 'units', 'pixels'); % reset
    else
      imwidpix = axpos(3);
      imhtpix = axpos(4);
    end
    
    [objs.raster objs.tcpeak objs.tunecurve objs.clhist out.cls_ranked out.clcolors] =...
        mkraster_multi('e',e,...
                       'epochi', a.plot.dat.epochi,...
                       'clnos', a.plot.dat.clnos,...
                 'clperm', a.plot.dat.clperm,...
                 'timewin', a.plot.dat.timewin,...
                 'parmwin', [],... % use full range of tunecurves
                 'imwidpix', imwidpix,... % make rast image to fill axis
                 'imhtpix', imhtpix,... % make rast image to fill axis
                 'tunecurveopt', a.plot.tunecurveopt,...
                 'parmdrawopt', a.plot.parmdrawopt,...
                 'tcpeakopt', a.plot.tcpeakopt,...
                 'rasteropt', a.plot.rasteropt,...
                 'whitebg', a.plot.whitebg,...
                 'cache', cache);
    
    
    % get behavior parameters object corresponding to the tunecurve used
    try
      objs.parm = ...
          mkparm_multi('e',e,...
                       'epochi', a.plot.dat.epochi,...
                       'timewin', a.plot.dat.timewin,...
                       'tunecurveopt', a.plot.tunecurveopt,...
                       'cache', cache);
    catch
      warning(['can''t get parm corresponding to tunecurve used to build ' ...
               'ordered rasters']);
      objs.parm = [];
    end
    
    hs = drawraster('raster', objs.raster,...
                    'opt', a.plot.rasterdrawopt,...
                    'parm', objs.parm,...
                    'parmdrawopt', a.plot.parmdrawopt,...
                    'whitebg', a.plot.whitebg,...
                    'ax', a.ax,...
                    'xlim', a.plot.dat.timewin_plot);
    
   case 'parm'
    
    e = subf_gete(a.plot.dat);

    objs.parm = ...
        mkparm_multi('e',e,...
                     'epochi', a.plot.dat.epochi,...
                     'timewin', a.plot.dat.timewin,...
                     'tunecurveopt', a.plot.tunecurveopt,...
                     'cache', cache);

    if isempty(a.plot.parmdrawopt.datalim),
      a.plot.parmdrawopt.datalim = [0 e.track.length];
    end
    
    h =  drawparm('parm', objs.parm,...
                  'opt', a.plot.parmdrawopt,...
                  'whitebg', a.plot.whitebg,...
                  'ax', a.ax,...
                  'xlim', a.plot.dat.timewin_plot);

    hs = [hs; h];
    
   case 'cont'
    
    objs.cont = ...
        mkcont('contdata', a.plot.dat.contdata,...
               'contvar', a.plot.dat.contvar,...
               'chans', a.plot.dat.chans,...
               'chanlabels', a.plot.dat.chanlabels,...
               'timewin', a.plot.dat.timewin,...
               'contopt', a.plot.contopt,... % currently empty
               'cache', cache);
   
    hs = drawcont('cont', objs.cont,...
                  'opt', a.plot.contdrawopt,...
                  'dpi', a.dpi,...
                  'plottypeno', a.plot.plottypeno,...
                  'whitebg', a.plot.whitebg,...
                  'ax', a.ax,...
                  'xlim', a.plot.dat.timewin_plot);
    
   case 'specgram'

        
    [objs.specgram objs.cont] = ...
        mkspecgram('contdata', a.plot.dat.contdata,...
                   'contvar', a.plot.dat.contvar,...
                   'chans', a.plot.dat.chans,...
                   'chanlabels', a.plot.dat.chanlabels,...
                   'timewin', a.plot.dat.timewin,...
                   'contopt', a.plot.contopt,... 
                   'specgramopt', a.plot.specgramopt,...
                   'cache', cache);
    
    hs = drawspecgram('sg', objs.specgram,...
                      'opt', a.plot.specgramdrawopt,...
                      'whitebg', a.plot.whitebg,...
                      'ax', a.ax,...
                      'xlim', a.plot.dat.timewin_plot);
    
   case 'none'
    
   otherwise
    error('unrecognized plot type');
  end
  
  % store any new data in the cache
  for otype = fieldnames(objs)';
    otype = otype{:};
    obj = objs.(otype);
    if ~isempty(obj) && ~obj.cache_hit,
      cache = mkcache('cache', cache, 'add_obj', obj);
    end
  end

  % invert y-axis if requested
  if a.plot.yflip,
    set(a.ax,'ydir','reverse');
  else
    set(a.ax,'ydir','normal');
  end

%%%%%%%%%%%%%%%
%% SUBFUNCTIONS

function e = subf_gete(dat)
  
  if ~isempty(dat.e),
    e = dat.e;
  elseif ischar(dat.ename) && ...
        evalin('base', ['exist(''' dat.ename ''', ''var'')']),
    e = evalin('base',dat.ename);
  end

  if ~exist('e','var') || ~isstruct(e),
    error('no ''plot.dat.e'' struct or ''plot.dat.ename'' provided');
  end
  
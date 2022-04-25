function [r pks tc clh cls_ranked clcolors] = mkraster(varargin)
% MKRASTER make a spike raster plot image
%
% $Id$
  
  % initialize outputs;
  [r pks tc clh cls_ranked clcolors] = deal([]); %#ok
  
  r = struct(...
      'type', 'raster',...
      ...
      'e',[],... % inputs
      'edesc', [],...
      'clnos', [],...
      'clperm', [],...
      'timewin',[],... % time range to compute
      'parmwin', [],... % parameter range to compute
      'imwidpix',[],... % width of requested image in pixels
      'imhtpix', [],... % height of requested image in pixels
      'tunecurveopt', [],...
      'parmdrawopt', [],...
      'tcpeakopt', [],...
      'rasteropt', [],...
      'whitebg', [],...
      'dim', 1,...
      ...
      'rastimg',[],... % outputs
      'timehistedgecenters', [],...
      'parmhistendedges', [],...
      'label',[],...
      ...
      'template', [],...
      'cache',[],...
      'cache_hit', false);
  
  r = parseArgsLite(varargin,r);
  
  % handle 'template' object (add it to cache, parse other args)
  r = obj_reparse(r,varargin);
  
  if islogical(r.clnos),
    r.clnos = find(r.clnos);
  end
  
  % clnos = [] means use all clusters (do before cachesearch);
  if isempty(r.clnos),
    r.clnos = 1:length(r.e.cl);
  end
  
  % find object from cache, if possible
  r = obj_cachesearch(r);

  if r.rasteropt.tickwidpix < 1,
    warning(['raster tick width less than 1 can lead to lost spikes on ' ...
             'screen and in printouts! Use with caution.']);
  end
    
  if r.cache_hit,
    r = obj_cleanup(r);
    return
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % we didn't get a cache hit, recompute.

  %%% setup 

  %% colormap stuff
  
  % 'double' raster images get pretty big in the cache, so use more
  % compact datatype.
  %
  % use int16 (1-indexed) rather than uint16 for now because uint16 colormap
  % indices are zero-indexed , which is a pain!

  imtype = 'int16';

  % we assume that the first 4 colors in the cmap are:
  % 1: plotbg_black : the color of the plot background when ~whitebg
  % 2: rastbg_black : the color of the raster band bg when ~whitebg
  %
  % 3: plotbg_white : the color of the plot background when whitebg
  % 4: rastbg_white : the color of the raster band bg when whitebg

  if ~r.whitebg;
    plotbgcol = 1;
    rastbgcol = 2;
  else
    plotbgcol = 3;
    rastbgcol = 4;
  end
  
  nclcolors = r.rasteropt.nclcolors;
  nbgcolors = 4;

  % because of later reshaping of colormap in 'electrode' rastcolormode, below
  nclcolors_factor = 4;
  if rem(nclcolors,nclcolors_factor)
    error(['Sorry, but the number of cluster colors must be evenly divisible ' ...
           'by ' num2str(nclcolors_factor) '. It just does.']);
  end
  
  %%% get tuning curves
  % do this before setting parmwin
  tc = mktunecurve('e', r.e,...
                   'edesc', r.edesc,...
                   'clnos', r.clnos,...
                   'tunecurveopt', r.tunecurveopt,...
                   'cache', r.cache);
  
  % save for use in mktcpeaks (hackish, ick!)
  r.cache = mkcache('cache', r.cache, 'add_obj', tc);
  
  % if no parmwin requested, use full range of tuning curves
  if isempty(r.parmwin),
    parmwin = tc.tunecurveopt.histedges{r.dim}([1 end]);
  else
    parmwin = r.parmwin;
  end
  
  % pixel/data scaling factors
  dataperpix_ht = abs(diff(parmwin))/r.imhtpix;
  dataperpix_wid = abs(diff(r.timewin))/r.imwidpix;

  % pixel values may be fractions, draw image inside axes
  % (x-axis pixel alignment taken care of by choosing rtimebin to be a multiple
  % of pixel width, below)
  rastimght = floor(r.imhtpix);
  % +/-0.5 pixel offset for pixel centers for image plotting. 
  % Since rastimght may be rounded from requested parmwin, use rastimght
  % rather than r.imhtpix for top edge.
  r.parmhistendedges = parmwin;

  % amount (in pixels) to extend each tick above/below actual peak
  tickoffset = floor((r.rasteropt.tickhtpix-1)/2);

  % shortcut
  ncls = length(r.clnos);

  
  %%% get peaks of tuning curves
  pks = mktcpeak('e', r.e,...
                 'edesc', r.edesc,...
                 'clnos', r.clnos,...
                 'tunecurveopt', r.tunecurveopt,...
                 'tcpeakopt', r.tcpeakopt,...
                 'cache', r.cache);

  
  %%% cluster colors
  % 
  % Have to wait until we have pks to do this...
  % nclcolors = available cluster colors 
  
  clcolors = zeros(1,ncls,imtype); % zero will be invalid

  switch r.rasteropt.colormode
   
   case 'electrode'

    % assign colors by electrode
    for j = 1:length(r.e.electrode);
      clcolors(r.e.electrode(j).cls) = ...
          mod(j,nclcolors) + 1;
    end

    % shuffle the colors, so that nearby electrodes don't have nearby colors (if
    % there are only a few tetrodes, e.g.) Do this deterministically so it will
    % look the same with each redraw.
    shuffi = reshape(1:nclcolors, [], nclcolors_factor)';
    clcolors = shuffi(clcolors) + nbgcolors;

   case 'clno'
    % assign colors in cluster order
    for j = 1:ncls,
      clcolors(j) =...
          mod(j,nclcolors)...
          +1;
    end
    
    % shuffle as for 'electrode'
    shuffi = reshape(1:nclcolors, [], nclcolors_factor)';
    clcolors = shuffi(clcolors) + nbgcolors;
    
   case 'firstpeak'
    % assign colors cyclically to clusters in order of their first peak in
    % 'place' (helps to avoid case of 2 cells with overlapping peaks having
    % same color, at least in one-peak case)
    
    firstpeaks = zeros(size(clcolors)); 
    for j = 1:ncls;
      if isempty(pks.pkdats{j}),
        firstpeaks(j) = Inf;
      else
        firstpeaks(j) = pks.pkdats{j}(1,1);
      end
    end
    [discard peakrank] = sort(firstpeaks); %#ok
    
    for j = 1:length(peakrank),
      clcolors(peakrank(j)) =...
          mod(j,nclcolors)...
          +nbgcolors+1;
    end
    
   case 'gray'
    % black ticks on whitebg, white ticks on blackbg
    if r.whitebg,
      clcolors(:) = 1;
    else
      clcolors(:) = 3;
    end
      
   case {'parm' 'parmshades'}
    % assign all clusters the color of the parm being reconstructed
    % (i.e. for inbound, make all clusters red)
    
    % major HACK, but fun: the raster color is specified in indexed color, drawn
    % from a user-provided colormap, but ptcolor is an rgb value. So we take
    % advantage of the fact that the cmap is required to be linear in HSV to
    % find the closest indexed color
    
    % first get hue.
    ptcolor_hsv = rgb2hsv(r.parmdrawopt.ptcolor(1,:));
    ptcolor_hue = repmat(ptcolor_hsv(1),ncls,1);
    
    % jitter hue if requested, so we can still tell some cells apart
    if strcmp(r.rasteropt.colormode, 'parmshades'),
      ptcolor_hue = mod(ptcolor_hue+linspace(-0.1, 0.1,ncls)',1);
    end    

    % use the hue as a lookup into the cmap
    clcolors = round(ptcolor_hue*(nclcolors-1)+1)...
        +nbgcolors; % correct for reserved colors

   otherwise
    error('unknown ''rastcolormode'' type');
  end

  %%% find spikes in time window / generate raster histograms
  
  % set up histogram bins for  rasters
  rtimebin = r.rasteropt.tickwidpix * dataperpix_wid;
  
  clh = mkclhist(...
      'e', r.e,...
      'clnos', r.clnos,...
      'timewin', r.timewin,...
      'timebinsize', rtimebin,...
      'binarize', true,...
      'cache', r.cache);
  
  % reorder the spike records according to clperm
  if ~isempty(r.clperm)
    clh.clhists = clh.clhists(r.clperm,:);
  end
  
  % save for output
  r.timehistedgecenters = clh.histedgecenters;
  
  %%% initialize raster image to bgcolor 
  r.rastimg = repmat(cast(plotbgcol,imtype),rastimght,clh.nbins);
    
  if isempty(r.clnos),
    return,
  end
  
  %%% build / draw raster image

  % do ranked plots
  if r.rasteropt.ranked,
    
    % find all clusters with a first peak
    pk1il = []; 
    for j = 1:ncls,
      if ~isempty(pks.pkdats{j})
        % append a row with the cl index and the first peak location
        pk1il = [pk1il; j pks.pkdats{j}(1,1)];
      end
    end
    
    % get the cl indexes in the 
    pk1il_sorted = sortrows(pk1il,2);
    cls_ranked = pk1il_sorted(:,1);
    ncls_ranked = length(cls_ranked);

    if ncls_ranked > 0,
      r.parmhistendedges = [0.5 ncls_ranked+0.5];
     else
      r.parmhistendedges = [0 1];
    end
    
    if ~isempty(cls_ranked),
      ncls_ranked = length(cls_ranked);
      
      rank_spacing = rastimght ./ ncls_ranked;
      
      for k = 1:ncls_ranked
        cli = pk1il_sorted(k,1);
        r.rastimg(floor((k-1).* rank_spacing)+1:floor(k.*rank_spacing), ...
                  clh.clhists(cli,:)) = ...
            clcolors(cli);
      end
    end

  else
    %%% Draw gray bg for all raster bands (even those w/ no spikes).
    %
    % Had to separate out this step so that bg doesn't overwrite previous
    % rasters (i.e. so that all bgs are 'behind' all spike ticks).
    if r.rasteropt.rasterbg,
      tickcol = false(rastimght,1);

      for j = 1:ncls,
        for pktop = 1:r.rasteropt.maxpks,
          if ~isempty(pks.pkdats{j}) && length(pks.pkdats{j}(:,1)) >= pktop,
            
            % get peak to plot
            pk = pks.pkdats{j}(pktop,1);
            
            % convert from data to pixel coords
            pkpix = round((pk - parmwin(1)) / dataperpix_ht);
            
            % crop to rasterimage
            ticktop = pkpix + tickoffset;
            tickbot = pkpix - tickoffset;

            if tickbot > rastimght,
              continue, % don't draw this band
            else
              tickbot = max([1 tickbot]);
            end
            
            if ticktop < 1,
              continue, % don't draw this band
            else
              ticktop = min([rastimght ticktop]);
            end

            tickcol(tickbot:ticktop) = true;
            
          end
        end
      end
      % fill image array with raster bg color (2)
      r.rastimg(tickcol,:) = rastbgcol;
    end
    
    %%% now draw ticks for spikes

    for j = 1:ncls,
      if any(clh.clhists(j,:)), % any spikes to plot?
        tickcol = false(rastimght,1);
        for pktop = 1:r.rasteropt.maxpks,
          if ~isempty(pks.pkdats{j}) && length(pks.pkdats{j}(:,1)) >= pktop,
            
            % get peak to plot
            pk = pks.pkdats{j}(pktop,1);
            
            % convert from data to pixel coords
            pkpix = round((pk - parmwin(1)) / dataperpix_ht);
            
            % crop to rasterimage
            ticktop = pkpix + tickoffset;
            tickbot = pkpix - tickoffset;
            if tickbot > rastimght,
              continue, % don't draw this band
            else
              tickbot = max([1 tickbot]);
            end
            if ticktop < 1,
              continue, % don't draw this band
            else
              ticktop = min([rastimght ticktop]);
            end
            tickcol(tickbot:ticktop) = true;
          end
        end
        % fill image array
        r.rastimg(tickcol,clh.clhists(j,:)) = clcolors(j);
      end
    end
  end  
  %%% create label for plotting from tunecurve name
  r.label = tc.tunecurveopt.name;
  
  r = obj_cleanup(r);

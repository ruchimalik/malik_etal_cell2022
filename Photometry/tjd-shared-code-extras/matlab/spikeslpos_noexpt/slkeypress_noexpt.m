function slkeypress(fig)
% SLKEYPRESS - keypress callback function for spikeslpos 
%   
% $Id: slkeypress.m 481 2008-04-28 04:26:05Z tjd $
  
% set small factor to move on shift
sf = 0.2;

% default zoom factor
zf = 2;

% extra 'shift' zoom factor
zfmult = 5;

% # of screens to pan
movef = 0;

% zoom factor
zoomf = 1;

% get last key pressed
cc = get(fig,'CurrentCharacter');


% get panel containing last object clicked on
panels = findobj(gcbf,'type','uipanel');
if isempty(panels),
  return
end

for panel = panels(:)';
  sla = getappdata(panel,'slargs');
  sld = getappdata(panel,'sldata');

  switch lower(cc)

   case [] % as when shift key is pressed
    return
    
    
    %%% Pan

   case ']' % one screen right
    movef = +1;

   case '[' % one screen left
    movef = -1;

   case '}' % smallfactor screens right
    movef = sf;
    
   case '{' % smallfactor screens left
    movef = -sf;

    
    %%% Zoom  
   case '.' % zoom in by zoomfactor
    zoomf = zf;

   case '>' % shift-zoom in by big zoomfactor
    zoomf = zfmult*zf;
    
   case ',' % zoom out by zoomfactor
    zoomf = 1/zf;
    
   case '<' % shift-zoom out by big zoomfactor
    zoomf = 1/(zfmult*zf);

    
    %%% Viewlists  
   case 28 % left arrow
    if ~isempty(sla.viewlist) && sla.viewlisti > 1,
      sla.viewlisti = sla.viewlisti-1;
      % we use viewlist only if no timewin is provided
      sla.timewin = [];
    end

   case 29 % right arrow
    if ~isempty(sla.viewlist) && sla.viewlisti < size(sla.viewlist,1);
      sla.viewlisti = sla.viewlisti+1;
      % we use viewlist only if no timewin is provided
      sla.timewin = [];
    end
    
    %%% Raster Stuff
    
   case {'b', 'B'} % toggle raster Background
    for k = 1:length(sla.plots),
      if isfield(sla.plots{k}, 'rasteropt'),
        sla.plots{k}.rasteropt.rasterbg = ~sla.plots{k}.rasteropt.rasterbg;
      end
    end
% $$$     sla.rasterbg = ~sla.rasterbg;
        
    
    %%% Other
    
    % toggle showing/hiding all segments
   case 's'
    nsl = numel(sla.seglists);
    if isempty(sla.plotsegs)||any(sla.plotsegs)
      sla.plotsegs = false(nsl,1);
    else
      sla.plotsegs = true(nsl,1);
    end

    
   case 'h' % full crosshair toggle
    if strcmp(get(fig,'pointer'), 'fullcrosshair'),
      set(fig,'pointer','arrow');
    else
      set(fig,'pointer','fullcrosshair');
    end
    
   otherwise
       % unrecognized key
    return
    
  end

  setappdata(panel,'sldata',sld); % set here in case movie is canceled

  sla = slargsmovezoom(sla,sld,'movef',movef,'zoomf',zoomf);
  
  % redraw
  spikeslpos_noexpt('argstruct', sla);
end

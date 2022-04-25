function slresize_multi(fig)
% SLRESIZE - redraw spikeslpos panels on resize
% $Id: slresize.m 419 2007-08-29 23:16:05Z tjd $

% get panel containing last object clicked on
panels = findobj(fig,'type','uipanel');

for panel = panels';
  
  % get slargs
  sla = getappdata(panel,'slargs');

  if isempty(sla) || ...
        strcmp(get(panel,'units'),'pixels'); % weird race condition with
                                             % resize function interrupts
    continue
  end
  
  if isstruct(sla),
    spikeslpos_multi('argstruct', sla, 'panel', panel);
  end
end
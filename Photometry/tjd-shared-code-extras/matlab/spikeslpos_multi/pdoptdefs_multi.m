function pdopt = pdoptdefs_multi(varargin)
  
  a = struct('e', [],...
             'epochi', []);
  a = parseArgsLite(varargin, a);
  
  if ~isempty(a.e),
    
    tracklim = [0 a.e.track(a.e.epoch(a.epochi).tracki).length];
  else
    tracklim = [];
  end
  
% parmdrawopt for plotting behavior
pdopt.generic = mkparmdrawopt('datalim', tracklim);

pdopt.outbound = mkparmdrawopt('ptcolor', [0 0 0.8],...
                               'datalim', tracklim);
pdopt.inbound = mkparmdrawopt('ptcolor', [0.8 0 0],...
                              'datalim', tracklim);
pdopt.place = mkparmdrawopt('ptcolor', [0.8 0 0.8],...
                            'datalim', tracklim);
pdopt.pos_x_vel = mkparmdrawopt('ptcolor', [0 0 0.8; 0.8 0 0 ],...
                                'datalim', tracklim);

function opt = mkrasteropt(varargin)
% MKRASTEROPT make an options struct for mkraster
  opt = struct(...
      'ranked', false,...
      'tickwidpix', 1,...
      'tickhtpix', 7,...
      'maxpks', 1,...
      'rasterbg', true,...
      'colormode','parmshades',... % or 'parm', 'clno', 'electrode', 'firstpeak'
      'nclcolors',[]);
  
  opt = parseArgsLite(varargin, opt);

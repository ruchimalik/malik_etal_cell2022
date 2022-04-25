function opt = mktcpeakopt(varargin)
% MKPARMESTOPT - make an 'options' structure for parmest (parameter reconstruction)
%
% 'preset's available: 'bayes_zhang', 'basis', 'basis_pdf'
%
  
  opt = struct(...
      'template', [],...
      ...
      'peakfindgw', [] ,... % cm
      'peakfindth', 0); % Hz
  
  opt = parseArgsLite(varargin,opt);
  
  if ~isempty(opt.template),
    if ~strcmp(varargin{1},'template'),
      error('if pre-existing object is passed in as template, it must be first arg');
    end
    opt.template.template = [];
    % remove the 'template' from the args and use it
    opt = parseArgsLite(varargin(3:end),opt.template);
  end
  opt.template = [];
  
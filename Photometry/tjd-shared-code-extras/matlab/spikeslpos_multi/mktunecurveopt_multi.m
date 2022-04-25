function opt = mktunecurveopt_multi(varargin)
% MKTUNECURVEOPT - make an 'options' structure for tunecurve
  
  opt = struct(...
      'template', [],...
      ...
  'name', [],...
      'tuneparams', {{}},...
      'histedges', {{}},...
      'keepbins', {{}},...
      'runtimes',[-Inf +Inf],...
      'validate', 'useall',... % or 'useeven' or 'useodd'
      'smooth_sds', {{[]}}); % size of stdev of gaussian window in data
                           % units (e.g. cm, one per dimension)

  opt = parseArgsLite(varargin,opt);

  if ~isempty(opt.template),
    if ~strcmp(varargin{1},'template'),
      error('if pre-existing object is passed in, it must be first arg');
    end
    opt.template.template = [];
    % use pre-existing object as new template
    opt = parseArgsLite(varargin(3:end),opt.template);
  end

% $$$   if isempty(opt.tuneparams),
% $$$     error('no ''tuneparams'' specified');
% $$$   end

  if length(opt.tuneparams) ~= length(opt.histedges)
    error('There must be the same # of ''histedges'' as ''tuneparams''');
  end
  
  if ~isempty([opt.smooth_sds{:}]) && length(opt.tuneparams) ~= length(opt.smooth_sds)
    error('If provided, there must be the same # of ''smooth_sds''  as ''tuneparams''');
  end
  
  for k = 1:numel(opt.keepbins)
    opt.keepbins{k} = logical(opt.keepbins{k});
  end
  
  % when validating, we use even or odd seconds
  if ~any(strcmp(opt.validate, {'useall', 'useeven', 'useodd'}))
    error('Validate arg must be: useall, useeven, or useodd');
  end
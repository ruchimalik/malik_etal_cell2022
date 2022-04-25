function opt = mkparmestopt(varargin)
% MKPARMESTOPT - make an 'options' structure for parmest (parameter reconstruction)
%
% 'preset's available: 'bayes_zhang', 'basis', 'basis_pdf'
%
  
  a = struct(...
      'template', [],...
      ...
      'preset','bayes_zhang',...
      ...  
      'name',[],...
      'method', [],... % 'bayes' or 'basis'
      'history_gausswin_sd_bins', [],... % cm of gaussian to convolve with
          ...                          % prev estimate
      'occnorm',false,... % prior for position (occupancy correction)
      'columnnorm',false, ... % normalize each estimate
      'pdfnorm',false,... % normalize place-rate maps
      'expnorm',false,... % exp term in bayesian estimator
      'bsmallf',[],... % small factor to add to place-rate in bayesian est
      'binarizect',0, ... % count each cell only once per time bin
      ...
      'ratenorm',false, ... % normalize each cells' firing to their ave f.r.
      'decim_frac', [], ... % decimate a fraction of all spikes (cf mkclhist/mkparmest)
      'decim_type', 'dropspikes', ... % type of decimation ('dropspikes'/'scalerate')
      'basisnorm',false, ... % scale output by sum of tuning curves
      'infoscale', false,... % scale cells' firing by spatial info score
      ...                       % of tuning curves
      'tcfudgef', [],... % scale tuning curves
      'cellcountnorm',false, ... % scale each estimate by # of cells in timebin
      'spikecountnorm',false); % scale each estimate by # of spikes in timebin
          
  
  a = parseArgsLite(varargin,a);
  
  if ~isempty(a.template),
    if ~strcmp(varargin{1},'template'),
      error('if pre-existing object is passed in as template, it must be first arg');
    end
    opt.template.template = [];
    % remove the 'template' from the args and use it
    a = parseArgsLite(varargin(3:end),a.template);
  end
  
  % if 'preset' is given, set new defaults then re-parse arg list 
  switch a.preset
   case 'bayes_zhang'
    % simple memoryless Bayesian estimator, from Zhang et al, jneurophys 1998
    a.method = 'bayes';
    a.occnorm = false; % uniform prior
    a.expnorm = true; 
    a.columnnorm = true;
    a.columnnormfr = true; % normalize outb/inb maps together
   case 'basis'
    a.method = 'basis';
   case 'basis_pdf'
    a.method = 'basis';
    a.pdfnorm = true;
    a.pdfnormfr = true; % normalize outb/inb pdfs together
   case ''
    % no preset, just use args
   otherwise
    error('unrecognized preset');
  end
  
  % re-parse args, with new defaults set by preset
  % (allows us to override things set in presets)
  opt = parseArgsLite(varargin,a);
  
  % set name, if none given
  if isempty(opt.name) && ~isempty(opt.preset),
    if length(varargin) == 2,
      opt.name = opt.preset;
    else
      opt.name = [opt.preset '_mod'];
    end
  end

  opt.preset = '';
  opt.template = [];
  
  if isempty(opt.method),
    error('no ''method'' specified');
  end
  
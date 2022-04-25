function [tc occ] = dotunecurve(varargin)
% DOTUNECURVE - gets histograms of cells according to behav parameter
% return values in Hz. accepts multiple cls.
%
% [dotunecurve, histbins, cli] = 
%      'e', exp struct
%      'clnos', clusters to use
%      'opt', tunecurve options object cf mktunecurveopt.m
%
% $Id$

%INITIALIZE VARS
  
% constants

  % interval between position records, in seconds.
  posrecinterval = 1/30;

  % validate args 
  a = struct(...
      'e',[],...
      'clnos',[],...
      'opt',[]);

  a = parseArgsLite(varargin,a);
  
  if islogical(a.clnos)
    a.clnos = find(a.clnos);
  end
  
  %% use mktunecurveopt just for arg-checking
  a.opt = mktunecurveopt('template',a.opt);
  
  ndim = length(a.opt.histedges);
  ncls = length(a.clnos);
  
  % cell array of histedges to be passed into histn
  histedges = a.opt.histedges;

  % get eventual size of tc array for pre-allocation
  for dim = 1:ndim
    % there will be histedges-1 bins
    if ~isempty(a.opt.keepbins) && ~isempty(a.opt.keepbins{dim}),
      histblen{dim} = sum(a.opt.keepbins{dim});
    else
      histblen{dim} = size(histedges{dim},2)-1;
    end
  end
  
  % preallocate tuning curves (cl x n-d)
  tc = zeros(ncls,histblen{:}); 
  
  allpos = getepd(a.e.pos, 'time', a.opt.tuneparams{:});
  % select only pos records during requested time window
  allpos = allpos(allpos(:,1)> a.opt.runtimes(1) & allpos(:,1) < a.opt.runtimes(2),:);
  
  if ~strcmpi(a.opt.validate,'useall')
    even_pos_timesi = mod(allpos(:,1),2)<1;
    switch a.opt.validate
     case 'useeven',
      allpos = allpos(even_pos_timesi,:);
     case 'useodd',
      allpos = allpos(~even_pos_timesi,:);
     otherwise
      error('invalid value for tunecurveopt.validate');
    end
  end
  
  occ = histn(allpos(:,2:end),...
              a.opt.histedges, ...
              a.opt.keepbins, ...
              'noedgebin');
  
  if sum(occ) == 0,
    warning('No position records match requested time / param filtering')
    % tc already pre-allocated to be all zeros, so we can just return
    return;
  end
  
  for j = 1:ncls,
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % make tuning curves (uncorrected parameter histograms)
    
    % get relevant pos data (from e.pos) for each spike in e.cl
    cl_data = getepdi(a.e.pos, ...
                      getepd(a.e.cl(a.clnos(j)),'posi'), ...
                      'time', a.opt.tuneparams{:});

    % select only spikes during requested time window
    cl_data = cl_data(cl_data(:,1)> a.opt.runtimes(1) & cl_data(:,1) < a.opt.runtimes(2),:);
    % if validating, use only cls where floor(spiketime) is even or odd
    % (i.e. 1 second bins)
    if ~strcmpi(a.opt.validate,'useall')
      even_cl_timesi = mod(cl_data(:,1),2)<1;
      switch a.opt.validate
       case 'useeven',
        cl_data = cl_data(even_cl_timesi,:);
       case 'useodd',
        cl_data = cl_data(~even_cl_timesi,:);
      end
    end
    
      
    if any(cl_data),
      % we don't know in advance how many dimensions the histogram will
      % have, so we have to generate the index string dynamically
      colons = repmat({':'},1,ndim);
      tc(j,colons{:}) = histn (cl_data(:,2:end), ...
                               a.opt.histedges,...
                               a.opt.keepbins,...
                               'noedgebin');
      % else leave as zeros
    end

  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth tuning curves and occupancy hists 

% keep a copy of unsmoothed occ
occ(occ==0)=NaN;
  
% use Fabian's smoothn on both occ and spike counts
  
% set up args for smoothn (0 means no filtering)
ndim = length(a.opt.histedges);
for k = 1:ndim
  if isempty(a.opt.smooth_sds{k})
    filtsd(k) = 0;
  else
    filtsd(k) = a.opt.smooth_sds{k};
  end
  parm_dx(k) = mean(diff(a.opt.histedges{k}));
end
  
if any(filtsd)
  [occ kernel] = smoothn(occ, filtsd, parm_dx, 'nanexcl', 1, 'correct', 1);
  
  % dimension 1 is clusters, we don't want to smooth across that!
  [tc kernel] = smoothn(tc, [0 filtsd], [1 parm_dx], 'nanexcl', 1, 'correct', 1);
  
end

  
% $$$   
% $$$   % note HACK: we use filtfilt, which will apply the filter twice. To get
% $$$   % appropriate filter width, use requested sd*(sqrt(0.5)). (Determined
% $$$   % empirically to have same effect as applying full filter once). Magnitude
% $$$   % unaffected by applying filter twice, since fB always sums to 1 from
% $$$   % gausswinsd.
% $$$   
% $$$   tc_col = tc;
% $$$   
% $$$   for k = 1:ndims(tc)-1 % first dim is clusters
% $$$       % filter only operates columnwise, so shift each param into columns by
% $$$       % turns.
% $$$ 
% $$$       tc_col = shiftdim(tc_col,1);
% $$$ 
% $$$       if ~isempty(a.opt.smooth_sds{k}),
% $$$ 
% $$$         % calculate the sampling frequency for this parameter
% $$$         parm_fs = 1/mean(diff(a.opt.histedges{k}));
% $$$         % the sd size hack here (see above)
% $$$         fB = gausswinsd(a.opt.smooth_sds{k}*sqrt(0.5), parm_fs);
% $$$       
% $$$         tc_col = filtfilt(fB,1,tc_col);
% $$$         
% $$$       end
% $$$   end
% $$$   
% $$$   tc_col = shiftdim(tc_col,1);
% $$$   
% $$$   tc = tc_col;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % normalize position histogram by occupancy
  
  oldwarn = warning('off', 'Matlab:DivideByzero');
  tc = tc ./ repmat(shiftdim(occ,-1),ncls,1);
  warning(oldwarn);  

  % locations the rat never visted are NaNs, convert them to 0
  tc(isnan(tc)) = 0;
  %  tc(repmat(shiftdim(occ_unsmoothed==0,-1),ncls,1)) = 0;

  % convert to spikes/second by multiplying by position records/second
  tc = tc ./ posrecinterval;

  % convert occ to seconds
  occ = occ .* posrecinterval;
  


function [c cdat] = mkcorr(varargin)
% MKCORR make a 'corr' object: do an xcorr between 2 cdat chans, 2 cdat
% structs, or acorr on 1 chan/struct. Handles case of several segs within
% these cdats. Provides several different normalizations in output.
%
% -NaNs in data treated as missing values; correction takes account of this
% -overlapping/unordered segs OK.
  
% TODO:
% -deal with short windows? 
% -speed up calc of xc_correct/_count by not duplicating for events of
% same length
  
  c = struct(...
      'type', 'corr',...
      ...
  'cdats', [],... % inputs
      'samplerate', [],... % in samples/sec (empty means use rate from A)
      'refsig', '',...
      'maxlag_t', [],...
      'segs',[],...
      'segclip_t', [0 -0],... % clip segments (can also be single-valued)
      'seg_mindur', 0,...
      'xcov', false,... % calculate x-covariance (detrend = 'constant')
      'randperm', false,... % randomize datapoints within each seg
      ...% outputs
  'acorr', [],... % boolean
      'xc_raw', [],... % array of uncorrected xcorrs for each seg
      'xc_raw_mean', [],...
      'xc_raw_concat', [],...
      ...
      'raw_correct', [],... % # of comparisons per seg at each lag
      'raw_correct_sum', [],... % # of comparisons contributing at each lag
      'raw_correct_count_sum', [],... % # of segments contributing at each lag
      ...
      'xc_unbiased', [],... % array of unbiased xcorrs for each seg
      'xc_unbiased_mean', [],... % mean of unbiased xcorrs across segs 
      'xc_unbiased_concat', [],...% unbiased by # of comparisons
          ...
      'xc_coeff', [],... % array of 'coeff' normalized xcorrs (cf. 'xcorr')
      'xc_coeff_mean', [],...
      'xc_coeff_concat', [],...
          ...
      'xc_unbiasedcoeff', [],...
      'xc_unbiasedcoeff_mean', [],...
      'xc_unbiasedcoeff_concat', [],...
      'lags_t', []....
      );
  
  c = parseArgsLite(varargin,c);

  %%% set up/validate inputs

  if isempty(c.maxlag_t)
    if isempty(c.segs)
      error('if no segs, must provide ''maxlag_t''');
    end
    c.maxlag_t = max(diff(c.segs,[],2));
    warning(['no ''maxlag_t'' provided, using max seg duration: ' ...
             num2str(c.maxlag_t)]); 
  end
  
  %%% select data
  
  % find enclosing time window (saves time when combining/interpolating cdats)
  timewin = [max([c.cdats.tstart]) min([c.cdats.tend])];
  if isempty(c.segs),
    % if no segs given,assume corr over all overlapping range
    c.segs = timewin;
  else

    if any(c.segclip_t)
      if numel(c.segclip_t) == 1,
        c.segclip_t = c.segclip_t * [1 -1];
      end
      % clip segments
      c.segs = gplus(c.segs, c.segclip_t);
    end

    % remove any segments with no duration after clipping
    c.segs(diff(c.segs,[],2)< c.seg_mindur,:) = [];
    
    if isempty(c.segs)
      error('no segs that meet criteria');
    end
      
    % avoid clipping on resample, but don't make timewin larger than cdat overlap
    pad_t = 2./min([c.cdats.samplerate c.samplerate]);
    timewin = [max([timewin(1) min(c.segs(:))-pad_t]) ....
               min([timewin(2) max(c.segs(:))+pad_t])];
  end
  
  if diff(timewin)<=0, 
    error('no overlap between cdat inputs, segs; can''t compute xcorr');
  end

  switch length(c.cdats)
   case 1,
    cdat = contwin(c.cdats, timewin);
    if all(size(cdat.data,2) ~= [1 2])
      error('Exactly 1 or 2 channels must be specified (use contchans?)')
    end
   
    if ~isempty(c.samplerate),
      cdat = continterp(cdat, 'samplerate', c.samplerate);
    end
    
   case 2,
    if size(c.cdats(1).data,2) ~= 1 || size(c.cdats(2).data,2) ~= 1,
      error(['If 2 cdats are provided, they must have single data channels ' ...
             '(use contchans on inputs)']);
    end

    % Align samples by interpolation so that xcorr makes sense 

    % do it here, so that we don't trigger unnecessary interpolation
    c.cdats(1) = contwin(c.cdats(1), timewin);
    c.cdats(2) = contwin(c.cdats(2), timewin);

    if isempty(c.samplerate),
      % check for refsig
      switch c.refsig
       case 'A'
        c.samplerate = c.cdats(1).samplerate;
        timewin = [c.cdats(1).tstart c.cdats(1).tend];
       case 'B'
        c.samplerate = c.cdats(2).samplerate;
        timewin = [c.cdats(2).tstart c.cdats(2).tend];
       otherwise
        error('no samplerate or refsig provided');
      end
    end
    
    cdat = contcombine(c.cdats(1), c.cdats(2),...
                       'samplerate', c.samplerate,...
                       'timewin', timewin);
   otherwise
    error('Only 1 or 2 cdat structs can be cross-correlated');
  end
      
  if any(c.segs(:)>cdat.tend) || any(c.segs(:)<cdat.tstart),
    warning('some segs outside of cdat range, ignoring');
    c.segs = c.segs(inseg([cdat.tstart cdat.tend], c.segs));
  end
  
  
  %%% set up for xcorr
  
  % are we doing an acorr (1 channel) or an xcorr (2 channels)?
  c.acorr = size(cdat.data,2) == 1;
  
  % convert lags to samples
  maxlags = ceil(c.maxlag_t .* cdat.samplerate);
  
  nsegs = size(c.segs,1);
  
  % preallocate
  c.xc_raw = zeros(nsegs,2*maxlags+1);
  c.raw_correct = c.xc_raw;
  cAA0 = zeros(nsegs,1);
  cBB0 = cAA0;
  
  for k = 1:nsegs

    % get data for this seg
    cdat_win = contwin(cdat,c.segs(k,:));
    A = cdat_win.data(:,1);
    if c.randperm
      A = A(randperm(numel(A)));
    end
    
    if ~c.acorr % xcorr
      B = cdat_win.data(:,2);
      if c.randperm
        B = B(randperm(numel(B)));
      end

      if c.xcov
        A = A - nanmean(A);
        B = B - nanmean(B);
      end

      % how many comparisons at each lag for this seg?
      c.raw_correct(k,:) = ...
          c.raw_correct(k,:) + ...
          round(xcorr(~isnan(A),~isnan(B), maxlags, 'none')'); % sometimes ~=0
      
      % NaN treated as missing values: don't increment correction
      % counter, no correlation since value = 0;
      A(isnan(A)) = 0;
      B(isnan(B)) = 0; 
      
      [c.xc_raw(k,:) lags] = xcorr(A, B, maxlags, 'none');

      % Compute autocorrelations at zero lag for 'coeff' scaling below (borrowed
      % from xcorr). NaNs are zeroed, so don't contribute to sum
      cAA0(k) = sum(abs(A).^2);
      cBB0(k) = sum(abs(B).^2);
      
    else % acorr

      if c.xcov
        A = A - nanmean(A);
      end

      % NaN treated as missing values: don't increment correction
      % counter, no contribution to correlation since value = 0;
      A(isnan(A)) = 0;
      
      [c.xc_raw(k,:) lags] = xcorr(A, maxlags, 'none');

      c.raw_correct(k,:) = ...
          c.raw_correct(k,:) + ...
          round(xcorr(~isnan(A), maxlags, 'none')'); % sometimes ~=0
      
    end
    
  end

  %% post-process various: correction strategies, etc

  lag0 = maxlags+1;
  
  % # of comparisons at each time lag
  c.raw_correct_sum = sum(c.raw_correct, 1);
  
  % # of segments contributing at each time lag
  c.raw_correct_count_sum = sum(c.raw_correct>0, 1);
  
  %%% 'raw', (AKA Matlab's 'none')
  
  % raw values are already in xc_raw
  
  % mean of per-segment raw estimates
  c.xc_raw_mean = nanmean(c.xc_raw,1);
  
  % concat not meaningful (would be same as mean)
  
  %%% 'unbiased'

  % unbiased estimate for each segment
  c.xc_unbiased = c.xc_raw ./ c.raw_correct;
  c.xc_unbiased(isinf(c.xc_unbiased)) = NaN;
  
  % mean of unbiased estimates
  c.xc_unbiased_mean = nanmean(c.xc_unbiased,1);
  
  % unbiased estimate as if one long signal with only within-seg
  % comparisons used
  c.xc_unbiased_concat = nansum(c.xc_raw,1)./c.raw_correct_sum;

  
  %%% 'coeff'
  % (like Matlab: normalize lag0 to corr coeff, but don't take account of # of
  % comparisons)

  if c.acorr
    % value at t=0 is 1 for autocorrs
    scale = c.xc_raw(:,lag0);
  else
    % borrowed from matlab xcorr
    scale = sqrt(cAA0.*cBB0);
  end

  % coeff estimate for each segment
  c.xc_coeff = grdivide(c.xc_raw, scale);

  % mean of coeff estimates
  c.xc_coeff_mean = nanmean(c.xc_coeff,1);
  
  % coeff estimate as if one long signal
  if c.acorr
    scale = sum(c.xc_raw(:,lag0));
  else
    scale = sqrt(sum(cAA0)*sum(cBB0));
  end
  c.xc_coeff_concat = sum(c.xc_raw,1)./scale;
  
  
  %%% 'unbiasedcoeff' 
  % unbiased estimate, with value at lag = 0 scaled to be corr coeff

  if c.acorr
    % value at t=0 is 1 for autocorrs
    scale = c.xc_unbiased(:,lag0);
  else
    % unbiased has already corrected for # of comparisons at lag0, so uncorrect
    scale = sqrt(cAA0.*cBB0) .* c.raw_correct(:,lag0);
  end
  c.xc_unbiasedcoeff = grdivide(c.xc_unbiased, scale(:));

  % mean of unbiased estimates
  c.xc_unbiasedcoeff_mean = nanmean(c.xc_unbiasedcoeff,1);

  % unbiasedcoeff_concat
  if c.acorr
    c.xc_unbiasedcoeff_concat = c.xc_unbiased_concat ./ c.xc_unbiased_concat(lag0);
  else
    c.xc_unbiasedcoeff_concat = c.xc_unbiased_concat ./ ...
        sqrt(sum(cAA0)*sum(cBB0)) * c.raw_correct_sum(lag0);
  end
  

  % calculate lags in seconds
  c.lags_t = lags ./ cdat.samplerate;
  
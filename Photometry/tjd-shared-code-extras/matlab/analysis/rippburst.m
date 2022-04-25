function [s cache] = rippburst(varargin)
% RIPPBURST detect ripples, bursts
  
% todo:
%  -low-pass filter for envelope rather than gausswin to be more explicit?
  
  a = struct(...
      'contdata', [],...
      'cdat_ripp', [],...
      'cdat_rippenv', [],...
      'envmean', [],... % use mean of envelope across many channels?
      'contvar', [],...
      'chans', [],...
      'chanlabels', [],...
      'timewin',[],...
      'ripp_F', [],...
      'filtopt', [],...
      'bouts',[],... % *stopped bouts*, used for determining mean/std thresh
      'env_method', 'hilbert', ...
      'rms_window_t', [],...
      'ripp_method', 'localmax',...
      'smooth_type', 'gausswin',...
      'smooth_sd_t', [],...
      'smooth_F', [],...
      'thresh_std', [],...
      'thresh_lo_std', [],...
      'thresh_hi_std', [],...
      'useSW', false,... % also apply a SW-band env threshold
      'SW_F', [],... % SW-bandpass frequencies (1x4)
      'SWenv_thresh_std', [],... % stdev above mean for SW amplitude
      'burst_iri', 0.1,...
      'cache', []);
  
  a = parseArgsLite(varargin,a);

  % save args, but not raw data
  s.args = a;
  s.args.contdata = [];
  s.args.cdat_ripp = [];
  s.args.cdat_rippenv = [];
  
  if ~isempty(a.smooth_sd_t) && strcmp(a.env_method, 'rms'),
    warning(['RMS envelope method *and* post-env smooth requested -- did ' ...
             'you mean this?']);
  end

  if ~isempty(a.timewin),
    error('''timewin'' not currently supported--use bouts/seginseg?');
  end
  a.timewin = [-Inf Inf];
  
  if isempty(a.bouts)
    warning('Searching for ripples in whole experiment');
    a.bouts = [-Inf Inf];
  end
  

  %% get unfiltered, ripple band filtered, and enveloped ripple signal, from
  %% cache if available

  %cache = a.cache;
  cache = [];
  
  % design ripple filter
  if ~isempty(a.filtopt),
    ripp_fopt = a.filtopt;
  else
    ripp_fopt = mkfiltopt('name', 'ripple', ...
                          'filttype', 'bandpass',...
                          'F', a.ripp_F);
  end
  
  % design env options
  envopt = mkenvopt('method', a.env_method,...
                    'rms_window_t', a.rms_window_t);
  
  contopt = mkcontopt();
  
  contopt_ripp = mkcontopt('filtopt', ripp_fopt);

  contopt_rippenv = mkcontopt('filtopt', ripp_fopt,...
                              'envopt', envopt);
  
  cont = mkcont ('contvar', a.contvar,...
                 'contdata', a.contdata,...
                 'chans', a.chans,...
                 'chanlabels', a.chanlabels,...
                 'contopt', contopt,...
                 'timewin', a.timewin,...
                 'cache', cache);
  
  %cache = mkcache('cache', cache, 'add_obj', cont);
  cache = [];
  
  [cont_ripp s.filt_ripp] = mkcont(...
      'contvar', a.contvar,...
      'contdata', a.contdata,...
      'chans', a.chans,...
      'chanlabels', a.chanlabels,...
      'contopt', contopt_ripp,...
      'timewin', a.timewin,...
      'cache', cache);
  
  % cache = mkcache('cache', cache, 'add_obj', cont_ripp);
  % cache = mkcache('cache', cache, 'add_obj', s.filt_ripp);
  cache = [];
  
  cont_rippenv = mkcont('contvar', a.contvar,...
                        'contdata', a.contdata,...
                        'chans', a.chans,...
                        'chanlabels', a.chanlabels,...
                        'contopt', contopt_rippenv,...
                        'timewin', a.timewin,...
                        'cache', cache);
  
  %cache = mkcache('cache', cache, 'add_obj', cont_rippenv);
  cache = [];

  cdat = cont.contdata;
  cdat_ripp = cont_ripp.contdata;
  cdat_rippenv = cont_rippenv.contdata;  

  % how many channels are we dealing with?
  nchans = size(cdat_rippenv.data,2);
  
  % hack in combined channel (only if more than one chan provided)
  if nchans >1, 
    cdat_rippenvmean = contmean(cont_rippenv.contdata);
    % no resampling should be necessary
    cdat_rippenv = contcombine(cdat_rippenv, {cdat_rippenvmean}, 'match_first', true);
    envmeanlast = true;
    nchans = nchans +1;
  else
    envmeanlast = false;
  end
  
  %%% CUT
  if ~isempty(a.cdat_ripp)
    cdat_ripp = a.cdat_ripp;
    cdat_rippenv = a.cdat_rippenv;
  end
    
  % save for output
  s.cdat = cdat;
  s.cdat_ripp = cdat_ripp;

  if envmeanlast
    s.cdat_rippenv = contchans(cdat_rippenv, 'chans', 1:nchans-1);
    s.cdat_rippenv_mean = contchans(cdat_rippenv, 'chans', nchans);
  else
    s.cdat_rippenv = cdat_rippenv;
    s.cdat_rippenv_mean = [];
  end
  
  switch a.smooth_type,
   case {'gausswin','hatwin'},
    envsmooth_fopt = mkfiltopt('name', ['rippenv_smooth_' a.smooth_type],...
                               'filttype', a.smooth_type, ...
                               'sd_t', a.smooth_sd_t);
    
   case {'rectwin'},
    envsmooth_fopt = mkfiltopt('name', ['rippenv_smooth_' a.smooth_type],...
                               'filttype', a.smooth_type, ...
                               'length_t', a.smooth_sd_t);    
   case 'lowpass'
    envsmooth_fopt = mkfiltopt('name', 'rippenv_smooth_lowpass',...
                               'filttype', 'lowpass', ...
                               'F', a.smooth_F);
    
   case 'none'
    envsmooth_fopt = [];
   otherwise
    error('unrecognized ''smooth_type''');
    
  end
  envsmooth_contopt = mkcontopt('filtopt', envsmooth_fopt);
  
  if ~isempty(envsmooth_fopt)
    cont_rippenv_smooth = mkcont('contdata', cdat_rippenv,...
                                 'chans', [],...
                                 'chanlabels', [],...
                                 'contopt', envsmooth_contopt,...
                                 'timewin', a.timewin,...
                                 'cache', cache);
    
    %    cache = mkcache('cache', cache, 'add_obj', cont_rippenv_smooth);
    cache = [];
    
    cdat_rippenv_smooth = cont_rippenv_smooth.contdata;
    if envmeanlast
      s.cdat_rippenv_smooth = contchans(cdat_rippenv_smooth, 'chans', 1:nchans-1);
      s.cdat_rippenv_smooth_mean = contchans(cdat_rippenv_smooth, 'chans', nchans);
    else
      s.cdat_rippenv_smooth = cdat_rippenv_smooth;
      s.cdat_rippenv_smooth_mean = [];
    end
    
  else
    s.cdat_rippenv_smooth = cdat_rippenv;
  end
  
  switch(a.ripp_method)
   case {'localmax', 'localmax2'}
    % get local maxima of smoothed ripple envelope (one cell per column)
    [lindex s.ripps_all] = contlocalmax(cdat_rippenv_smooth);
    for k = 1:nchans,
      rippenvf_data = contsegdata(contchans(cdat_rippenv_smooth,'chans', k), a.bouts);
      s.rippenvf_mean(k) = mean(rippenvf_data(:));
      s.rippenvf_std(k) = std(rippenvf_data(:));
      s.rippenvf_thresh(k) = s.rippenvf_mean(k) + (a.thresh_std .* s.rippenvf_std(k));
      
      % save out signal stats
      s.thresh(k) = s.rippenvf_thresh(k);
      
      if ~isempty(s.ripps_all{k})
          ripp_times_all = s.ripps_all{k}(:,1);
          %ripp_indexes = s.ripps_all{k}(:,2);
          ripp_vals = s.ripps_all{k}(:,3);
      else
          ripp_times_all = false(0,1);
          %ripp_indexes = false(0,1);
          ripp_vals = false(0,1);
      end
      
      inbouts_i = ...
          inseg(a.bouts, ripp_times_all);
      
      aboveth_i = ...
          ripp_vals > s.rippenvf_thresh(k);
      
      s.ripps_goodi{k} = inbouts_i & aboveth_i;
      
      s.ripps_good{k} = s.ripps_all{k}(s.ripps_goodi{k},:);
      if ~isempty(s.ripps_good{k}),
        s.ripp_times{k} = s.ripps_good{k}(:,1);
      else
        s.ripp_times{k} = [];
      end
      
    end
    
    
   case 'thresh'
    % ripple is at mean of ripple start/end times (threshold crossing times)
    for k = 1:nchans,
      % take stdev of smoothed data (mean should be 0);
      bout_data = contsegdata(s.cdat_rippenv_smooth, a.bouts);
      s.rippenvf_mean(k) = mean(bout_data(:,k));
      s.rippenvf_std(k) = std(bout_data(:,k));
      s.rippenvf_thresh_lo(k) = s.rippenvf_mean(k) + (a.thresh_lo_std .* s.rippenvf_std(k));
      s.rippenvf_thresh_hi(k) = s.rippenvf_mean(k) + (a.thresh_hi_std .* s.rippenvf_std(k));

      s.ripp_times{k} = contbouts(s.cdat_rippenv_smooth,...
                                  'datargunits', 'data',...
                                  'thresh', s.rippenvf_thresh_lo(k),...
                                  'minevpeak',s.rippenvf_thresh_hi(k),...
                                  'window', 0.0,...
                                  'mindur', 0.01);
      
      s.ripps_good{k} = mean(s.ripp_times{k},2);
    end
     
   otherwise
    error('unrecognized ''ripp_method''');
    
  end

  %%% apply Sharp-wave amplitude criterion, if requested
  if a.useSW,
% $$$     'SWenv_thresh_std', [],... % stdev above mean for SW amplitude
% design ripple filter
    
    SW_fopt = mkfiltopt('name', 'SW', ...
                        'filttype', 'bandpass',...
                        'F', a.SW_F);
    
    % design env options
    envopt = mkenvopt('method', 'hilbert');
    contopt_SW = mkcontopt('filtopt', SW_fopt);
    contopt_SWenv = mkcontopt('filtopt', SW_fopt,...
                              'envopt', envopt);

    cache = [];
    [cont_SW s.filt_SW] = mkcont(...
        'contvar', a.contvar,...
        'contdata', a.contdata,...
        'chans', a.chans,...
        'chanlabels', a.chanlabels,...
        'contopt', contopt_SW,...
        'timewin', a.timewin,...
        'cache', cache);

    cache = [];
    cont_SWenv = mkcont('contvar', a.contvar,...
                        'contdata', a.contdata,...
                        'chans', a.chans,...
                        'chanlabels', a.chanlabels,...
                        'contopt', contopt_SWenv,...
                        'timewin', a.timewin,...
                        'cache', cache);
                      
                      
    cache = [];
    cdat_SW = cont_SW.contdata;
    cdat_SWenv = cont_SWenv.contdata;  
    
    warning('Hacked in SWenv smoothing using same params as rippenv');
    cdat_SWenv = contfilt(cdat_SWenv,'filtopt', envsmooth_fopt);
    
    if envmeanlast
      cdat_SWenvmean = contmean(cont_SWenv.contdata);
      % no resampling should be necessary
      cdat_SWenv = contcombine(cdat_SWenv, {cdat_SWenvmean}, 'match_first', true);
    end
    
    % save for output
    s.cdat_SW = cdat_SW;
    
    if envmeanlast
      s.cdat_SWenv = contchans(cdat_SWenv, 'chans', 1:nchans-1);
      s.cdat_SWenv_mean = contchans(cdat_SWenv, 'chans', nchans);
    else
      s.cdat_SWenv = cdat_SWenv;
      s.cdat_SWenv_mean = [];
    end

    for k = 1:nchans,
      %%% get mean/std of SW
      cdat_SWenv_k = contchans(cdat_SWenv,'chans',k);

      SWenv_data = contsegdata(cdat_SWenv_k, a.bouts);
      s.SWenv_mean(k) = mean(SWenv_data(:));
      s.SWenv_std(k) = std(SWenv_data(:));
      s.SWenv_thresh(k) = s.SWenv_mean(k) + (a.SWenv_thresh_std .* s.SWenv_std(k));

      %% look up SWenv value at ripp times
      ripp_times_all = s.ripps_all{k}(:,1);
      ripp_SWenv_vals = contlookup(cdat_SWenv_k, 'xt', ripp_times_all, ...
                                   'method', 'linear');

      s.ripps_all{k}(:,4) = ripp_SWenv_vals;
      
      aboveSWth_i = ...
          ripp_SWenv_vals > s.SWenv_thresh(k);

      s.ripps_goodi_SWth{k} = s.ripps_goodi{k} & aboveSWth_i;
      s.ripp_times_SWth{k} = s.ripps_all{k}(s.ripps_goodi_SWth{k},1);
    end    
  end
  
  for k = 1:nchans, 
    %% detect ripple bursts
    % each row of s.bursts{k} has burst start/end time, and # of ripps
    s.bursts{k} = [];
    
    % is next ripple within iri?
    if numel(s.ripps_good{k})>2
      burst_i = diff(s.ripps_good{k}(:,1))<a.burst_iri;
    else
      burst_i = [];
      s.bursts{k} = [];
      break
    end
      
    blen = 0;
    last_zero_i = 0;
    
    for m = 1:(length(s.ripps_good{k})-1)
      if burst_i(m)
        % == 1: we are in a burst: increment burst length counter, defer writing to
        % s.bursts{k} until we know how long the burst is
        blen = blen+1;
      else
        % == 0: burst is over/never started
        
        if blen >0
          % end of a burst, write it out
          s.bursts{k} = [s.bursts{k} ;
                         s.ripps_good{k}([last_zero_i+1 m],1)' blen+1];
          
        else
          
          % burst_i == 0 is also always a 'burst of one' ripple
          s.bursts{k} = [s.bursts{k} ;
                         s.ripps_good{k}([m m],1)' 1];
        end
        
        blen = 0;
        last_zero_i = m; % potential start  
      end
    end % burst finding loop
    
  end
  
  %% get ripple waveform peak nearest to smoothed amplitude peak (so that
  %%  we can average the signals later)
  
  for k = 1:nchans-1, % not valid for mean signal?

    % get every peak in the ripple-filtered lfp for this channel
    ripp_peaks_t = ((find(localmax(cdat_ripp.data(:,k)))-1) ./ ...
                    cdat_ripp.samplerate) + cdat_ripp.tstart;
    
    % find the peak time align this channel's rippenv peak 
    s.ripp_times_pks{k} = ripp_peaks_t(findnearest(ripp_peaks_t, ...
                                                   mean(s.ripp_times{k},2), ...
                                                   'sorted'));
    
    % find this channel's ripple-peak time nearest the peak in the *mean* ripple
    % envelope
    s.ripp_times_mean_pks{k} = ripp_peaks_t(findnearest(ripp_peaks_t,...
                                                      mean(s.ripp_times{end},2),...
                                                      'sorted'));

  end
  
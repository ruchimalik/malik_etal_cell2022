function [pmap] = ...
    doparmest(varargin)
% DOPARMEST - reconstruct a parameter from tuning curves, spiking
%
% INPUTS:
% clh: a histogram of cluster activity (i.e. spiking)
% histbins: the bin edges used to make the histogram
% tcdats: the tuning curves to be used
% occ: the animal's occupancy across the parameter to reconstruct
% clperm: a permutation of the cluster->tuning curve mapping
% opt: the 'parmestopt' structure with the reconstruction options
%
% OUTPUTS:
% pmap = 2D estimate through time (estimates in columns)
% See mkposestopt for reconstruction options

% todo:
%
% -exp term in zhang will break with non-uniform timebinsize
% -pass in clhists to save time? (use calling function with a.opt)
% -histbins hack OK? (since array is sorted)
%
% 'tc' = tuning curve (usu place-rate map)

a = struct('clh', [],...
           'histbins',[],...
           'tcdats', [], ...
           'occ', [], ...
           'clperm',[],...
           'opt',[]); % see mkparmestopt for options

a = parseArgsLite(varargin, a);

ncls = size(a.clh.clhists,1);

if ncls ~= 0 &&  isempty(a.tcdats)
  error('Must provide tcdats, use ''tunecurve''');
end

% to calculate n-dimensional tcdats, just convert each tc to a row vector at
% the beginning, then reshape at the end (since each bin is calculated
% independently)
tcsz = size(a.tcdats);
a.tcdats = reshape(a.tcdats,ncls,[],1);

if isempty(a.opt.bsmallf),
  a.opt.bsmallf = 0.01;
end

% reconstruction method
switch a.opt.method,
 case 'basis',
  % fn to use to combine outputs of cltc_comb, below
  clcl_fn = @plus;
  % fn to use to initialize the pmap
  mapinit_fn = @zeros;
 case 'bayes',
  clcl_fn = @times;
  mapinit_fn = @ones;
 otherwise
  error('requested ''method'' not supported');
end

if isempty(a.tcdats)
  pmap = mapinit_fn(prod(tcsz(2:end)),a.clh.nbins);
  pmap = pmap ./ prod(tcsz(2:end));
  tcsz = num2cell(tcsz);
  pmap = reshape(pmap',[],tcsz{2:end});
  warning('No clusters selected for parmest--returning uniform PDF');
  return
end

% number of columns (bins) in tuning curves; rows in pmap
[ntcs nrows] = size(a.tcdats);

% indexes of clusters with any non-zero response in tuning curves
tcanyi = find(any(a.tcdats,2))';


% $$$ if a.opt.ratenorm,
% $$$   % normalize spiking by each cell's mean firing rate during run
% $$$   for j = tcanyi,
% $$$     a.clh.clhists(j,:) = a.clh.clhists(j,:) ./ mean(tcs(j,:));
% $$$   end
% $$$ end

%%% preallocate output map(s)

% zeros or ones, depending on method
pmap = mapinit_fn(nrows,a.clh.nbins);

% no clusters, return empty arrays.
if ncls == 0,
  return;
end

tcs = double(a.tcdats);

if a.opt.pdfnorm
  
  % this is not really a pdf of spiking. Spiking is point
  % process. Likelihood of getting a certain # of spikes given a rate map
  % (and assuming a poisson process) is a different beast.
  
  % Note we mean-normalize so that the values don't get too small (float
  % underrun) with large numbers of clusters. Sensible?

  for j = tcanyi; % row-wise division loop
    tcs(j,:) = tcs(j,:) / mean(tcs(j,:));
  end
  
end

% term by which to multiply tuning curve magnitude
% (needs to go after pdfnorm)
if ~isempty(a.opt.tcfudgef) && all(a.opt.tcfudgef ~= false),
  if length(a.opt.tcfudgef) == 1,
    tcs = tcs .* a.opt.tcfudgef;
  else
    for j = 1:ntcs,
      tcs(j,:) = tcs(j,:) .* a.opt.tcfudgef(j);
    end
  end
end

% add a small value to tcs to avoid problem of 'zeroing' out estimate
% not necessary ?
if strcmp(a.opt.method, 'bayes'),
  bsmall = max(tcs(:)) * a.opt.bsmallf;
  tcs(tcs < bsmall) = bsmall;
end

% if requested, scale each tc pdf by the cell's 'place' spatial info score
if a.opt.infoscale,
  pinfo = spatialinfo(tcs);
  for j = 1:ncls; % row-wise mult loop
    tcs(j,:) = tcs(j,:) * pinfo(j);
  end
end

% after we've calc'd/scaled according to 'real' pdfs, we can now permute the
% cl->tc mapping, if requested
if ~isempty(a.clperm),
  tcs = tcs(a.clperm,:);
end

% only calculate map using cells with any tuning curve response
for j = 1:length(tcanyi),
  
  % for efficiency, and to avoid multiplying through by 0 in 'bayes'
  % (wow--speeds up run by ~50x in many cases)
  nonz = a.clh.clhists(j,:) ~= 0;
  
  % build the 2D histogram
  if any(a.clh.clhists(j,:)) && any(a.tcdats(j,:)), % since tcs is > 0 due to bsmallf
    % sum/multiply the products/power of firing rate and place-rate map
    pmap(:,nonz) = clcl_fn(pmap(:,nonz),...
                           cltc_comb(a.opt.method, ...
                                     tcs(j,:),...
                                     a.clh.clhists(j,nonz)));
  end
  
end

%%% normalizations/extra terms
rownorm = ones(nrows,1);
dorownorm = false;

elnorm = 1;
doelnorm = false;

colnorm = ones(1,a.clh.nbins);
docolnorm = false;


if a.opt.expnorm,
  % the exp term from eq 36 of Zhang et al 1998 J. Neurophys
  tcsum = sum(tcs)';

% $$$   % avoid exp(zero)
% $$$   tcsum(~tcsum) = -Inf;
  
  if length(a.clh.histbinsize) == 1 || ~any(diff(a.clh.histbinsize))
    dorownorm = true;
    rownorm = rownorm .* exp(-a.clh.histbinsize(1) * tcsum);
    rownorm_zeroi = (rownorm == 0);
    if any(rownorm_zeroi)
      warning(['float underrun in ''exp'' term--parmest unreliable--try ' ...
               'smaller time window']);
      rownorm(rownorm_zeroi) = realmin;
    end
  else 
    % if histbinsize is non-uniform (i.e. if we have given histbins),
    % then we can't just do a columnnorm/rownorm
    arrnorm = ones(nrows,a.clh.nbins);
    for k = 1:a.clh.nbins
      arrnorm(:,k) = exp(tcsum);
    end
    for j = 1:nrows,
      arrnorm(j,:) = arrnorm(j,:) .^ -a.clh.histbinsize;
    end
    pmap = pmap .* arrnorm;
  end
end

%%% test for float underruns/overruns
if any(isnan(pmap(:))),
  warning(['un-normalized parmest map contains NaN''s -- estimate unreliable ' ...
           '-- try smaller timewin']);
end

if any(isinf(pmap(:))),
  warning(['un-normalized parmest map contains Inf''s -- estimate unreliable ' ...
           '-- converting to realmax -- try smaller timewin']);
  pmap(isinf(pmap)) = realmax;
end


%%% row-wise normalizations:

if a.opt.occnorm,
  dorownorm = true;
  
  % convert to P(x) - unnecessary if scaling later, but whatevs.
  a.occ = a.occ ./ sum(a.occ);
  a.occ(~a.occ) = Inf;
  rownorm = rownorm ./ a.occ(:);

end

if a.opt.basisnorm,
  dorownorm = true;
  
  % normalize maps by sum of basis functions (pdfs or raw dats? dats looks way
  % nicer, does it also include a rate normalization?) to correct for 'bias' -
  % maybe more valid if you want to interpret output as an estimate?
  
  tcsum = sum(a.tcdats)';

  % mean normalize so that the scaling can stay the same
  tcsum = tcsum./mean(tcsum);

  % avoid div by zero -> nan's
  tcsum(~tcsum) = Inf;

  % divide each row by average
  rownorm = rownorm / tcsum;

end

% do row-wise normalization:
if dorownorm,
  for k = 1:nrows,
    pmap(k,:) = pmap(k,:) * rownorm(k);
  end
end


%%% element-wise (scalar) normalizations

if a.opt.infoscale,
  doelnorm = true;
  
  % normalize maps by mean spatial info
  elnorm = elnorm / (sum(pinfo) / ncls);
end

if strcmp(a.opt.method,'basis');
  doelnorm = true;
  
  % normalize by # of cells used to make estimate
  elnorm = elnorm /  ncls;
end

if doelnorm,
  pmap = pmap.*elnorm;
end


%%% Column-wise normalizations

% do last because columnnorm sums across existing maps

if a.opt.columnnorm,
  docolnorm = true;
  % make column sum to 1 (pdf) (add eps to avoid small values >1)
  colnorm = colnorm ./ sum(pmap);
end

% these are weird options that attempt to use # of cells/spikes in each
% bin as a stand-in for 'confidence', and skew the visualization accordingly.

if a.opt.cellcountnorm,
  docolnorm = true;
% $$$   sharpenf = 0.5;
  colnorm = colnorm .* cellcount.^sharpenf;
end

if a.opt.spikecountnorm
  docolnorm = true;
% $$$   sharpenf = 0.5;
  colnorm = colnorm .* spikecount.^sharpenf;
end

% do column-wise normalization
if docolnorm,
  colnorm = colnorm + eps(colnorm);
  for k = 1:a.clh.nbins, % column-wise division loop
    pmap(:,k) = pmap(:,k).* colnorm(k);
  end
end



%%% convert back to nd array
% time in first dimension: n timebin rows by 100 x 2, e.g.
tcszc = num2cell(tcsz);
pmap = reshape(pmap',[],tcszc{2:end});

%% history term 
if ~isempty(a.opt.history_gausswin_sd_bins)
  warning('Woop--untested history stuff');
  % implement history term from Zhang 1998 J. Neurophys
  % -prev estimate * empirically determined gaussian
  % -skip first timebin as no prev est is available
  % -convert pmap back to N-D first
  tcszc = num2cell(tcsz);
  if numel(a.opt.history_gausswin_sd_bins) ~= numel(tcszc)-1
    warning('wrong ndims for history gausswin');
  end
  
  arrnorm = smoothn(pmap,...
                    [0 a.opt.history_gausswin_sd_bins],...
                    'correct', 1, 'nsds', 8);
  pmap(2:end,:) = pmap(2:end,:) .* arrnorm(1:end-1,:);

  % make column sum to 1 (pdf) (add eps to avoid small values >1)
  colnorm = ones(1,a.clh.nbins);
  colnorm = colnorm ./ sum(pmap');
  colnorm = colnorm + eps(colnorm);
  for k = 1:a.clh.nbins, % column-wise division loop
    pmap(k,:) = pmap(k,:).* colnorm(k);
  end

end



%%% SUBFUNCTIONS

function map = cltc_comb(method, tc, clhist)
% subfunction to combine a spike histogram and a place-rate map

% both tc and clhist are passed in as row vectors
  
  switch method
    
   case 'basis',
    map = tc' * clhist;
   
   case 'bayes',
    % there is no matrix-matrix power operation
    map = repmat(tc', 1, length(clhist)) .^ ...
          repmat(clhist, length(tc), 1);
  
  end
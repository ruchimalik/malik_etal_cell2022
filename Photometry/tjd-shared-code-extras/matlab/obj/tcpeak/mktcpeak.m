function [pks tc] = mktcpeak(varargin)
% MKTCPEAK - make a tuning curve peaks structure
  
  pks = struct(...
      'type', 'tcpeak',...
      ...
      'e', [],... % tc inputs
      'edesc', [],...
      'clnos', [],...
      'tunecurveopt',[],... 
      'tcpeakopt',[],...
      'dim', 1,...
      ...
      'pkdats',{{}},... % outputs
      ...
      'template', [],... % caching, etc
      'cache', [],...
      'cache_hit', false,...
      ...
      'plot',false);


% $$$   if isempty(pks.peakfindgw),
% $$$     error('no ''peakfindgw'' specified');
% $$$   end
% $$$   
% $$$   if isempty(pks.peakfindth),
% $$$     error('no ''peakfindth'' specified');
% $$$   end

      
  % handle 'template' object arg
  pks = obj_reparse(pks,varargin);

  if islogical(pks.clnos),
    pks.clnos = find(pks.clnos);
  end
  
  pks = obj_cachesearch(pks);
  
  if pks.cache_hit,
    pks = obj_cleanup(pks);
    tc = [];
    return;
  end
  
  %%%% no cache hit, make new peaks

  pks.pkdats = [];

  tc = mktunecurve('e', pks.e,...
                   'edesc', pks.edesc,...
                   'clnos', pks.clnos,...
                   'tunecurveopt', pks.tunecurveopt,...
                   'cache', pks.cache);

    
  % tcdats may be a superset; select tuning curves
  ncls = length(pks.clnos);
  [dummy tcdatsi] = ismember(pks.clnos,tc.clnos); %#ok
  tcdats = tc.tcdats(tcdatsi,:);
  
  for j = 1:ncls;
    % get pf centers
    [pks_j tcdatsf_j] = subf_findpeaks(tcdats(j,:),...
        pks.tcpeakopt.peakfindgw,...
        pks.tcpeakopt.peakfindth,...
        tc.tunecurveopt.histedges{pks.dim});
    
    pks.pkdats=[pks.pkdats; {pks_j}];
    
    if pks.plot,
      figure;
      subf_peaksplot(tcdats(j,:), pks_j, tcdatsf_j, tc.tunecurveopt.histedges{a.dim}, pks.peakfindth);
    end
  end
  
  pks = obj_cleanup(pks);

  
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SUBFUNCTIONS

function [pks nf] = subf_findpeaks(n,gw,th,histedges)
  
pks = [];

% smooth that ass...
if ~isempty(gw)
  fa = gausswinsd(gw,1/mean(diff(histedges)));
  nf = filtfilt (fa, sum(fa), n);
else
  nf = n;
end

% find peaks
[dummy pks_val] = localmax(nf');
pks = pks_val{1};

if isempty(pks),
  return;
end

% lookup peak
histctrs = histedges(1:end-1) + diff(histedges);
pks(:,1) = histctrs(pks(:,1));

% select peaks above threshold
pks = pks(pks(:,2) > th, :);

% reorder peaks according to 'strength' (height of filtered signal at
% peak)
if ~isempty(pks);
  pks = sortrows(pks,-2);
end

% $$$ 
% $$$ % approximate derivatives
% $$$ nfd = diff(nf);
% $$$ 
% $$$ % find zero crossings (where smoothed signal is above threshold).
% $$$ % Also get value at peak
% $$$ for j = 2:size(nfd,2)-1
% $$$   if nfd(j) < 0 && nfd(j-1) >= 0 && nf(j) > th
% $$$     pks = [pks; mean(histedges(j:j+1)) nf(j)];
% $$$   end
% $$$ end
% $$$ 
% $$$ if ~isempty(pks);
% $$$   pks = flipud(sortrows(pks,2));
% $$$ end


function h = subf_peaksplot(tcdats, peaks, tcdatsf, histedges, th)

h = gca;

histbinctrs = mean([histedges(1:end-1); histedges(2:end)]);

% $$$ % last bin no use
% $$$ tcdats(end) = [];
% $$$ tcdatsf(end) = [];

% forward bar plot (blue)
barh = bar  (histbinctrs, tcdats);
set(barh, 'facecolor', 'b');
set(barh, 'edgecolor', 'w');

% plot smoothed 
plot (histbinctrs, tcdatsf, 'k--');

 maxy = max([tcdats tcdatsf]);
 if maxy 
   ylim([0 1.5.*maxy]);
 end

% $$$ ylim([-200 200]);
% $$$ maxy = 25;

if ~isempty(peaks),
  for j = peaks(:,1)'
    plot([j j], [1.4*maxy 1.45*maxy ], 'k', 'linewidth', 3);
    plot([j j], [0 1.2*maxy], 'k--', 'linewidth', 0.5 );
  end
end

plot([histedges(1) histedges(end)], [th th], 'm:');






  
  
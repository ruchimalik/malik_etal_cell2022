function vci = tcvalidcells(varargin)
% TCVALIDCELLS Select clusters according to mean/max firing rate criteria (logical)
%
% Warning: assumes clusters are in first dimension, linear position in
% next dimension.
  
  a = struct(...
      'tc', [],...
      'min_peakfr', 0,...
      'max_meanfr', Inf);
  
  a = parseArgsLite(varargin,a);

  % default use all cls
  vci = true(1,length(a.tc.clnos));
  
  if ~isempty(a.max_meanfr),
    vci = vci & mean(a.tc.tcdats(:,:),2)' < a.max_meanfr;
  end
  
  if ~isempty(a.min_peakfr)
    vci = vci & max(a.tc.tcdats(:,:),[],2)' > a.min_peakfr;
  end

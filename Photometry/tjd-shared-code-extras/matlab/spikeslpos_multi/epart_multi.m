function newe = epart_multi(varargin)
% IM_PARTITION - create exp with partitioned clusters from exp with parm clusters
  
  a = struct('e', [],...
             'timewin', [-Inf Inf],...
             'uvthresh', '0',...
             'r_splits', [0 100],... % default is no splits
             'ang_splits', [0 100]);
  
  a = parseArgsLite(varargin,a);

  newe = a.e;

  neps = numel(a.e.epoch);
  for epi = 1:neps
    newe.epoch(epi).cl = [];
  end
  
  %%%  partitions each parm cl
  for epi = 1:neps
    disp(['partitioning tets from epoch : ' a.e.epoch(epi).name]);
    for k = 1:length(a.e.epoch(epi).cl)
      
      disp(['partitioning parms: ' a.e.epoch(epi).cl(k).name]);
      
      part_cls = clpart('cl', a.e.epoch(epi).cl(k),...
                        'timewin', a.timewin,...
                        'uvthresh', a.uvthresh,...
                        'r_splits', a.r_splits, ...
                        'ang_splits', a.ang_splits);
      
      newe.electrode(k).cls{epi} = length(newe.epoch(epi).cl) + (1:length(part_cls));
      
      newe.epoch(epi).cl = [newe.epoch(epi).cl part_cls];
      
    end
  end


function [posi] = posifilt(varargin)
% POSIFILT select position records meeting behavioral criteria (returns logical index)
%   
%  'ep', e.pos struct
%  'posi', initial posi list (e.g. from cl) to filter. If empty, use all
%  'timewin', time window to analyze
%  'filtparamname', name of paramater to filter on (e.g. 'lvel', optional)
%  'filtparamvals', less-than/greater-than values to be selected
%  
%      to select:
%                           
%         >-10 & <10, use [ -10   10]
%         <-10 | >10, use [  10  -10]
%         >5,        use  [   5  Inf]
%         <5,        use  [-Inf    5]
%         ==5             [   5    5]
  

  
  a = struct(...
      'ep',[],...
      'posi',[],...
      'timewin',[-Inf Inf],...
      'filtparamname',[],...
      'filtparamval',[],...
      'filtparamfn', []);


  a = parseArgs(varargin, a);

  p_time = getepd(a.ep, 'time');

  % if no initial posi provided, use all
  if isempty(a.posi)
    posi = true(size(a.ep.dat,1),1);
  else
    if islogical(a.posi),
      posi = a.posi;
    else
      posi = false(size(a.ep.dat,1),1);
      posi(a.posi) = true;
    end
  end
  
    
  % select on timewindow
  posi = posi & ...
         (p_time > a.timewin(1) & ...
          p_time < a.timewin(2));
  
  if ~isempty(a.filtparamname), % do param filtering as well,

    % fn is stored as a string, because fn handles do not play nicely with
    % saving/being stored in userdata, in my experience
    filtfn = eval(a.filtparamfn);
    
    p_dat = getepd(a.ep, a.filtparamname);
    posi = posi & filtfn(p_dat, a.filtparamval);
    
% $$$     
% $$$     filtparamgt = a.filtparamvals(1);
% $$$     filtparamlt = a.filtparamvals(2);
% $$$     
% $$$     p_dat = getepd(a.ep, a.filtparamname);
% $$$     
% $$$     % select on param values
% $$$     if filtparamgt < filtparamlt, % use &
% $$$       posi = posi & ...
% $$$              (p_dat > filtparamgt &...
% $$$               p_dat < filtparamlt);
% $$$     elseif filtparamgt > filtparamlt, % use |
% $$$       posi = posi & ...
% $$$              (p_dat > filtparamgt |...
% $$$               p_dat < filtparamlt);
% $$$     else % gt == lt
% $$$       posi = posi & ...
% $$$              (p_dat == filtparamgt);
% $$$     end  
  end % if paramname provided
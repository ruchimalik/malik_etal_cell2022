function [obj] = obj_cachesearch(obj)
% OBJ_CACHESEARCH function to do cache searches on objects

  otype = obj.type;
  
  % default to no cache hit
  obj.cache_hit = false;

  % initialize cache, if necessary
  if isempty(obj.cache),
    obj.cache = mkcache();
  end
  
  % initialize cache field
  if ~isfield(obj.cache, otype),
    obj.cache.(otype) = [];
  end

  % add input template object (if any) to cache list, for arg-checking
  % purposes
  if ~isempty(obj.template),
    obj.cache.(otype) = [obj.template obj.cache.(otype)];
  end
  
  if isempty(obj.cache.(otype)),
    return
  end
  
  % get correct equality test fn for object type
  switch otype
   case 'parmest',
    obj_eqfn = @(old, new)...
        structcmp(old,new,...
                  {'edesc' 'clnos' 'clperm' 'timewin' 'timebinsize' 'overlapsperbin'...
                   'histbins' 'firstspiketriggers' 'parmestopt' 'tunecurveopt'});
    
   case 'tunecurve',
    obj_eqfn = @(old, new)...
        structcmp(old,new,...
                  {'edesc' 'tunecurveopt'}) &&...
        all(ismember(new.clnos,old.clnos)); % old *contains* requested clnos

   case 'tcpeak',
    obj_eqfn = @(old, new)...
        structcmp(old,new,...
                  {'edesc' 'clnos' ...
                   'tunecurveopt' 'tcpeakopt'});

   case 'raster',
    obj_eqfn = @(old, new)...
        structcmp(old,new,...
                  {'edesc' 'clnos' 'clperm' 'timewin' 'parmwin'...
                   'imwidpix' 'imhtpix' 'whitebg'...
                   'tcpeakopt' 'tunecurveopt' 'rasteropt'});
    
   case 'clhist'
    obj_eqfn = @(old,new)...
        structcmp(old,new,{'edesc' 'clnos' 'binarize' 'rate' 'firstspiketriggers'})...
        &&...
        (~isempty(new.histbins) && ...
         structcmp(old,new, {'histbins'})...
         ||... 
         structcmp(old,new, {'timewin' 'timebinsize'}));
    
   case 'parm',
    % don't use structcmp(tunecurve) b/c we don't care if bins change
    obj_eqfn = @(old,new)...
        structcmp(old,new,{'edesc', 'timewin', 'params'}) &&...
        ((isempty(old.tunecurveopt) && isempty(new.tunecurveopt)) ||...
         structcmp(old.tunecurveopt, new.tunecurveopt, ...
                   {'tuneparams' 'histedges', 'keepbins'}));
    
   case 'cont',
    obj_eqfn = @(old,new)...
        ((~isempty(old.contvar) && structcmp(old,new,'contvar')) ||...
         (isempty(old.contvar) && isempty(new.contvar) && ...
          structcmp(old.contdata,new.contdata,...
                    {'name', 'chanlabels', 'units'})))...
        &&...
        structcmp(old,new,{'chans', 'chanlabels'})...
        &&...
        structcmp(old.contopt,new.contopt,{'autoresample'})...
        &&...
        (isempty(old.contopt.envopt) && isempty(new.contopt.envopt) ||...
         structcmp(old.contopt.envopt, new.contopt.envopt))...
        &&...
        (isempty(old.contopt.filtopt) && isempty(new.contopt.filtopt) ||...
         structcmp(old.contopt.filtopt, new.contopt.filtopt,...
                   {'filttype' 'F' 'atten_db' 'ripp_db' 'datatype' 'sd_t' ...
                   'length_t'}));

   case 'specgram',
    obj_eqfn = @(old,new)...
        (structcmp(old,new,'contvar') ||...
         structcmp(old,new,'contname'))...
        &&...
        structcmp(old,new,{'timewin' 'chans', 'chanlabels', 'contopt' 'specgramopt'});

   case 'filt',
    obj_eqfn = @(old,new)...
        structcmp(old.filtopt, new.filtopt,...
                  {'filttype' 'F' 'sd_t' 'length_t' 'atten_db' 'ripp_db' 'datatype'});    

   case 'seg', 
    obj_eqfn = false; % don't bother testing, for now
    
   otherwise
    error('unrecognized object type in cache, no eqfn');
  end

  %%%%%% test if we can reuse obj from cache
  for oldobj = obj.cache.(otype),
    
    % deal with e/edesc
    if all(isfieldmult(oldobj,{'e' 'edesc'})) && isstruct(oldobj.e),
      oldobj.edesc = oldobj.e.desc;
    end
    
    if obj_eqfn(oldobj, obj),

      % deal with e/edesc
      if all(isfieldmult(obj,{'e' 'edesc'})),
        oldobj.e = obj.e;
      end
      
      obj = oldobj;
      obj.cache_hit = true;
      return;
    end
  end
  

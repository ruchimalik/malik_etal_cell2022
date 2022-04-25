function [e args] = imexpt_multi(varargin)
% IMEXPT - import full experiment and clusters
%
% $Id: imexpt_multi.m 529 2010-02-09 00:18:41Z tjd $
%
% function [e defaults] = imexpt_multi(
%
%                 'session_import_prefs', [],...
%
%                 'datadir', [],
%                 'desc', [],
%                 'spikeparams', {'time' 'id'}.
%
%                 'd1d2distwin', [],
%                 'xyvelth', [],
%                 'xyangvmagth,[],
%                 'veliters',[],
%
%                 'clusterstats', true,
%
%                 'parm_uvthresh', [],
%                 'parm_adoffset', [],
%
%                 'epoch', [],
%                 'epochfile', [],
%                 'trodenames', [],
%                 'trodenameregexp', ''.
%                 'trodegains', {},...
%                 'troderates', {},...
%
%                 'region', [],
%
%                 'eposstruct', [],
%                 'eposfile', [],
%                 'velfsec',0.25,
%
%                 'etrackstruct', [],
%                 'etrackfile', [],
%                 'etrackpts', [],
%                 'etrackptsfile', [],
%
%                 'clfileregexp', '^cl-\d+$'
%                 'clnumregexp', '^cl-(\d+)$'
%
%             flag args:
%                  'noprompt', true)
%
%   NOTES: if a file is not found on the matlab path, we will look for it
%          in the 'datadir' directory.
%
% Calling imexpt with 'getargs' as the first argument, will return only the
% list of default arguments in 'args'
%
%
% ***** Basic arguments 
%
% 'prefs_struct' = struct containing prefs 
%
% 'datadir' = directory containing adextract'd data to import for an
%          experiment. (Often data/<animal>/<day>/extract) default: [] =
%          current directory.
%
% 'desc' = unique name for experimental struct. by default, we construct one
%          from the data directory structure and the selected epochs
%          (e.g. paul_05_run)
%
% 'spikemode' = how to import spikes into cl structs:
%          -'cl-across': match clusters by number across epochs, only
%           include clusters occurring in all requested epochs. cl-name
%           is <trodename>/<clname>. (Default; see 'HOW TO USE', below)
%          -'cl-all': no matching across epochs; import each cl-file as a
%           new cluster. trodename is <trodename>/<epochname>/<clname>
%          -'parm': import all spikes (meeting parm_uvthresh) from each
%           trode.
%
% 'spikeparams' = parameters to load from the cl or parm files. Cell array
%          of strings. Empty means load all params in cl-file. ({'time', 'id'})
%
% 'parm_uvthresh' = threshold (in microvolts) for spikes to be imported
%          from parmfile if spikemode is 'parm'
%
% 'parm_adoffset' = offset for spike amplitude in parm files (usually 0
%          or 2048) depending on A-D card used. Only used when spikemode is
%          'parm'.
%
%
% ***** information about behavior, epochs to import
%
% 'epoch' = an array of epoch structs, usually generated in a session_prefs
%          file, using the 'imepoch' and 'addepoch' functions.  
%          Each epoch struct contains fields:
%            -name (must be unique, and correspond to clsubdirs if
%             applicable)
%            -times (must be non-overlapping)
%            -type : behavior type ('field': just use xy,
%             'splinetrack': do spline-based linearization)
%            -cmperpix: cm per camera pixel in .pos file
%            -headdiroffset: angle in degrees to be added to xydir to
%             correct for LED bar orientation.  Convention is front-back
%             positioning, with diode1 in front.(i.e. When animal is running
%             from start to end of track--a +ve lvel--relative angle between
%             diodes and track should be 0 degrees). [] means choose from [0,
%             90, 180, 270] degrees as appropriate.
%            -im_cl: whether to import clusters from clsubdir named for
%             this epoch
%            -im_parm: whether to import spikes from parm (.pxyabw) file
%             for this epoch
%
% 'epochfile' = .mat file containing an epoch struct named 'epoch'
%
%
% ***** information about electrodes, recording locations
%
% 'trodenameregexp' = the regexp to use to find electrode directories. E.g.
%           '^t\d{2}$', or '^[a-i][12]\d{2}$'.  Use '[]', with the
%           apostrophes, to select no electrodes. 
%
% 'trodenames' = cell array of electrode names/electrode dir names
%
%
% 'region' = an array of region structs, usually generated in a
%          session_prefs file. Each region struct contains fields:
%            -name (e.g. 'ca1', 'thal', must be unique)
%            -electrodes {cell array of trode names}
%            -im_eeg: create an eeg cdat for these tets (boolean)
%
% **** Params for linear interpolation ('posinterp.m') 
%
% 'd1d2distwin' = range of valid interdiode distances. ( [] = accept all
% points, try [3 20] )
%
% 'xyvelth' = velocity threshold, in cm/s, for posinterp: pts discarded
% until overall velocity below this. ( [] = accept all points, try 85)
%
% 'xyangvmagth' = angular velocity threshold, in degrees/s, for
% posinterp. Same treatment as xyvelth. ([] = accept all points, try 1000)
%
% 'veliters' = max # of iterations of posinterp to get rid of high
% vel/angvel points.([] use posinterp default, try 8)
%
% **** Filtering position records
%
% 'velfsec' = smooth linear/xy velocity/speed with a gaussian kernel with
% specified SD (0.25 s)
%
% **** Cluster statistics
% 
% 'clusterstats': compute cluster statistics: l-ratio, isolation
% distance, firing rate. (true)
%
% **** Pre-calculated values for track and pos structs.
%
% Zero or One of:
%
%   'eposstruct': an existing array of e.epoch.pos. structs (1 per epoch)
%   'eposfile': a .mat file containing such a struct
%
% Zero or One of: (for n tracks)
%
%   'etrackstruct' = an array (n x 1) of existing e.track structures
%   'etrackfile' = path to file containing an etrackstruct, in a variable 'etrk'
%   'etrackpts' = spline control points ({[m x 2] x n} coordinates;
%                 (usually only used for debugging).
%
% **** Tweaky stuff
%
% 'clfileregexp' = the regexp to use to find cluster files. 
% ('^cl-\d+$' finds 'cl-1', 'cl-15'.)
%
% 'clnumregexp' = the regexp to use to find the cluster number from a
% cluster file name. If empty, do not reorder clusters.
%('^cl-(\d+)$' finds '15' from 'cl-15'.)
%
%
% **** Flag arguments
%
% 'noprompt' = try to run in batch mode without prompting user (true)
%
%
%
% HOW TO USE:
%
% Designed to make it easy to have different # of cells in some epochs
% (e.g. more in 'run' than can be followed in 'sleep')
%
% Sample workflow in xclust3
%  -Load all pts, create 'pre-sleep', 'run', and 'post-sleep' epochs
%  -Save, e.g., ../paul05.epoch file
%  -for each electrode:
%    -Load 'run' epoch points, cluster e.g. 8 cells
%    -save cl-files in run/ subfolder
%    -load 'post-sleep' epoch, modify/delete/add clusters, **keeping mapping
%     of 8 clusters to same cells**, adding any new cells starting at cl-9,
%     save clfiles to post-sleep/ subfolder.
%    -repeat for other epochs.
%
% Key point: When clustering, use the cl-# as a way of associating clusters
% across epochs. When you import 'run', you will get all the run cells. When
% you import 'run' and 'post-sleep', you will only get the cluster #s that
% are preserved across both epochs.
%
% (The above assumes you are using the default 'cl-across' spikemode. If you
% use 'cl-all', then no attempt is made to match clusters across epochs, and
% each cl-file gets its own 'cl' struct)

% datadir "data directory" structure ([] = optional)
%    -*.p
%    -[etrk.mat] e.track structure
%    -[track.txt] track points
%    -[*.epoch]
%    -[trodedir] ... "electrode directory"
%       -parmfile
%       -ttfile or *.header file for amp gains, etc
%       -epochdirs ... cluster subdirectory (named for an epoch)
%         -[cl-*] ...
%
%
% Tom Davidson (tjd@mit.edu)
% distantly related to Linus Sun's (linus@mit.edu) imexpt.m version 1.4

% TODO:
%  -region

%start timers
tic;
eatclocks = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants / Prefs

% posinterp should use single-diode entries to update interpolation
useonegood = true;

% Automatically updated by SVN, included in metadata below
svn_id_string = ('$Id: imexpt_multi.m 529 2010-02-09 00:18:41Z tjd $');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Arguments / check for file existence

defaults = struct( ...
    'datadir', [],...
    'desc', [],...
    'epoch',[],'epochfile',[], ...
    'spikemode','cl-across',...
    'trodenames', {{}}, 'trodenameregexp', '[]', ...
    'trodegains', {{}},...
    'troderates', {{}},...
    'region', [],...
    'spikeparams', {{'time', 'id'}},...
    ...
    'd1d2distwin', [], ...
    'xyvelth', [], ...
    'xyangvmagth', [], ...
    'veliters', [], ...
    'velfsec', 0.25,...
    'clusterstats', [],...
    'parm_uvthresh', [],...
    'parm_adoffset', [],...
    'eposstruct',[], 'eposfile',[],...
    'etrackstruct',[],'etrackfile',[], 'etrackpts', [], ...
    'clfileregexp', '^cl-\d+$', ...
    'clnumregexp', '^cl-(\d+)$', ...
    'noprompt', 1);

% pass out the default argstruct, if requested, on imexpt('getargs');
if strcmp(varargin{1},'getargs')
  e = defaults; % put it in both args
  args = defaults;
  return;
end

if strcmp(varargin{1},'argstruct'),
  if length(varargin) > 2 || ~isstruct(varargin{2}),
    error(['If ''argstruct'' is provided, a struct must be provided as the ' ...
           '2nd (and only other) argument.']);
  end
  args = varargin{2};
else  
  % parse arguments
  args = parseArgs(varargin, defaults);
end

% datadir required
if isempty(args.datadir),
  error('must provide a ''datadir''.');
else
  datadir = args.datadir;
end

if ~iscell(args.etrackpts) && ~isempty(args.etrackpts)
  args.etrackpts = {args.etrackpts};
end

% use defaults if empty array passed in (override parseArgs)
if isempty(args.trodenameregexp) && isempty(args.trodenames)
  warning(['no ''trodenameregexp'' or ''trodenames'' provided, no electrodes ' ...
           'will be imported']); %#ok
  args.trodenameregexp = '[]';
end

if isempty(args.clfileregexp)
  args.clfileregexp = defaults.clfileregexp;
end

if isempty(args.clnumregexp)
  args.clnumregexp = defaults.clnumregexp;
end

if ~isempty(args.eposfile) && ~ischar(args.eposfile),
  error('Argument ''eposfile'' must be a string.');
end

if ~isempty(args.eposstruct) && ~isstruct(args.eposstruct),
  error('Argument ''eposstruct'' must be a struct.');
end

if ~isempty(args.etrackstruct) && ~isstruct(args.etrackstruct),
  error('Argument ''etrackstruct'' must be a struct.');
end

if ~isempty(args.etrackfile) && ~ischar(args.etrackfile),
  error('Argument ''etrackfile'' must be a string.');
end

% did user give us appropriate # of pos args?
switch (sum([...
    ~isempty(args.eposstruct) ...
    any(args.eposfile)]))
 case 0
  % OK, we will generate from .pos file
 case 1
  % OK, but check for arg conflicts

  if ~isempty(args.d1d2distwin) || ~isempty(args.xyvelth) ||  ...
        ~isempty(args.xyangvmagth) || ~isempty(args.veliters),
    error('Can''t provide posinterp params and a pre-calculated pos struct');
  end

 otherwise
    error('Only one of ''eposstruct'' and  ''eposfile'' can be provided');
end

% Validate eposfile arg
if args.eposfile,
  if ~(exist(args.eposfile, 'file') == 2),
    % not a full path; see if it's a file in the data directory
    if exist([datadir '/' args.eposfile], 'file') == 2,
      disp(['File found in ''datadir'':' args.etrackfile]);      
      args.eposfile = [datadir '/' args.eposfile];
    else
      error (['''eposfile'': ' args.eposfile ' missing.']);
    end
  end
end


% did user give us appropriate # of track args?
switch sum([~isempty(args.etrackstruct) ...
            ~isempty(args.etrackfile) ...
            ~isempty(args.etrackpts)]),
 case 0
  if args.noprompt,
    warning(['One of ''etrackstruct'', ''etrackfile'', ''etrackpts'', ' ...
             'and ''etrackptsfile'' must be provided to run in ''noprompt'' ' ...
             'mode.']); %#ok
    args.noprompt = false;
  end
 case 1
  %ok
 otherwise
  error('Only one of ''etrackstruct'', ''etrackfile'', and ''etrackpts'' can be provided');
end

% Validate etrackfile arg
if args.etrackfile,
  if ~(exist(args.etrackfile, 'file') == 2),
    % not a full path; see if it's a file in the data directory
    if exist([datadir '/' args.etrackfile], 'file') == 2,
      disp(['File found in ''datadir'':' args.etrackfile]);      
      args.etrackfile = [datadir '/' args.etrackfile];
    else
      error (['''etrackfile'': ' args.etrackfile ' missing.']);
    end
  end
end



%%%%%%%%%%%%%%%%%%%%%
% add some metadata %
%%%%%%%%%%%%%%%%%%%%%
[status hostname] = unix('hostname');
if status == 0, % success
  e.imexpt.hostname = hostname;
end
e.imexpt.svn_id = svn_id_string;
e.imexpt.path = datadir;
e.imexpt.date = datestr(now);

% save the args used to make this, but not any large structs.
e.imexpt.args = args;
if ~isempty(e.imexpt.args.eposstruct), 
  e.imexpt.args.eposstruct = 'struct provided on command line';
end
if ~isempty(e.imexpt.args.etrackstruct), 
  e.imexpt.args.etrackstruct = 'struct provided on command line';
end


%%%%%%%%%%%%%%%%%%%%%
% Load epoch times  %
%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n'));
disp('======================================');
disp('= Loading Epoch Data              =');
disp('======================================');
disp(sprintf('\n'));

if ~isempty(args.epoch)
  disp('Using epoch struct from args');
else
  disp(['Loading Epoch file: ' args.epochfile]);
  % load epoch
  load(args.epochfile, 'epoch', '-MAT');
  if ~exist('epoch', 'var'),
    error(['Error loading variable ''epoch'' from file: ' args.epochfile]);
  else
    args.epoch = epoch;
    clear('epoch');
  end
end

neps = numel(args.epoch);
e.epoch = repmat(struct(), 1, neps); % with no fields yet

[e.epoch.name] = deal(args.epoch.name);
[e.epoch.times] = deal(args.epoch.times);
[e.epoch.behav] = deal(args.epoch.behav); % 'run', 'sleep'
[e.epoch.envt] = deal(args.epoch.envt); % 'field', 'splinetrack'
[e.epoch.tracki] = deal([]);
[e.epoch.pos] = deal([]);

if any(isinf([e.epoch.times])) || any([e.epoch.times]<0)
  error('Infinite/negative epoch times not supported');
end

% determine # of tracks we need to build

ntracks = 0;
for k = 1:neps
  if strcmp(args.epoch(k).envt, {'splinetrack'})
    ntracks = ntracks+1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load/process position info  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n'));
disp('======================================');
disp('= Loading Position Data              =');
disp('======================================');
disp(sprintf('\n'));

%%% make sure we have a .pos file to import
posfilename = dir([datadir '/*.pos']);
switch size(posfilename,1)

  case 1
    posfilename = [datadir '/' posfilename.name];
    
  case 0
    error('No *.pos file in data directory.');
 
  otherwise,
    error('More than one *.pos file in data directory.');

end

if ~isempty(args.eposfile) || ~isempty(args.eposstruct)
  %% use eposfile/struct if provided
  if ~isempty(args.eposfile),
    finfo = dir(args.eposfile);
    disp(['Using file: ' args.eposfile ' for pos structure']);
    disp(['Last saved: ' finfo.date]);
    
    % load epos
    load(args.eposfile, 'epos', '-MAT');
    if ~exist('epos', 'var'),
      error(['Error loading variable ''pos'' from file: ' args.eposfile]);
    end
  else
    % use epos struct, if provided
    epos = args.eposstruct;
  end

  % check that we have pos structs for each epoch
  if numel(epos)~=neps
    error('Provided eposstruct/eposfile does not match number of epochs');
  end
  
  % distribute 
  for k = 1:neps
    e.epoch(k).pos = epos(k);
  end
  clear('epos');

  
else
  %% import pos data from scratch
  for epochi = 1:neps,
    
    disp(sprintf('\n'));
    disp(['Importing position data for epoch ''' args.epoch(epochi).name '''']);
    disp(['Epoch times: ' num2str(args.epoch(epochi).times)]);
    
    if isempty(args.epoch(epochi).cmperpix),
      cmppstring = '[] (No conversion)';
    else
      cmppstring = num2str(args.epoch(epochi).cmperpix);
    end

    disp(['Epoch cmperpix: ' cmppstring]);
    disp(sprintf('\n'));

    % inefficient in that we read the whole pos-file each time, but not a
    % bottleneck, scientifically !!!!! Stop working on this and do science!
    % Yes, you!

    ep = impos([], posfilename, args.epoch(epochi).times, args.epoch(epochi).cmperpix);

    % Interpolate missing diodes, calculate head dir, xy-velocity, etc
    %   *note we don't add linearized data until after we've made/loaded the
    %    track. See call to posaddlin, below.
    
    disp(sprintf('\n'));
    disp('Processing / interpolating position data (posinterp.m)');
    
    ep = posaddxy (ep);
    
    ep = posinterp(ep,...
                   'useonegood',useonegood,...
                   'd1d2distwin', args.d1d2distwin,...
                   'xyvelth', args.xyvelth, ...
                   'xyangvmagth', args.xyangvmagth,...
                   'veliters', args.veliters, ...
                   'verbose', true);

    % correct head direction
    ep = posaddheaddir(ep, args.epoch(epochi).headdiroffset);

    % copy epoch fields to pos struct
    ep.epochname = args.epoch(epochi).name;
    ep.epochtimes = args.epoch(epochi).times;
    ep = reorderstructure(ep, 'epochname', 'epochtimes');

    e.epoch(epochi).pos = ep;
    
    % save epos in a temp file so we don't have to recalculate if errors
    % encountered later
    if ~exist('temp_epos', 'var'),
      temp_epos = [tempname '.mat'];
    end
    disp(['Saving temporary e.epoch.pos structs to ' temp_epos]);
    epos = [e.epoch.pos]; %#ok
    save(temp_epos, 'epos');
    clear epos;
    
  end
  
end


% plot position points for visual inspection
mapfig = figure;
rc = ceil(sqrt(neps));
for epochi = 1:neps,
  mapax(epochi) = subplot(rc,rc,epochi, 'parent', mapfig); %#ok
  trackmap('ax', mapax(epochi),...
           'epos',e.epoch(epochi).pos,...
           'plots',{'centxy' 'interp'});
  title(mapax(epochi), e.epoch(epochi).name, 'interpreter', 'none');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If epoch is of type 'splinetrack', call MKETRK to create a track model and
% spline, evaluate the spline, and add all to e.track. Uses provided
% track or spline control points if requested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n'));
disp('=======================');
disp('= track linearization =');
disp('=======================');

if args.etrackfile,
  finfo = dir(args.etrackfile);
  disp(['Using file: ' args.etrackfile ' for track structure']);
  disp(['Last saved: ' finfo.date]);
  
  % load etrk
  load(args.etrackfile, 'etrk', '-MAT');
  if ~exist('etrk', 'var'),
    error(['Error loading variable ''etrk'' from file: ' args.etrackfile]);
  else
    args.etrackstruct = etrk; %#ok
    clear('etrk');
  end
end

% use track if provided
if ~isempty(args.etrackstruct),
  
  if numel(args.etrackstruct)~=ntracks
    error('provided etrackstruct contains wrong number of tracks');
  end
  e.track = args.etrackstruct;
  
else
  tracki = 0;
  for k = 1:neps
    if ~strcmp(args.epoch(k).envt, 'splinetrack')
      e.epoch(k).tracki = [];
    else
      
      tracki = tracki+1;
      disp(['Creating track #: ' num2str(tracki)]);
      
      % UGLY HACK
      % Figure restacking is completely borked in R2007a under ubuntu: calling
      % 'figure(h)' does not restack, despite what documentation says. Even
      % modal windows stay buried. Instead, copy to a new figure:
      
      trackfig(tracki) = figure('color', 'w');  %#ok
      copyax = copyobj(mapax(k), trackfig(tracki));
      set(copyax, 'position', [0 0 1 0.9]);
      
      figpos = get(trackfig(tracki), 'position');
      set(trackfig(tracki), 'position', [0 0 figpos(3:4)]);

      if ~isempty(args.etrackpts)
        newtrk =  mketrk (copyax, args.etrackpts{k}, args.noprompt);
      else
        newtrk = mketrk (copyax);
      end

      % update pointers from epoch <-> track
      e.epoch(k).tracki = tracki;
      newtrk.epochi = k;
      
      e.track(tracki) = newtrk;
      clear newtrk;
    end
  end
end
% so we don't have to recalculate if errors encountered later
if isempty(args.etrackfile),
  temp_etrk = [tempname '.mat'];
  disp(['Saving e.track to ' temp_etrk]);
  etrk = e.track; %#ok
  save(temp_etrk, 'etrk');
  clear etrk;
end

% distribute tracki indexes to epoch, pos (may not be present in saved file)
[e.epoch.tracki] = deal([]);
for k = 1:ntracks
  e.epoch(e.track(k).epochi).tracki = k;
end


for k = 1:neps
  if ~isempty(e.epoch(k).tracki)
    try %could fail if axes are deleted
      trackmap('ax',mapax(k),...
               'etrack',e.track(e.epoch(k).tracki),...
               'plots',{'track', 'voronoi'});
      drawnow;
    catch
    end
  end
end

% same x/y scale for all envts in mapax plot
if numel(mapax)>1
  xls = cell2mat(get(mapax, 'xlim'));
  yls = cell2mat(get(mapax, 'ylim'));
else
  xls = get(mapax, 'xlim');
  yls = get(mapax, 'ylim');
end
set(mapax, 'xlim', [min(xls(:,1)) max(xls(:,2))]);
set(mapax, 'ylim', [min(yls(:,1)) max(yls(:,2))]);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add linearized position data to e.epoch.pos %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(sprintf('\n'));
disp('Adding linearized track data to e.epoch.pos')

if ~isempty(args.eposstruct),
  disp('(re-linearizing pos data in user-provided pos struct)');
end

for k = 1:neps
  % only linearize when there is an associated track
  if ~isempty(e.epoch(k).tracki)
    e.epoch(k).pos = posaddlin (e.epoch(k).pos, e.track(e.epoch(k).tracki));
  
    e.epoch(k).pos = posfiltsec (e.epoch(k).pos, 'lvel', args.velfsec);
    e.epoch(k).pos = posfiltsec (e.epoch(k).pos, 'lspeed', args.velfsec);
    
  end
  e.epoch(k).pos = posfiltsec (e.epoch(k).pos, 'xyvel', args.velfsec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use description provided by user, or build default description from
% datadir /data500/pierre/12 -> pierre_12 and prompt user:

if ~isempty(args.desc),
  e.desc = args.desc;
  
else
  descdef = [];
  [tok] = regexp (datadir, '/+(\w+)/+(\w+)/*$', 'tokens'); % get last two dir names
  if size(tok{:},2) == 2,
    descdef = [tok{1}{1} '_' tok{1}{2}];
  end

  desc = [];
  if ~args.noprompt,
    desc = input(['A unique name for this experiment''s dataset? [' descdef ']: '], 's');
  end

  if isempty(desc),
    e.desc = descdef;
  else
    e.desc = desc;
  end
end


%FIXME regions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next load all clusters   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(sprintf('\n'));
disp('====================================');
disp('= Loading Spike/Electrode Data from Expt =');
disp('====================================');


clustidx = zeros(1,neps);
trodei = 0;

% init electrode field (in case we load no clusters)
e.electrode(1) = struct('name',[],...
                        'cls', [],...
                        'region', [],...
                        'type', [],...
                        'ad', [], ...
                        'rate',[], ...
                        'gain', []);

%% get electrode directores
% (use either trodenames, or trodenameregexp)

if ~isempty(args.trodenames)
  for tn = args.trodenames
    if ~isdir([datadir '/' tn{:}])
      error(['Requested electrode directory doesn''t exist: ' tn{:}])
    end
  end
  trodedirs = args.trodenames;

else
  % use trodenameregexp
  % Find all directories in datadir
  dirlist = dir(datadir);
  dirlist = dirlist([dirlist.isdir]);
  trodedirs = {};
  for dirname = {dirlist.name};
    if regexp(dirname{:}, args.trodenameregexp);
      trodedirs(end+1) = dirname; %#ok
    end
  end
end

%iterate over all electrode directories
for tdir = trodedirs
  
  tdir = tdir{:}; %#ok
  tname = tdir; 
  tdir = [datadir '/' tdir]; %#ok
  
  trodei = trodei + 1;
  
  disp(['Processing electrode directory: ' tdir]);
  
  % Look for .tt file or headerfile
  ttfile = dir([tdir '/*.tt']);
  headerfile = dir([tdir '/*.header']);
  switch size(ttfile,1)
   case 1
    if isempty(headerfile),
      unix (['header ' tdir '/' ttfile.name ' > ' tdir '/' ttfile.name '.header']);
      headerfile = dir([tdir '/*.header']);
    end
   case 0,
    if isempty(headerfile),
      % no need to warn--we either have gains or we don't (error later)
      %      warning (['No ttfile or header file found in directory ' tdir]); %#ok
    end
   otherwise,
    error (['More than one ttfile found in directory ' tdir]);
  end
  headerfile = [tdir '/' headerfile.name];
  
  % set up e.tetrode record
  e.electrode(trodei).name = tname;
  e.electrode(trodei).cls = {};
  e.electrode(trodei).region = [];
  e.electrode(trodei).regioni = [];
  e.electrode(trodei).type = 'tetrode';
  e.electrode(trodei).ad   = import_adinfo(headerfile);
  
  
  adhdr = e.electrode(trodei).ad;
  if isempty(adhdr) || (isstruct(adhdr) && isempty(fieldnames(adhdr))),

    % no header file, get rate/gain either from args, or, as a last
    % resort, manually
    
    % can only provide trodegains if there is a matching set of trodenames
    if ~isempty(args.trodegains) && ~isempty(args.trodenames)
      e.electrode(trodei).gain = args.trodegains{trodei};
    else
      e.electrode(trodei).gain = input(['Gain for electrode ' e.electrode(trodei).name '?:']);
    end

    if ~isempty(args.troderates) && ~isempty(args.trodenames)
      e.electrode(trodei).rate = args.troderates(trodei);
    else        
      defaultrate = 250000;
      e.electrode(trodei).rate = input(['Rate for this electrode? (' ...
                          num2str(defaultrate) '): ']);
    end
    
  else

    e.electrode(trodei).rate = adhdr.rate;
    
    % get electrode gain
    fldmatch = regexp(fieldnames(adhdr)', 'chan\d+ampgain','match','once');
    trodegain = [];
    for fld = fldmatch,
      if ~isempty(fld{1}),
        trodegain = [trodegain adhdr.(fld{1})]; %#ok
      end
    end
    e.electrode(trodei).gain = trodegain;
  end

  switch args.spikemode
   case {'cl-across', 'cl-all'}

    % build requested cldir list from epochs.im_cl:
    cldirs = {args.epoch([args.epoch.im_cl]).name};
    ncldirs = size(cldirs,2);

    % get index into e.epoch for each cldir
    cldir_epi = find([args.epoch.im_cl]);
    
    % Get list of cl names for each epoch directory
    for k = 1:ncldirs
      fnames = dir([tdir '/' cldirs{k}]);
      fnames = {fnames.name};
      
      % select only cl-files (using regexp)
      fmatch = regexp(fnames, args.clfileregexp, 'match');
      clnames{k} = [fmatch{:}]; %#ok
      clnames{k} = substrnumsort(clnames{k},args.clnumregexp); %#ok
    end
    
    % if spikemode is 'cl-across', only import clusters that appear in all
    % requested epoch dirs 
    if strcmp(args.spikemode, 'cl-across')
      if ncldirs > 1
        clacross = clnames{1};
        for k = 2:ncldirs
          clacross = intersect(clacross, clnames{k});
        end
        clacross = substrnumsort(clacross,args.clnumregexp);
        [clnames{1:ncldirs}] = deal(clacross); %#ok
      end
    end
    
    for k = 1:ncldirs
      
      cldir = cldirs{k};
      epi = cldir_epi(k);
      
      e.electrode(trodei).cls{epi} = {};

      % Load cluster data
      for clname = clnames{k},
        clname = clname{:}; %#ok
        clfile = [tdir '/' cldir '/' clname];
        disp(['Importing Cluster File: ' clfile]);
        clustidx(epi) = clustidx(epi) + 1;

        e.epoch(epi).cl(clustidx(epi)) = imclust(clfile, args.spikeparams, e.epoch(epi).pos);
        e.epoch(epi).cl(clustidx(epi)).electrode = trodei;
        e.epoch(epi).cl(clustidx(epi)).name = [cldir '/' clfile];
        
        % build list of cls on this electrode
        e.electrode(trodei).cls{epi} = [e.electrode(trodei).cls{epi} clustidx(epi)];
        
      end
    end
    
    % for epochs without im_cl, make e.cl empty
    for k = find(~[args.epoch.im_cl])
      e.epoch(k).cl = [];
    end
    
   case 'parm'

    % build requested cldir list from epochs.im_parm:
    parm_epi = find([args.epoch.im_parm]);
    
    parmfile = dir([tdir '/*.pxyabw']);
    switch numel(parmfile)
     case 0
      error(['No parmfile found for requested electrode--can''t import ' ...
             'spikeparams: ' tdir]);
     case 1
      parmfile = [tdir '/' parmfile.name];
     otherwise
      error(['Multiple parmfiles found for requested electrode: ' tdir]);
    end

    % import spikes as a one cluster per tetrode, per epoch
    for epi = parm_epi
      
      disp(['processing parmfile: ' parmfile]);
      
      e.epoch(epi).cl(trodei) = ...
          parm2cl(...
              'epos', e.epoch(epi).pos,...
              'parmfile', parmfile,...
              'tetname', tname,...
              'clno', [],...
              'trodeno', trodei,...
              'gain', e.electrode(trodei).gain,...
              'rate', e.electrode(trodei).rate,...
              'uvthresh', args.parm_uvthresh,...
              'timewin', e.epoch(epi).times,...
              'offset', args.parm_adoffset);
      
      e.electrode(trodei).cls{epi} = trodei;
      
    end

    % for epochs without im_parm, make e.cl empty
    for k = find(~[args.epoch.im_parm])
      e.epoch(k).cl = [];
    end
    
   otherwise
    error('unsupported spikemode');
    
  end
end

disp('Done loading all clusters...');
disp(sprintf('\n'));

%% parse regions
disp('Loading region information...');

nregions = numel(args.region);

if nregions==0
  e.region = [];
else
  rnames = {args.region.name};
  if numel(unique(rnames))<numel(rnames)
    error('repeated region names')
  end

  for k = 1:nregions
    rname = args.region(k).name;
    e.region(k).name = rname;

    trodenames = args.region(k).trodenames;
    if numel(unique(trodenames))<numel(trodenames)
      error(['repeated trodename in ''region'': ' rname]);
    end
    e.region(k).trodenames = trodenames;

    [dum trodei] = ismember(trodenames, {e.electrode.name});
    if any(trodei == 0)
      error(['no match for trodename in region: ' rname]);
    end
    e.region(k).trodei = trodei;
    
    % add region/i to each electrode
    [e.electrode(trodei).region] = deal(rname);
    [e.electrode(trodei).regioni] = deal(k);
  end
  
  % now that we have trodei, we can look up names from e.electrode
  e.region = rmfield(e.region, 'trodenames');
  
end
disp('Done loading region information...');
disp(sprintf('\n'));


%%% add cluster stats (l-ratio, isolation distance, firing rate...)
if ~isempty(args.clusterstats) && args.clusterstats && ~strcmp(args.spikemode, 'parm')
  disp(sprintf('\n'));
  for k = find([args.epoch.im_cl]),
    disp(['Computing cluster statistics for epoch: ' e.epoch(k).name '...']);
    e.epoch(k).cl = addclusterstats(...
        'e', e,...
        'epochi', k,...
        'datadir', datadir);
  end;
  disp('Done computing all cluster stats.');
  disp(sprintf('\n'));
end


%% clean up and go home

if isempty(args.etrackfile),
  disp(['Deleting temporary etrk file: ' temp_etrk]);
  delete(temp_etrk);
end

if isempty(args.eposfile),
  disp(['Deleting temporary epos file: ' temp_epos]);
  delete(temp_epos);
end


disp([sprintf('\n') 'Time Elapsed to Import Data: ' num2str(toc) ]);
disp([ 'CPU Time Elapsed to Import Data: ' num2str(cputime-eatclocks) ]);
disp(sprintf('\n'));
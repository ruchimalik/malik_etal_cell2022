%% Sample code to import, demodulate, and normalize raw photometry signals recorded with a 2-color, 2-frequency rig 

% written by Tom Davidson (tjd@alum.mit.edu) April 2016 
% modified by Clio Korn (clio.korn@ucsf.edu) spring 2018
% modified by Ruchi Malik (ruchi.malik@ucsf.edu) April 2021

% NOTE: in order to run this script, Tom Davidson's collection of Matlab functions (in folders called 'tjd-shared-code' and 'tjd-shared-code-extra) must be in your Matlab path


addpath(genpath('/Users/Ruchi/photometry scripts'));


%% STEP 1: set up for analysis 
% see Tom D's function documentation for information on what each parameter means

% define analysis parameters

tanksdir = ['/Users/Ruchi/photometry scripts/Data/'];          % 'tanksdir' is the name of the path to the folder containing the photometry data tank and Matlab scripts (including this one)
if ~exist(tanksdir, 'dir'),
    error('DataTanks folder not found at: %s', [path_to_Example_Data tanksdir])
end 
tankname = 'FP-RZ5P-170421-164150';                             % 'tankname' is the name of the photometry data tank (which contains data from all sessions for all animals in your experiment)
blockname = 'RM_M6_rD2-090101-053914';                          % 'blockname' is the name of the folder containing data from a particular session in a particular animal

savename = 'RM_M6_rD2.mat';                                                              
 
Raw1_chanlabels = {'Ref1X', 'Ref2X', 'Det1', 'blank'};                      % define the different channels/components of the raw data
signal_labels = {'470nm', '405nm'};                                         % define the names of the 470 and 405 components of the signal


% define time window for full session
t_1 = 0;                                                     % t_1 is the FIRST timestamp for behavioral event that BEGINS each trial
t_2 = 1320;                                       % t_2 is the LAST timestamp for the behavioral event that ENDS each trial            
FP_timewin = [(t_1 - 5) (t_2 + 5)];                                       % define the session time window to start 30sec before t_1 and end 30sec after t_2 (adjust buffer time as appropriate)
                                                                             

% define demodulation parameters (see contdemodulate.m for documentation)
cfg.demod_BW_F = [1 3];                                                     % bandwidth (Hz) 
cfg.demod_ripp_db = 0.1;                                                    % filter design parameter: bandpass ripple
cfg.demod_atten_db = 50;                                                    % filter design parameter: stopband rejection

% define normalization parameters (see FP_normalize.m for documenation)
cfg.FPnorm_norm_type = 'fit';                                               % type of post-fit normalization to perform 'fit' or 'const'
cfg.FPnorm_control_LP_F = [1 3];                                            % low-pass filter transition band (Hz)
cfg.FPnorm_dFF_zero_prctile = 1;                                            % arbitrarily re-zero data to prctile

% define baseline rig fluorescence for each channel with animal not plugged in (units = Volts) - use [] if you didn't measure this 
cfg.rig_baseline_V = [];                                                    % if you did measure this, code for rig baseline 480/405 looks like this:  cfg.rig_baseline_V = [130 62];



%% STEP 2: load in raw FP signals (detector output, carrier frequencies)
% see Tom D's function documentation for further information

if ~exist('cache', 'var'), cache = mkcache(); end;
fprintf('Loading store: ''Fi1r''\n'); 
S = TDT_Import(tanksdir, tankname, blockname, 'Fi1r'); 
c_FP_Raw = imcont('tdtwave', S, ...
    'chans', [],...
    'chanlabels', Raw1_chanlabels,...
    'dataunits', 'V');
fprintf('-->Done!\n');



%% STEP 3: demodulate raw detector signal
% see Tom D's function documentation for further information

fprintf('Demodulating raw photometry signal ...\n');
[c_Mag, FP_Ref_F, FP_PSDs, cache] = ...
    contdemodulate(c_FP_Raw, ...
    'nsignals', 2,...
    'signal_labels', signal_labels,...
    'bandwidth_F', cfg.demod_BW_F,...
    'ripp_db', cfg.demod_ripp_db,...
    'atten_db', cfg.demod_atten_db,...
    'cache', cache);
fprintf('-->Done!\n');



%% STEP 4: normalize demodulated signal
% see Tom D's function documentation for further information

t_1 = 60;                                                     % t_1 is the FIRST timestamp for behavioral event that BEGINS each trial
t_2 =  c_FP_Raw.tend-20;                                       % t_2 is the LAST timestamp for the behavioral event that ENDS each trial            
FP_timewin = [(t_1) (t_2)];
fprintf('Normalizing photometry signal ...\n');
[c_dFF, c_Regress, bls] = ...
    FP_normalize(c_Mag, ...
    'timewin', FP_timewin,...
    'norm_type', cfg.FPnorm_norm_type ,...
    'control_LP_F', cfg.FPnorm_control_LP_F, ...
    'rig_baseline_V', cfg.rig_baseline_V,...
    'dFF_zero_prctile', cfg.FPnorm_dFF_zero_prctile);
fprintf('-->Done!\n');



%% STEP 5: save outputs for use in Plotting_photoM.m and Quantification_photoM.m scripts

save (savename, 'c_FP_Raw', 'c_dFF', 'c_Mag', 'c_Regress', 'bls', 't_1', 't_2')
% save('Rand_dFF','c_dFF');                                               % this is the key output: the data structure containing the dF/F values
% save('Rand_c_Mag','c_Mag');                                               % this is required for generating some of the plots in Plotting_photoM.m
% save('Rand_c_Regress','c_Regress');                                       % this is required for generating some of the plots in Plotting_photoM.m
% save('Rand_bls','bls');                                                   % this is required for generating some of the plots in Plotting_photoM.m





%%  Quantify photometry data recorded during NOE behavior

% written by Ruchi Malik (ruchi.malik@ucsf.edu) April 2021

% run Processing_photoM.m to import raw data into Matlab, demodulate it, and normalize it (i.e. calculate dF/F)

%% STEP 1: load c_dFF data structure and extract dF/F from it
load('RM_M2.mat');     % this was saved at end of Processing_photoM.m
data = c_dFF.data;     % extract dF/F values from c_dFF and put into vector called 'data'
load('M2_results_noNOE1.mat');
clearvars -except Times c_dFF data

choice =3; %choice is 1 for all NOE, choice is 2 for early NOE, choice is 3 for late NOE

seg = 10;  % segment length is the time (in sec) that will be cut around (i.e. before and after) the event of interest


if choice ==1
    savname = 'NOE_quant_results_all.mat'; % all NOE bouts
elseif choice ==2
    savname = 'NOE_quant_results_early.mat'; % only early NOE bouts
elseif choice ==3
    savname = 'NOE_quant_results_late.mat'; % only late NOE bouts
    
end
%% STEP 2: import timestamp spreadsheet and match to photometry data
% read in timestamps of behavioral events - each column in the spreadsheet is a type of behavioral event, each row is a single trial
times_raw_NOE(:,1) = Times.NOE_start;           % xlsread('NAME OF SPREADSHEET', 'NAME OF TAB', 'CELLS TO BE IMPORTED')
times_raw_NOE(:,2) = Times.NOE_end; 
times_raw_noNOE(:,1) = Times.Obj_added;
times_raw_noNOE(:,2) = Times.Obj_removed;
times_raw_base(:,1) = Times.baseline_start;
times_raw_base(:,2) = Times.baseline_end;


%% STEP 3: extract time windows for analysis
adjustment = c_dFF.tstart;                                           % from first timestamp, subtract the buffer time (in sec) used to define FP_timewin in Step 2 of Processing_photoM.m
times_adjs_noNOE = times_raw_noNOE - adjustment;                                    % subtract the result from the entire matrix of raw timestamps
times_noNOE = round(c_dFF.samplerate * times_adjs_noNOE);                           % multiple the resulting matrix by the photometry signal sample rate; round because the sample rate is not a whole number

times_adjs_NOE = times_raw_NOE - adjustment;                                    % subtract the result from the entire matrix of raw timestamps
times_NOE = round(c_dFF.samplerate * times_adjs_NOE);   

times_adjs_base = times_raw_base - adjustment;                                    % subtract the result from the entire matrix of raw timestamps
times_base = round(c_dFF.samplerate * times_adjs_base);   


% create separate vectors for behavioral events of interest
%trial_num = times(:,1); 
start_interaction = times_NOE(:,1);
end_interaction = times_NOE(:,2);  

start_base = times_base(:,1);
end_base = times_base(:,2);

start_noNOE = times_noNOE(:,1);
end_noNOE = times_noNOE(:,2);

sz = size(c_dFF.data(:,1),1);
t = linspace(1,sz, sz); 

%% STEP 4: convert time windows to frames

% first for no NOE time windows
noNOE_timefrms = [];
% for i = 1: size(start_noNOE,1)
    i = 1;
 j = size(start_noNOE,1);
     if choice == 1
    indx = find(t>start_noNOE(i) & t<end_noNOE(j) | round(t)==start_noNOE(i) | round(t)==end_noNOE(j));
    elseif choice ==2
        indx = find(t>start_noNOE(i) & t<start_interaction(i) + (180*60) | round(t)==start_noNOE(i));
    elseif  choice ==3
        indx = find(t>(start_interaction(i) + (180*60)) & t<end_noNOE(j) | round(t)==(start_interaction(i) + (180*60)) | round(t)==end_noNOE(j));
    end
       
     noNOE_frames = indx;



early_idx = find(start_interaction <= start_interaction(1)+180*60);
earlystart_times = start_interaction(early_idx);
earlyend_times = end_interaction(early_idx);

late_idx = find(start_interaction > start_interaction(1)+180*60);
latestart_times = start_interaction(late_idx);
lateend_times = end_interaction(late_idx);

if choice == 2
    start_interaction = earlystart_times;
    end_interaction = earlyend_times;
elseif choice ==3
    start_interaction = latestart_times;
    end_interaction = lateend_times;
end

% second for NOE time windows
NOE_timefrms = [];

for i = 1: size(start_interaction,1)
    
    indx = find(t>start_interaction(i) & t<end_interaction(i) | round(t)==start_interaction(i) | round(t)==end_interaction(i));
    
        NOE_timefrms{i} = indx;
     
end 

NOE_frames = [];
NOE_frames = cell2mat(NOE_timefrms);



allnoe_frames = horzcat(NOE_timefrms{:});
[a, a1, a2] = intersect(allnoe_frames, noNOE_frames);

noNOE_frames(a2) = [];

dFF_noNOE = c_dFF.data(noNOE_frames);

mean_noNOE = mean(dFF_noNOE);
std_noNOE = std(dFF_noNOE);



base_timefrms = [];
 i = 1;
 j = size(start_base,1);
    
    base_frames = find(t>start_base(i) & t<end_base(j) | round(t)==start_base(i) | round(t)==end_base(j));
     
     
dFF_base = c_dFF.data(base_frames);

mean_base = mean(dFF_base);
std_base = std(dFF_base);


%% Step 5: define the time points of interest  and cut data vector (dF/F for whole session) around behavioral time points of interest, average, and send to Excel (and subsequently Prism) for further perusal and plotting


time_point = start_interaction;

r = randi([base_frames(1)+150 base_frames(size(base_frames,2))-150],1,20);
time_point2 = r';  % random time points during baseline period
 

% re-run below code for each behavioral time point of interest by commenting in/out definitions of 'time_point' and 'totab' in Step 6 above


% create timeline against which to plot cross-trial average signal at behavioral events of interest; event will be centered at time 0 on the x axis

timeline_neg = -seg:(1 / c_dFF.samplerate):0;                               % the negative part of the timeline goes from the negative of the segment value to 0 in increments the size of the sample rate
timeline_pos = fliplr(timeline_neg * -1);                                   % generate the positive part of the timeline by multiplying the negative part by -1 and flipping it
timeline_center = 0;                                                        % the center of the timeline is 0
Timeline = horzcat(timeline_neg, timeline_center, timeline_pos)';           % concatenate these three vectors to get the full timeline


% cut dF/F around behavioral time point of interest : NOE 

Segs = zeros((round(c_dFF.samplerate * seg) + round(c_dFF.samplerate * seg) + 1),length(time_point));    % matrix of zeros, # of rows = length of timeline, # of columns = # of trials
for i = 1:length(time_point);                                               % going sequentially through the timestamp matrix (i.e. trial by trial)
    indx = time_point(i,1);                                                % indx is value of timestamp for trial i
    lower = indx - round(c_dFF.samplerate * seg);                           % define lower bound of cut i.e. indx minus the segment value - multiplied by the sample rate to be in photometry signal units
    upper = indx + round(c_dFF.samplerate * seg);                           % define upper bound of cut i.e. indx plus the segment value - multiplied by the sample rate to be in photometry signal units
    cut = data(lower:upper,1);                                              % select dF/F from the lower to the upper bound
    Segs(:,i) = cut;                                                        % put cut piece of data in Segs matrix
end

 zscore_Segs = (Segs-mean_noNOE)./std_noNOE; % eas using noNOE before not base
 
 
% calculate mean and SEM signal at behavioral time point of interest

Mean = nanmean(Segs,2);                                                     % average across trials (columns in Segs matrix) to get mean signal
SEM = nansem(Segs,2);                                                       % calculate SEM across trials

Mean_zs = nanmean(zscore_Segs,2);

zscore_Segs = movmean(zscore_Segs,30);




% cut dF/F around behavioral time point of interest : baseline

Segs_base = zeros((round(c_dFF.samplerate * seg) + round(c_dFF.samplerate * seg) + 1),length(time_point2));    % matrix of zeros, # of rows = length of timeline, # of columns = # of trials
for i = 1:length(time_point2);                                               % going sequentially through the timestamp matrix (i.e. trial by trial)
    indx = time_point2(i,1);                                                % indx is value of timestamp for trial i
    lower = indx - round(c_dFF.samplerate * seg);                           % define lower bound of cut i.e. indx minus the segment value - multiplied by the sample rate to be in photometry signal units
    upper = indx + round(c_dFF.samplerate * seg);                           % define upper bound of cut i.e. indx plus the segment value - multiplied by the sample rate to be in photometry signal units
    cut = data(lower:upper,1);                                              % select dF/F from the lower to the upper bound
    Segs_base(:,i) = cut;                                                        % put cut piece of data in Segs matrix
end

 zscore_Segs2 = (Segs_base - mean_base)./std_base;
% calculate mean and SEM signal at behavioral time point of interest

Mean_base = nanmean(Segs_base,2);                                                     % average across trials (columns in Segs matrix) to get mean signal
SEM_base = nansem(Segs_base,2);                                                       % calculate SEM across trials

Mean_zs_base = nanmean(zscore_Segs2,2);

zscore_Segs2 = movmean(zscore_Segs2,30);


numean_base = movmean(Mean_zs_base, 30);
figure
plot_Mean_flexBounds(zscore_Segs','r', 'SEM', Timeline');
hold on
plot_Mean_flexBounds(zscore_Segs2','k', 'SEM', Timeline');

%% STEP 6: quantify INDIVIDUAL TRIAL dF/F signals at behavioral time points of interest, then send to Excel
% re-run below code for each behavioral time point of interest by commenting in/out definitions of 'time_point' and 'totab' in Step 6 above


% define time window (relative to behavioral event) within which quantification is to be performed (change as appropriate)

window = seg/2;                                                                 % set window length (in both positive and negative directions from behavioral event) in seconds
window_adjusted = round(c_dFF.samplerate * window);                         % express time window in units of photometry signal
center_indx = round(size(zscore_Segs,1)/2);                                        % find index of central point of averaged signal (i.e. at time of behavioral event)
window_start = center_indx - window_adjusted;
window_end = center_indx + window_adjusted;


% calculate peak and/or trough of signal during time window selected for assessment

zspeaks = zeros(1,size(zscore_Segs,2));                                              % peak is maximum signal size
zspeakspre =  zeros(1,size(zscore_Segs,2));                                            % trough value is minimum signal size
for i = 1:size(zscore_Segs,2);
    zspeaks(1,i) = max(zscore_Segs(center_indx:window_end,i));
    zspeakspre(1,i) = mean(zscore_Segs(window_start:center_indx,i));
end;                             

peaks = zeros(1,size(Segs,2));                                              % peak is maximum signal size
peakspre =  zeros(1,size(Segs,2));                                            % trough value is minimum signal size
for i = 1:size(zscore_Segs,2);
    peaks(1,i) = max(Segs(center_indx:(window_end),i));
    peakspre(1,i) = mean(Segs((window_start):center_indx,i));
end;                             


% calcualte AUC during time window

                                           % AUCpre is sum of signal values from the start of the time window up to the behavioral event
AUC = zeros(1,size(zscore_Segs,2));                                           % AUCpost is sum of signal values from the behavioral event up to the end of the time window
for i = 1:size(zscore_Segs,2);
    
    AUC(1,i) = sum(zscore_Segs(center_indx:(window_end),i));
end;


peaks_base = zeros(1,size(zscore_Segs2,2));                                              % peak is maximum signal size
                                          
for i = 1:size(zscore_Segs2,2);
    zspeaks_base(1,i) = mean(zscore_Segs2(center_indx:window_end,i));
    
end;                             

% calculate AUC during time window
% AUCpre is sum of signal values from the start of the time window up to the behavioral event

AUC_base = zeros(1,size(zscore_Segs2,2));                                           
for i = 1:size(zscore_Segs2,2);
    
    AUC_base(1,i) = sum(zscore_Segs2(center_indx:window_end,i));
end;



%% STEP 7: save the data

Trials.NOE = Segs;
Trials.NOEmean = Mean;
Trials.zscoreNOE = zscore_Segs;
Trials.zscoreNOEmean = Mean_zs;

Trials.Base = Segs_base;
Trials.Basemean = Mean_base;
Trials.zscoreBase = zscore_Segs2;
Trials.zscoreBaseMean = Mean_zs_base;
Trials.Timeline = Timeline;
Trials.SEM_zsNOE = nansem(Trials.zscoreNOE,2);
Trials.SEM_zsBase = nansem(Trials.zscoreBase,2);


Stats.PeaksNOE = peaks;
Stats.PeakspreNOE = peakspre;
Stats.AUCNOE = AUC;
Stats.zs_peakNOE = zspeaks;
Stats.zs_peakspreNOE = zspeakspre;
Stats.zs_PeaksBase = zspeaks_base;
Stats.AUCBase = AUC_base;

save(savname, 'Trials', 'Stats')

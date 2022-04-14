function [allpower, avglogspect, fs] = expPower(eegfile, timefile, freqbands, channels, uniqtypes)
% Power calculation for a single file, output average spectrogram across
% all conditions and avg power measures for different trial types in
% discrete frequency bands

% assume there are 3 columns in timefile-- the first two represent the beginning / end
% of each interval, and the third tags the type of interval (e.g. trial vs.
% inter-trial vs. baseline, closed arm vs. center vs. open arm, baseline
% vs. task-social vs. task-non-social, etc.) - eg "107.2067 115.6367 0"

% allpower is a cell array, indexed by channel, trial type, and frequency
% band with the power in each frequency band for every segement analyzed

% avglogspect is an averaged (log) spectrogram - averaged across all
% conditions for each electrode, fs is corresponding frequencies

fid = fopen(timefile);

% making array where each column is time point, first row is starts, second
% row is stops, third row is type

timepts = fscanf(fid, '%f', [3,inf]);
nchannels = length(channels);
nbands = size(freqbands,1);
ntypes = length(uniqtypes);
ntrials = size(timepts,2);

% initialize summary arrays
allpower = cell(nchannels,ntypes,nbands);
avglogspect = cell(nchannels,1);
nspect = zeros(nchannels,1);


% first read in all the data
[data, header] = readedf5(eegfile);
% sampling interval
delt = header.duration / header.nsamples(channels(1));

% number of samples in filter window
fftdur = 2; % 2 second window for FFT
% number of points for filter window
L = round(fftdur/delt); % fftdur = window in seconds, delt = time step,

for channelIDX=1:length(channels),
    % stepping through each trial
    for trialIDX = 1:ntrials,
        % entire length of current event
% %         if timepts(3,trialIDX) == 2
% %             start = round(timepts(1,trialIDX)/delt); % start time converted to samples
% %             stop = round((timepts(1,trialIDX) + 4)/delt);
% %         else
            start = round(timepts(1,trialIDX)/delt); % start time converted to samples 
           % stop = start + round(3/delt);
             stop = round(timepts(2,trialIDX)/delt); % stop time converted to samples
% %         end
        % checks if outside bounds of data, prints warning
        if start < 1 || stop > length(data{channels(channelIDX)}),
            disp('Event outside bounds of data')
            break;
        end
        
        % checks if event is shorter than required filter window, prints
        % warning
        if stop-start < L,
            disp('Event shoter than sampling window-power')
            break;
        end
        
        % spectrogram analysis of power using short time fourier transform
        % L - window for analysis
        % 0 - no overlap between segments
        % 1/delt - sampling rate of data to determine frequencies
        
        % scurr - current spectrogram - rows = frequency, columns = time points
        % each entry is complex number - phase & magnitude 
        
        % fs - frequencies corresponding to rows of scurr/power
        
        
%        d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 59.85, 'halfpowerfrequency2', 60.15, 'samplerate', 2000);
%     data{channels(channelIDX)}= filter(d,data{channels(channelIDX)});
%      d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 19.9, 'halfpowerfrequency2', 20.2, 'samplerate', 2000);
%     data{channels(channelIDX)}= filter(d,data{channels(channelIDX)});
%      d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 39.9, 'halfpowerfrequency2', 40.2, 'samplerate', 2000);
%     data{channels(channelIDX)}= filter(d,data{channels(channelIDX)});
        
        
        [scurr,fs,~,power] = spectrogram(data{channels(channelIDX)}(start:stop),L,L/2,[],1/delt);
        
        % if avg spect has been initiated for that channel, take the log of
        % the spectrogram and add the mean to the current mean
        if ~isempty(avglogspect{channelIDX}),
            avglogspect{channelIDX} = avglogspect{channelIDX} + mean(log(abs(scurr')),1);
        % if not initiated, create avglogspect for that channel with the
        % mean of the log value of scurr
        else
            avglogspect{channelIDX} = mean(log(abs(scurr')),1);
        end
        
        
        Ntmp = size(scurr); % size(scurr,1) = number of frequencies measured
        % nspect - number of spectra computed for that channel
        nspect(channelIDX) = nspect(channelIDX) + Ntmp(1);  
        
        % identify which kind of event it is
        type = find(uniqtypes == timepts(3,trialIDX));
        
        for bandIDX=1:nbands,  % for each frequency band
            % find indicies of all frequences in freqspect between
            % frequencies of interest
            freqin = find(fs > freqbands(bandIDX,1) & fs <= freqbands(bandIDX,2));
            
            % use indices to sum all power within that band
            powers = sum(abs(power(freqin,:)));
            powers = mean(powers);
            % adds single number to cell for that type of event
            % end result: vector in each cell with number of entries
            % corresponding to number of segments that match that type
            allpower{channelIDX,type,bandIDX} = [allpower{channelIDX,type,bandIDX}, powers];
        end
    end
end

% divide sum in avglogspect by total number of spectra
% avglogspect - average for each channel pooling all event types
for channelIDX=1:length(channels),
    avglogspect{channelIDX} = avglogspect{channelIDX} / nspect(channelIDX);
end
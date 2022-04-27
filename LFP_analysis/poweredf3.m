function [allpower, difftypes, avglogspect, fs, t,indvtrials,nspect] = poweredf3(eegfile, timefile, freqbands, channels)

fid = fopen(timefile);

% assume there are 3 columns -- the first two represent the
% beginning / end of each interval, and the third tags the type of
% interval (e.g. trial vs. inter-trial vs. baseline, closed arm vs.
% center vs. open arm, baseline vs. task-social
% vs. task-non-social, etc.)

% allpower is a cell array, indexed by channel, trial type, and
% frequency band, that contains the power within each frequency band
% for each channel and trial type

% difftypes is a vector that contains the different trial types

% avglogspect is an averaged (log) spectrogram, and fs is the
% corresponding vector of frequencies

% number of points / size of window to use for computing FFTs
%L = 2^11;




timepts = fscanf(fid, '%f', [3,inf]);
N = size(timepts);

ntrials = N(2)

fftdur = 2;

difftypes = unique(timepts(3,:));
ntypes = length(difftypes);

% first read in all the data
[data, header] = readedf5(eegfile);

%  data =cellfun(@(x) x/1000,data,'un',0); %making all data in microvolts

delt = header.duration / header.nsamples(1)
fftdur = 2;
% L=2^11;
L = round(fftdur/delt);

M = size(freqbands)

avglogspect = [];
 indvtrials = [];

% Lev = ceil(3*(1/59)/delt+1); %filter order
% B2 = fir1(Lev, [59 61]*delt*2, 'stop'); %

for i=1:length(channels),

    for p=1:length(difftypes),
        for k=1:M(1),
            allpower{i,p,k} = [];
        end
        nspect(i,p) = 0;
        avglogspect{i,p} = [];
        indvtrials{i,p} = [];
        
    end
    
   for currtrial = 1:ntrials,
        
        start = round(timepts(1,currtrial)/delt);
        stop = start + round(4/delt);
%         stop = round(timepts(2,currtrial)/delt);
        
        
        if start < 1 || stop > length(data{channels(i)}),
            continue;
        end
        
        if stop-start < L,
            continue;
        end


        
%        data{channels(i)} = filtfilt(B2,1,data{channels(i)});
%        
%     bandstop filter to remove electrical noise
   d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 59.85, 'halfpowerfrequency2', 60.15, 'samplerate', 2000);
     data{channels(i)}= filter(d,data{channels(i)});
%       d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 19.85, 'halfpowerfrequency2', 20.15, 'samplerate', 2000);
%      data{channels(i)}= filter(d,data{channels(i)});
%       d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 39.85, 'halfpowerfrequency2', 40.15, 'samplerate', 2000);
%      data{channels(i)}= filter(d,data{channels(i)});



        [scurr,fs,t,P] = spectrogram(data{channels(i)}(start:stop),L,L/2,[0:120],1/delt);
        
        p = find(difftypes == timepts(3,currtrial));
            
%         if ~isempty(avglogspect{i,p}),
%             avglogspect{i,p} = avglogspect{i,p} + mean(log(abs(scurr')),1);
%         else
%             avglogspect{i,p} = mean(log(abs(scurr')),1);
%         end
        
        if ~isempty(avglogspect{i,p}),
            avglogspect{i,p} = avglogspect{i,p} + mean(log(abs(scurr')),1);
        else
            avglogspect{i,p} = mean(log(abs(scurr')),1); 
        end
        
       
        
             indvtrials{i,p}(size(indvtrials{i,p},1)+1,:) = mean(log(abs(scurr')),1);
      
        Ntmp = 1;%size(scurr); 
        nspect(i,p) = nspect(i,p) + Ntmp(1);
        
        for k=1:M(1),
            freqin = find(fs > freqbands(k,1) & fs <= freqbands(k,2));

            powers =log(sum(abs(P(freqin,:))));
            
            allpower{i,p,k} = [allpower{i,p,k}, powers];
        end
    end
end

for i=1:length(channels),
    for p=1:length(difftypes),
        avglogspect{i,p} = avglogspect{i,p} / nspect(i,p);
         avglogspect{i,p} = smooth(avglogspect{i,p}); % added smoothin RM 111719
    end
end

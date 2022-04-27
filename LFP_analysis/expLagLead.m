function [avgCovar, avgCovarPhase, avgPLS, avgPlsPhase, avgWpli, avgCoh,avgPli,avgPSI, allCovar, allCovarPhase, allPLS, allPlsPhase, allWpli,allCoh, allPli,allPSI,allCorr, allTimes] = expLagLead_MC(eegfile, timefile,freqbands, channels,uniqtypes)

% Modified by Margaret Cunniff from lagedf5 (Sept 2017)
% Modified RM 12/04/19 to add Cross-corr (allCorr and allTimes)
% Generates covariation, PLS, and WPLI measures
% Outputs both average values for given trial type & all values computed
% for each individual trial


% making array where each column is time point, first row is starts, second
% row is stops, third row is type
fid = fopen(timefile);
timepts = fscanf(fid, '%f', [3,inf]);

ntypes = length(uniqtypes);
ntrials = size(timepts,2);
nbands = size(freqbands,1);
nchannels = length(channels);
nelpairs = nchannels*(nchannels-1)/2;
elpairs = zeros(nelpairs,2);

% initialize summary arrays - average of all trials for a given condition
avgPLS = zeros(nelpairs, ntypes, nbands);
avgPlsPhase = zeros(nelpairs, ntypes, nbands);
avgCovar = zeros(nelpairs,ntypes, nbands);
avgCovarPhase = zeros(nelpairs,ntypes, nbands);
avgWpli = zeros(nelpairs, ntypes, nbands);
avgCoh = zeros(nelpairs, ntypes, nbands);
avgPli = zeros(nelpairs, ntypes, nbands);
avgPSI = zeros(nelpairs, ntypes, nbands);


% Individual data points for every trial analyzed
allCovar = cell(nelpairs, ntypes,nbands);
allCovarPhase = cell(nelpairs, ntypes,nbands);
allPLS = cell(nelpairs,ntypes,nbands);
allPlsPhase = cell(nelpairs,ntypes,nbands);
allWpli = cell(nelpairs,ntypes,nbands);
allCoh= cell(nelpairs,ntypes,nbands);
allPli = cell(nelpairs,ntypes,nbands);
allPSI = cell(nelpairs,ntypes,nbands);
allCorr = cell(nelpairs,ntypes,nbands);
allTimes = cell(nelpairs,ntypes,nbands);

% first read in all the data
[data, header] = readedf5(eegfile);

delt = header.duration / header.nsamples(channels(1));

% Event - individual entry for each calculation made - for each trial,
% events for every freq band/electrode pair/type
event = 0;
covar = [];
covarPhase = [];
pls = [];
plsphase = [];
wpli = [];
coh = [];
corr = [];
times = [];
pli = [];
PSI = [];

for ch1IDX=1:length(channels),
    for ch2IDX=1+ch1IDX:length(channels),
        % stepping through all possible electrode pairs, 
        for bandIDX=1:nbands,
            % individual calculations for each frequency band
            fmin = freqbands(bandIDX,1);
            fmax = freqbands(bandIDX,2);
            
                
            for trialIDX = 1:ntrials,
                % Stepping through each trail from time file
                start = round(timepts(1,trialIDX)/delt); % converting start & stop to sample numbers
               
                    stop = start + round(3.5/delt);
              
%                  stop = round(timepts(2,trialIDX)/delt);
               
                
                % Throw error if times invalid or too small for analysis
                % min length determined by length needed to create filter
                minEventLength = ceil(3*(1/fmin)/delt+1);
                
                if start < 1 || stop > length(data{channels(ch1IDX)}),
                    disp('Time outside bounds of data')
                    break;
                elseif stop-start < minEventLength, 
                    disp('Event smaller than sampling window')
                    disp(stop-start)
                    break;
                end
%                 
%                 d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 59.85, 'halfpowerfrequency2', 60.15, 'samplerate', 2000);
%     data{channels(ch1IDX)}= filter(d,data{channels(ch1IDX)});
%      d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 19.9, 'halfpowerfrequency2', 20.2, 'samplerate', 2000);
%     data{channels(ch1IDX)}= filter(d,data{channels(ch1IDX)});
%      d= designfilt('bandstopiir','filterorder', 10, 'Halfpowerfrequency1', 39.9, 'halfpowerfrequency2', 40.2, 'samplerate', 2000);
%     data{channels(ch1IDX)}= filter(d,data{channels(ch1IDX)});
%         
                % create reference arrays for each event with electrode
                % IDs, frequency band, and segment type
                event = event + 1; 
                
                % Event # acts as index
                el1(event) = channels(ch1IDX) ;% which channel is assigned as el1
                el2(event) = channels(ch2IDX); % which channel is assigned as el2
                freqband(event) = bandIDX; % which frequency band is analyzed
                type(event) = timepts(3,trialIDX); % what type is the segment
                el1_data = data{channels(ch1IDX)}(start:stop);
                el2_data = data{channels(ch2IDX)}(start:stop);
                
                [covar(event),covarPhase(event),pls(event),plsphase(event),wpli(event),coh(event),corr{event},times{event}, pli(event), PSI(event)] = ...
                    laglead(el1_data,el2_data,delt,fmin,fmax);
            end
        end
    end
end

try 
    length(el1);
catch
    error('No data generated by trialLagLead.')
end

% Now have several large vectors containing all of data and reference information
% Organizing into summary/output tables
elpair = 0;
for ch1IDX=1:length(channels),
    for ch2IDX=ch1IDX+1:length(channels),
        elpair = elpair+1;
        for bandIDX=1:nbands,
            for typeIDX=1:ntypes,
                % which electrode pair is being analyzed
                elpairs(elpair,:) = [ch1IDX ch2IDX];
                % find all values matching current electrodes, band, and
                % type
                findTrials = find(el1 == channels(ch1IDX) & el2 == channels(ch2IDX) & freqband == bandIDX & type == uniqtypes(typeIDX));
                if isempty(findTrials),
                    disp('No trials of that type.')
                    disp([ch1IDX, ch2IDX, bandIDX, typeIDX])
                    break;
                end
                
                % go through all the trials of that type and add them to
                % summary vectors for average values
                for trials=1:length(findTrials),
                    avgCovar(elpair,typeIDX, bandIDX) = avgCovar(elpair,typeIDX,bandIDX) + covar(findTrials(trials));
                    avgCovarPhase(elpair,typeIDX, bandIDX) = avgCovarPhase(elpair,typeIDX,bandIDX) + covarPhase(findTrials(trials));
                    avgPLS(elpair,typeIDX, bandIDX) = avgPLS(elpair,typeIDX, bandIDX) + pls(findTrials(trials));
                    avgWpli(elpair,typeIDX, bandIDX) = avgWpli(elpair,typeIDX, bandIDX) + wpli(findTrials(trials)) ;
                    avgPlsPhase(elpair,typeIDX, bandIDX) = avgPlsPhase(elpair,typeIDX, bandIDX) + plsphase(findTrials(trials));
                    avgCoh(elpair,typeIDX, bandIDX) = avgCoh(elpair,typeIDX, bandIDX) + coh(findTrials(trials));
                    avgPli(elpair,typeIDX, bandIDX) = avgPli(elpair,typeIDX, bandIDX) + pli(findTrials(trials));
                     avgPSI(elpair,typeIDX, bandIDX) = avgPSI(elpair,typeIDX, bandIDX) + PSI(findTrials(trials));
                end
                
                % Divide summed values by number of trials of that type for
                % individual mean for each condition 
                avgCovar(elpair, typeIDX, bandIDX) = avgCovar(elpair,typeIDX,bandIDX) / length(findTrials);
                avgPLS(elpair, typeIDX, bandIDX) = avgPLS(elpair,typeIDX,bandIDX) / length(findTrials);
                avgWpli(elpair, typeIDX, bandIDX) = avgWpli(elpair,typeIDX,bandIDX) / length(findTrials);
                avgPlsPhase(elpair,typeIDX, bandIDX) = avgPlsPhase(elpair,typeIDX, bandIDX) /length(findTrials);
                avgCoh(elpair,typeIDX, bandIDX) = avgCoh(elpair,typeIDX, bandIDX) /length(findTrials);
                avgPli(elpair,typeIDX, bandIDX) = avgPli(elpair,typeIDX, bandIDX) /length(findTrials);
                %avgCorr(elpair,typeIDX, bandIDX) = avgCorr(elpair,typeIDX, bandIDX) /length(findTrials);
                
                % add all individual values to cell arrays containing all
                % data
                for trials=1:length(findTrials),
                    allCovar{elpair,typeIDX, bandIDX} = [allCovar{elpair,typeIDX, bandIDX}, covar(findTrials(trials))];
                    allCovarPhase{elpair,typeIDX, bandIDX} = [allCovarPhase{elpair,typeIDX, bandIDX}, covarPhase(findTrials(trials))];
                    allPLS{elpair,typeIDX, bandIDX} = [allPLS{elpair,typeIDX, bandIDX}, pls(findTrials(trials))];
                    allPlsPhase{elpair,typeIDX, bandIDX} = [allPlsPhase{elpair,typeIDX, bandIDX}, plsphase(findTrials(trials))];
                    allWpli{elpair,typeIDX, bandIDX} = [allWpli{elpair,typeIDX, bandIDX}, wpli(findTrials(trials))];
                    allCoh{elpair,typeIDX, bandIDX} = [allCoh{elpair,typeIDX, bandIDX}, coh(findTrials(trials))];
                    allPli{elpair,typeIDX, bandIDX} = [allPli{elpair,typeIDX, bandIDX}, pli(findTrials(trials))];
                    allPSI{elpair,typeIDX, bandIDX} = [allPSI{elpair,typeIDX, bandIDX}, PSI(findTrials(trials))];
                    allCorr{elpair,typeIDX, bandIDX} = [allCorr{elpair,typeIDX, bandIDX}, corr{findTrials(trials)}];
                    
                    allTimes{elpair,typeIDX, bandIDX} = [allTimes{elpair,typeIDX, bandIDX}, times{findTrials(trials)}];
                end
            end
        end
    end
end

% Lags = linspace(-125, 125, 501);
% Lags = Lags';
% Corr = squeeze(allCorr);
% 
% X = Corr{1,1};
% Y = Corr{2,1};
% % Z = Corr{3,1};
% avgX = mean (X,2);
% avgY = mean (Y,2);
% % avgZ = mean (Z,2);
% 
% SEMX = std(X, 0, 2)/sqrt(size(X,2));
% %  ts = tinv([0.025  0.975],size(X,2)-1);
% %  CIX = ts(1,2)*SEMX;
% SEMY = std(Y, 0, 2)/sqrt(size(Y,2));
% % ts = tinv([0.025  0.975],size(Y,2)-1);
% %  CIY = ts(1,2)*SEMY;

% SEMZ = std(Z, 0, 2)/sqrt(size(Z,2));

% out = [Lags, avgX, SEMX, avgY, SEMY];
% csvwrite(strcat('out_', mouseID,'.csv'), out);
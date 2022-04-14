function [allPower, allCovar, allCovarPhase, allPLS, allPlsPhase, allWpli, ...
    allCoh,expLogPower, expCovar, expCovarPhase, expPLS, expPlsPhase, expWpli, ...
    expCoh,allPowerTypes, allPowerExps, allPowerTreatments, allPairTypes, allPairExps, allPairTreatments] = ...
    pairElectrodeBlockAnalysis_MC(filelist, fileprefix, channelIDs, types, treatments, animals, electrodeIDs)
% Returns power and coherence measures taken over blocks of time as specified by the user 
% Writes csv file for each measure
% Can be expanded to multiple electrodes but currently output optimized for
% 2 electrodes/1 pair, such as PFC/HPC recordings

% Modfied by Margaret Cunniff from lagleadscript (Sept 2016)
% Input arguments: 
    % filelist - txt file including EDF file & associated annotation file
    % fileprefix - prefix specifying name/location of output csv files
    % channelIDs - vector of recording electrodes for each animal
    % types - types includes in annotation files (eg 0:3)
    % treatments - vector for grouping variable, eg Virus or Genotype
    % animals - vector of IDs for all animals included
    % electrodeIDs - locations of electrode for output file labeling
    
    % EXAMPLE

%     filelist = 'pogZ_WT_EPM_blockAnnotations2.txt';
%     fileprefix = 'LFPanalysisTest_';
%     channelIDs = [2:3; 2:3; 2:3];
%     types = 0:1; % blocks
%     animals = [1, 2, 3];
%     treatments = ones(size(animals));
%     electrodeIDs = {'Hpc';'Pfc'};
    
    % SETUP - ESTABLISHING RESULTS VECTORS, ETC
    
    freqbands =[20 55; 65 90]; % frequency bands to use for analysis [4 12; 12 20;
    nbands = size(freqbands,1);
    nchannels = size(channelIDs,2); % number of recording electrodes
    nelpairs = nchannels*(nchannels-1)/2;
    ntypes = length(types); 
    nanimals = length(animals);
   
    % OUTPUT VARIABLES
    % allX - pooled individual data points from all experiments w/
    % accompanying arrays to describe type & experiment - ANOVA input
        % eg 1x10 array with all values, 1x10 array identifying which type
        % of block each data point belongs to, 1x10 array identifying which
        % experiment each data point came from, etc
    %expX - mean for each file for each condition

    % Individual Segment Values
    allPower = cell(nchannels, nbands);
    allCovar = cell(nelpairs, nbands);
    allCovarPhase = cell(nelpairs, nbands);
    allPLS = cell(nelpairs, nbands);
    allPlsPhase = cell(nelpairs, nbands);
    allWpli = cell(nelpairs, nbands);
    allCoh = cell(nelpairs, nbands);
    

    % References for grouping variables for power/single electrode measures
    allPowerTypes = cell(nchannels, nbands);
    allPowerExps = cell(nchannels, nbands);
    allPowerTreatments = cell(nchannels,nbands);
    
    %Grouping variables for all electrode pair measures - covar, PLS, WPLI
    allPairTypes = cell(nelpairs, nbands);
    allPairExps = cell(nelpairs, nbands);
    allPairTreatments = cell(nelpairs, nbands);

    % Mean measures for each experiment
    expPlsPhase = [];
    expCovarPhase = [];
    expCovar = [];
    expPLS = [];
    expWpli = [];
    expCoh = [];
    expLogPower = cell(nchannels, ntypes, nbands);

    % CALCULATING FILE MEASURES
    % stepping through files in the filelist
    nfile = 0;
    fid = fopen(filelist);
    eegfile = fscanf(fid, '%s', 1);

    while eegfile,
        nfile = nfile+1;
        disp(' ')
        disp(eegfile)
        timefile = fscanf(fid, '%s', 1);
        channels = channelIDs(nfile,:);
        
        if ~timefile,
            break;
        end

        disp('Calculating power')
        filePower = expPower(eegfile, timefile, freqbands, channels, types);
        % filePower - cell array, indexed by channel, trial type, and band
        % each cell contains a vector with average power for each segment
        % that matches the trial type (eg {1,1,1} will contain 1xN vector
        % where each entry is the average power of electrode 1 in the theta
        % band for a single segment in trials of type 1 and N is the
        % number of such segments
        
        % cell arrays of single electrode trial data organized by condition
        for channelIDX=1:nchannels,
            for bandIDX=1:nbands,
                for typeIDX=1:ntypes,
                    % log of all power values for condition
                    % 1xN array, N = # of segments meeting condition
                    fileLogPower = log(filePower{channelIDX,typeIDX,bandIDX});
                    
                    % Combining all log values from all files together in
                    % single vector - 1xN, = # of segments pooled across
                    % all files
                    allPower{channelIDX,bandIDX} = [allPower{channelIDX,bandIDX}, fileLogPower];
                    
                    % Reference arrays - same size as vectors in allPower
                    % noting Exp #, type, & treatment. Used to
                    % differentiated pooled values in allPower for stats
                    allPowerExps{channelIDX,bandIDX} = [allPowerExps{channelIDX,bandIDX}, nfile*ones(1,length(filePower{channelIDX,typeIDX,bandIDX}))];
                    allPowerTypes{channelIDX,bandIDX} = [allPowerTypes{channelIDX,bandIDX}, types(typeIDX)*ones(1,length(filePower{channelIDX,typeIDX,bandIDX}))];
                    allPowerTreatments{channelIDX,bandIDX} = [allPowerTreatments{channelIDX,bandIDX}, treatments(nfile)*ones(1,length(filePower{channelIDX,typeIDX,bandIDX}))];
                    
                    % Cell array where each cell is vector with mean log
                    % power for each file for that condition - 1xN
                    % array, N = number of files
                    expLogPower{channelIDX,typeIDX,bandIDX} = [expLogPower{channelIDX,typeIDX,bandIDX} mean(fileLogPower)]; 
                end
            end
        end


        disp('Calculating coherence measures')
        [fileCovar, fileCovarPhase, filePLS, filePlsPhase, fileWpli, fileCoh, trialCovar,trialCovarPhase, trialPLS, trialPlsPhase, trialWpli, trialCoh] =...
            expLagLead(eegfile, timefile, freqbands, channels, types);
        % fileX - vector of average values for given experiment for each trial type
        % trialX - cell arrays containing all individual values for every trial
        
        % covar - amplitude covariation
        % PLS - phase locking statistic (summing all phase diffs)
        % WPLI - weighted phase locking index
     
        % Vectors of average values from each individual file
        expPlsPhase = [expPlsPhase filePlsPhase];
        expCovar = [expCovar fileCovar];
        expCovarPhase = [expCovarPhase fileCovarPhase];
        expPLS = [expPLS filePLS];
        expWpli = [expWpli fileWpli];
        expCoh = [expCoh fileCoh];
     

        % cell arrays of electrode pair individual trial data organized by condition
        for elpairIDX=1:nelpairs,
            for bandIDX=1:nbands,
                for typeIDX=1:ntypes,
                    %creating summary variables - Lag/Lead
                    % combine all data together & create references for
                    % experiment and event type
                    allCovar{elpairIDX,bandIDX} = [allCovar{elpairIDX,bandIDX}, trialCovar{elpairIDX,typeIDX,bandIDX}];
                    allCovarPhase{elpairIDX,bandIDX} = [allCovarPhase{elpairIDX,bandIDX}, trialCovarPhase{elpairIDX,typeIDX,bandIDX}];
                    allPLS{elpairIDX,bandIDX} = [allPLS{elpairIDX,bandIDX}, trialPLS{elpairIDX,typeIDX,bandIDX}];
                    allPlsPhase{elpairIDX,bandIDX} = [allPlsPhase{elpairIDX,bandIDX}, trialPlsPhase{elpairIDX,typeIDX,bandIDX}];
                    allWpli{elpairIDX,bandIDX} = [allWpli{elpairIDX,bandIDX}, trialWpli{elpairIDX,typeIDX,bandIDX}];
                    allCoh{elpairIDX,bandIDX} = [allCoh{elpairIDX,bandIDX}, trialCoh{elpairIDX,typeIDX,bandIDX}];

                    % grouping arrays for posthoc ID of experiments/types
                    % for all pair measures
                    allPairExps{elpairIDX,bandIDX} = [allPairExps{elpairIDX,bandIDX}, nfile*ones(1,length(trialCovar{elpairIDX,typeIDX,bandIDX}))]; % exp - N = which file
                    allPairTypes{elpairIDX,bandIDX} = [allPairTypes{elpairIDX,bandIDX}, types(typeIDX)*ones(1,length(trialCovar{elpairIDX,typeIDX,bandIDX}))];   
                    allPairTreatments{elpairIDX,bandIDX} = [allPairTreatments{elpairIDX,bandIDX}, treatments(nfile)*ones(1,length(trialCovar{elpairIDX,typeIDX,bandIDX}))];
                end
            end
        end
        
        % Advance to next EEG file from file list
        eegfile = fscanf(fid, '%s', 1);
    end

    % reformatting output matricies to more easily copy into spreadsheet
    % only valid for 2 electrodes

    % power labeling assumes 1st el = HPC, 2nd = PFC
    % creating series of tables with columns for each frequency band and
    % animal and type labeled -> can be fed into pandas/copied to spreadsheet 
   
    % creating column vectors for labeling type & animal for each measure
    % repeated so that all measures from one animal are in consecutive
    % rows, same types spread out across animals
    types = reshape(types,ntypes,1);
    typeArray = repmat(types,nanimals,1);
    animalArray = zeros(nanimals*ntypes,1);
    treatmentArray = zeros(nanimals*ntypes,1);
    z = 1;
    for animalIDX=1:nanimals
        for typeIDX = 1:ntypes
            animalArray(z,1) = animals(animalIDX);
            treatmentArray(z,1) = treatments(animalIDX);
            z = z+1;
        end
    end
     % EG if Types are 0,1 and animals are A, B:
     % typeArray = 0 1 0 1
     % animal Array = A A B B 
   
    % writing average values for each animals to CSV files
    % 6 column array - Animal, Type, Theta, Beta, Low Gamma, High Gamma
    % squeeze coverts 3D output array to 2D (first dim is # of elpairs,
    % only one in the PFC/HPC case)
    expCovar = [animalArray typeArray treatmentArray squeeze(expCovar)];
    csvwrite(strcat(fileprefix,'avgCovar.csv'), expCovar);

    expCovarPhase = [animalArray typeArray treatmentArray squeeze(expCovarPhase)];
    csvwrite(strcat(fileprefix,'avgCovarPhase.csv'), expCovarPhase);

    expPLS = [animalArray typeArray treatmentArray squeeze(expPLS)];
    csvwrite(strcat(fileprefix,'avgPLS.csv'), expPLS);

    expWpli = [animalArray typeArray treatmentArray squeeze(expWpli)];
    csvwrite(strcat(fileprefix,'avgWpli.csv'), expWpli);

    expPlsPhase = [animalArray typeArray squeeze(expPlsPhase)];
    csvwrite(strcat(fileprefix,'avgPlsPhase.csv'), expPlsPhase);
    
    expCoh = [animalArray typeArray treatmentArray squeeze(expCoh)];
    csvwrite(strcat(fileprefix,'avgCoh.csv'), expCoh);

    % same concept for power but input data formatted differently 
    % seperating out power measures from different electrodes 
    pwrEl1 = expLogPower(1,:,:);
    pwrEl2 = expLogPower(2,:,:);

    el1Array = squeeze(cell2mat(pwrEl1));
    el2Array = squeeze(cell2mat(pwrEl2));

    pwrAnimals = repmat(animals',ntypes,1);
    pwrTreatments = repmat(treatments',ntypes,1);
    pwrTypes = [];

    m=0;
    for k=1:ntypes
        pwrTypes = [pwrTypes; m*ones(nanimals,1)];
        m = m +1;
    end
    
    el1Pwr = [pwrAnimals pwrTypes pwrTreatments el1Array];
    csvwrite(strcat(fileprefix,'avg',electrodeIDs{1},'Power.csv'), el1Pwr);

    el2Pwr = [pwrAnimals pwrTypes pwrTreatments el2Array];
csvwrite(strcat(fileprefix,'avg',electrodeIDs{2},'Power.csv'), el2Pwr);
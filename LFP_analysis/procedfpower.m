function [avgpowers, avglogpowers, allpowers, expnumb, exptype, ...
          avglogspect, fs, pwr, pwrtypes, pwrexp, ppwr, ppwrind,corr1, corr2, indlogspect,allspects]  = procedfpower(filelist, freqbands, channels)

% INPUT ARGUMENTS

% filelist is the list of files to analyze

% freqbands is a matrix; each row has two elements which define a
% range of frequencies (in Hz) to analyze

% channels is a vector of the channels to analyze


% OUTPUT ARGUMENTS

% avgpowers is a cell array indexed by channel and frequency band.
% Each element is a matrix, in which each row corresponds to a
% different trial type, and each column corresponds to a different
% experiment

% avglogpowers is similar, but it the average of the log powers

% allpowers is a 2-dim cell, also indexed by channel and frequency
% band. Each element contains a vector which are power measurements
% for that channel / frequency band from a particular epoch type
% and experiment.  expnumb and exptype are corresponding arrays
% which contain the experiment number and epoch type respectively

% avglogspect is a 2-dim cell array, indexed by channel and trial
% type.  Each element is a spectrogram, averaged across all
% experiments.  fs is the corresponding vector of frequencies


fid = fopen(filelist);


filename = fscanf(fid, '%s', [1]);

n=0;

ntypes=0:2;
alltypes = length(ntypes);

M=size(freqbands);


for i=1:length(channels),
    for j=1:M(1),
        pwr{i,j} = [];
        pwrtypes{i,j} = [];
        pwrexp{i,j} = [];
        ppwr = [];
        ppwrind = []
    end
end

for i=1:length(channels),
    for j=1:M(1),
        allpowers{i,j} = [];
        expnumb{i,j} = [];
        exptype{i,j} = [];
    end
end

avglogspect = {};
indlogspect = {};
allspects = {};

while filename,
    filename
    timefile = fscanf(fid, '%s', [1]);

    n = n+1;
    [allpower,difftypes,avglogspecttmp,fs,t,indvtrials] = poweredf3(filename, timefile, freqbands, channels);
    
    N = size(allpower);
    
    for ch=1:N(1),
        for type=1:N(2),
            indlogspect{ch,type}(:,n) = 1*avglogspecttmp{ch,type}; %Rm replaced 2000 with 1 for multiplication factor
 
            if n ==1;
                size(avglogspecttmp{ch,type});
                avglogspect{ch,type} = avglogspecttmp{ch,type};
                nspect(ch,type) = 1;
                allspects{ch,type} = indvtrials{ch,type};
                
            else
                avglogspect{ch,type} = avglogspect{ch,type} + avglogspecttmp{ch,type};
                nspect(ch,type) = nspect(ch,type)+1;
                allspects{ch,type} = vertcat(allspects{ch,type},indvtrials{ch,type}) ;
                
            end
            
           
            for band=1:N(3),
                avgpowers{ch,band}(type,n) = mean(allpower{ch,type,band});
                avglogpowers{ch,band}(type,n) = mean((allpower{ch,type,band}));
                allpowers{ch,band} = [allpowers{ch,band}, allpower{ch,type,band}];
%                 avglogspect{ch,band}=log(allpowers{ch,band});
                expnumb{ch,band} = [expnumb{ch,band}, n*ones(1, ...
                                                             length(allpower{ch,type,band}))];
                exptype{ch,band} = [exptype{ch,band}, difftypes(type)*ones(1,length(allpower{ch,type,band}))];
                
                
                tmp = (allpower{ch,type,band});
                pwr{ch,band} = [pwr{ch,band}, tmp];
                pwrexp{ch,band} = [pwrexp{ch,band}, n*ones(1,length(allpower{ch,type,band}))];
                pwrtypes{ch,band} = [pwrtypes{ch,band}, difftypes(type)* ones(1,length(allpower{ch,type,band}))];
            end
        end
    end
    filename = fscanf(fid, '%s', [1]);
end

M = size(avglogspect);
for i=1:M(1),
    for j=1:M(2),
        %avglogspect{i,j} = smooth(2000*avglogspect{i,j})/ nspect(i,j);
        avglogspect{i,j} = smooth(1*avglogspect{i,j})/ nspect(i,j); %RM added smooth ; %Rm replaces 2000 with 1 (multiplication factor)
        
    end
end
PFC = (pwr{1,1})';
DHC = (pwr{2,2})';
types = (pwrtypes{1,1})';
new= [types PFC DHC];
ind1 = new(:,1) == 0;
ind2 = new(:,1) == 2;
A1 = new(ind1,:);
A2 = new(ind2,:);
corr1 = corrcoef(A1(:,2),A1(:,3))
%corr1 = corr1(1,2)
corr2 = corrcoef(A2(:,2),A2(:,3))


%%
% now for each frequency band and electrode, do ANOVA on all the
% power spectrum data
M=size(freqbands);
 for ch=1:N(1),
   for j=1:M(1)
%for i=1:length(channels),
 %   for j=1:M,
        [p,tbl] = anovan(pwr{ch,j}, {pwrexp{ch,j}, pwrtypes{ch,j}}, [1 0; 0 1]);
        % p-value for a main effect of interval type on power in
        % this frequency band
        ppwr(ch,j) = p(2);
        
        % now for each pair of conditions do the t-tests
%         pwr2{ch,j} = exp(pwr{ch,j});
        n = 0;
%         for k=1:alltypes,
%             for l=k+1:alltypes,
%                 n = n+1;
%                 ppwrind{ch,j}(n) = ranksum(pwr{ch,j}(find(pwrtypes{ch,j} == ...
%                                                      ntypes(k))), pwr{ch,j}(find(pwrtypes{ch,j} == ntypes(l))));
%                  ppwrtypes{ch,j}(n,1:2) = [difftypes(k), ntypes(l)];
%             end
%         end
    end
end
% 
%%
%%confidence interval added 081820 by RM
in = find(fs>57 & fs<63);
N = size(indlogspect);
spectSEM = {};
yCI95 = {};
for ch=1:N(1),
        for type=1:N(2),
            %indlogspect{ch,type} = indlogspect{ch,type}.*2000;
            
            a = size(indlogspect{ch,type},2);
            avg_ch{ch,type} = mean(indlogspect{ch,type},2);
            avg_ch{ch,type}(in) = NaN;
            spectSEM{ch,type} = (std(indlogspect{ch,type}'))/sqrt(a);
          
            
            CI95 = tinv([0.05  0.95], a-1); 
            yCI95{ch,type} = bsxfun(@times, spectSEM{ch,type}, CI95(:));
            spectSEM{ch,type} = spectSEM{ch,type}';
            %yCI95{ch,type} = mean(indlogspect{ch,type}') + (CI95*spectSEM{ch,type}); 
            
            yCI95{ch,type} = yCI95{ch,type}';
        end
end

 
 
% 
%     avglogspect{1,1}(in) = NaN;avglogspect{1,2}(in) = NaN;avglogspect{1,3}(in) = NaN;
% avglogspect{2,1}(in) = NaN;avglogspect{2,2}(in) = NaN;avglogspect{2,3}(in) = NaN;

%%
figure
    
plot (fs, avglogspect{1,1}, 'k','LineWidth',2)
hold
plot (fs, avglogspect{1,2}, 'b','LineWidth',2)

% plot (fs(1:150), avglogspect{1,3}(1:150), 'g','LineWidth',2)
%plot (fs(1:200), avglogspect{1,4}(1:200), 'r','LineWidth',2)
title(['Power Spectrum, PFC'])
    ylabel('Average Log Power')
    xlabel('Frequency (Hz)')
% 
%     
figure
    
plot (fs, avglogspect{2,1}, 'k','LineWidth',2)
hold
plot (fs, avglogspect{2,2}, 'b','LineWidth',2)

% plot (fs(1:150), avglogspect{2,3}(1:150), 'g','LineWidth',2)
%plot (fs(1:200), avglogspect{2,4}(1:200), 'r','LineWidth',2)
%plot (fs(1:200), avglogspect{2,5}(1:200), 'c','LineWidth',2)
title(['Power Spectrum, DHC'])
    ylabel('Average Log Power')
    xlabel('Frequency (Hz)')

    
    
    
% % figure  
% % plot (fs(1:150), avglogspect{3,1}(1:150), 'k','LineWidth',2)
% % hold
% % plot (fs(1:150), avglogspect{3,2}(1:150), 'r','LineWidth',2)
% % plot (fs(1:150), avglogspect{3,3}(1:150), 'b','LineWidth',2)
% % plot (fs(1:150), avglogspect{3,4}(1:150), 'g','LineWidth',2)
% % title(['Power Spectrum, VHC'])
% %     ylabel('Average Log Power')
% %     xlabel('Frequency (Hz)')
% % 
% 
% PFC = (pwr{1,1})';
% DHC = (pwr{2,5})';
% types = (pwrtypes{1,1})';
% new= [types PFC DHC];
% ind1 = new(:,1) == 0;
% ind2 = new(:,1) == 2;
% A1 = new(ind1,:);
% A2 = new(ind2,:);
% corr1 = corrcoef(A1(:,2),A1(:,3))
% %corr1 = corr1(1,2)
% corr2 = corrcoef(A2(:,2),A2(:,3))
% %corr2 = corr2(1,2)
%print 'corr2' [a
%%
for i = 1:size(allspects,1)
    for j = 1:size(allspects,2)
        indlogspect{i,j} = indlogspect{i,j};
   avg_ch{i,j} =  mean(allspects{i,j},1)'; 
    avg_ch{i,j} = movmean(avg_ch{i,j},10);
   stdev_ch{i,j} = std(allspects{i,j},1)';
sem_ch{i,j} = stdev_ch{i,j}/sqrt(size(allspects{i,j},1));
   sem_ch{i,j} = movmean(sem_ch{i,j},10);
% avg_ch{i,j} = avg_ch{i,j}';
% sem_ch{i,j} = sem_ch{i,j}';

% in = find(fs>57 & fs<63);
% 
%     avg_ch{i,j}(in) = NaN;
%     sem_ch{1,j}(in) = NaN;
    end
end

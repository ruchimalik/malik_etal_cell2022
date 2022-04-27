%function [maxval,lagtime,pls,plsphase,wpli] = laglead6(x,y,dt,fmin,fmax)
function [maxval,lagtime,pls,plsphase,wpli,coh,Ctmp,lags, pli, PSI] = laglead(x,y,dt,fmin,fmax)


% modified by RM on 5/21/19 -- added lagtime estimate

% first calculate the central frequency
fmid = 0.5*(fmin+fmax);

fftdur = 2; % 1 second window for mscoher
% number of points for filter window
L = round(fftdur/dt);

% the filter length should be 3 times the period corresponding to
% the minimum frequency
L = ceil(3*(1/fmin)/dt+1);

% now filter
B = fir1(L, [fmin fmax]*dt*2);

x2 = filtfilt(B,1,x);
y2 = filtfilt(B,1,y);
win = 1/(4*dt);
fseries = fmin:1:fmax;
[cxy,f] = mscohere(x2,y2,win,0.8*win,fseries, 1/dt);

coh = mean (cxy);

% now compute the hilbert transforms of both filtered signals and
% extract the amplitude envelopes

x3 = hilbert(x2);
y3 = hilbert(y2);

ampx = abs(x3);
ampy = abs(y3);

phasex = angle(x3);
phasey = angle(y3);

ampx2 = ampx- mean(ampx);
ampy2 = ampy - mean(ampy);

period = round((1/fmid)/dt);
rp2 = round(period/2);

m = 0;
for j=-rp2:rp2,
    m = m+1;
    if j < 0,
        acx(m) = sum(ampx2(1:length(ampx2)+j).* ampx2(1:length(ampx2)+j));
        acy(m) = sum(ampy2(-j+1:length(ampy2)).* ampy2(-j+1:length(ampy2)));
    else
        acx(m) = sum(ampx2(j+1:length(ampx2)) .* ampx2(j+1:length(ampx2)));
        acy(m) = sum(ampy2(1:length(ampy2)-j).* ampy2(1:length(ampy2)-j));
    end
end

[Ctmp,lags] = xcorr(ampx2,ampy2,rp2); %changed rp2 to period
%plot (lags, Ctmp)
lags= lags*dt*1000;%converting lags to ms 
%Ctmp = xcorr(ampx2,ampy2,rp2);

Ctmp = Ctmp ./ sqrt(acx'.*acy');
%plot (Ctmp)
%[maxval,lag] = max(Ctmp); % modified by RM 061119
[maxval,I] = max(Ctmp);

lagtime= lags(I); % time at which amp cross-corr. is max (in ms)

%converting lag index to actual time in ms
% [m,n] = size(Ctmp); 
% t = linspace(-(1000/fmid),(1000/fmid),m);
% lagtime = t(lag); % lag time in ms where Ctmp is max

tmp = exp(sqrt(-1)*(phasex-phasey)); % differences between phase at each time point
plstmp = (mean(tmp)); % mean of phase difference - normal PLS
plsphase = angle(plstmp); % angle of PLS


pls = abs(plstmp); % amplitude of PLS

X = ampx.*ampy.*tmp; % multiplying all phase differences by amplitude at that point
imX = imag(X);% taking the imaginary component
wpli = abs(mean(abs(imX) .* sign(imX))) / mean(abs(imX)); % 

% plitmp = phasex-phasey;
pli = mean(sign(imX)); %mean(sign(plitmp)); %Phase Lag Index

% PSItmp = X/(sqrt((ampx.ampx).*(ampy.ampy)));
denom = (sqrt((ampx.*ampx).*(ampy.*ampy))); %https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.20346
PSItmp = X./denom;
PSI = (imag(PSItmp));
PSI= mean(PSI); % Imaginary part of coherency
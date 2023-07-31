function [lat,coeff,rs,converge,yf] = latencyfit4AM(y,t,peak,fig)

%Fits a double Gaussian + cumulative Gaussian to neural data. latency is
%taken as the point at which the fit reaches 33% of its maximum.

%Input args:
%y = row vector of neural data, can be raw MUA or difference
%t = row vector of time in ms (same length as y)
%peak: 1 = chooses starting guesses appropriate for raw MUA (default), 2 (or any other value) = starting values for modulation curves
%fig = 1 (open a figure and show the fit), 0 = no figure (default)

%Output args
%lat = latency at 33% of the peak
%coeff = fit coefficients
%rs = adjusted r-squared value of fit
%converge - did fit converge?

%M.W.Self 2014

%Checks
if nargin < 3
    peak = 1;
end
if nargin < 4
    fig = 0;
end
if size(y,1)>size(y,2)
    y = y';
end
if size(t,1)>size(t,2)
    t = t';
end

if peak == 1
    %These are the limits and initial guesses of the fit.
    %These can be used to restrain the fit to exclude unlikely values and
    %prevent fitting of noise. 
    
    %Note also that the initial guesses for the gains will depend on your
    %units. This was designed to fit normalised MUA data, i.e. the data
    %falls mostly between 0 and 1.
    
    %The second Gaussian will sometimes fit noise components, but this
    %generally has little impact on the overall goodness of fit or the
    %measured latency.
    
    %PEAK FITTING LIMITS
    %1t Gaussian, mean = 20-100ms [50ms]
    %2nd Gaussian, mean = 40-200ms [100ms]
    %CumGauss, mean = 50-200ms [100ms]
    %first Gaussian
    g1g = [0.001 0.015 0.5];    %Gain of first Gaussian [min mean max]
    g1m = [0.02 0.05 0.1];      %Mean of first Gaussian
    g1s = [0.001 0.01 0.1];     %STD of first Gaussian
    %second Gaussian
    g2g = [0 0.01 0.5];         %Gain of second Gaussian
    g2m = [0.04 0.1 0.2];       %Mean of second Gaussian
    g2s = [0.001 0.02 0.2];     %STD of second Gaussian
    %cumulative gaussian
    cg = [0 0.1 0.5];           %Gain of cum. Gauss
    cm = [0.05 0.1 0.2];        %Mean of cum. Gauss
    cs = [0.001 0.05 1];        %Std of cum. Gauss
else
    %Fitting the modulation depends a little on the average shape of the
    %modulation. Some monkeys have more oscillatory modulation profiles than
    %others. In principle you could get rid of one of the two Gaussians for less
    %oscillatory profiles.
    
    %The fits depend a lot on the limits, if you have very similar
    %modulation shapes across electrodes then you can use very tight
    %limits and get well behaved fits. The more variability you have the
    %wider the limits will have to be and the greater the chance that you
    %get some weird fits. Check adjusted r-squared values to reject bad fits.
    
    %MODULATION FITTING LIMITS
    %1t Gaussian, mean = 50-125ms [100ms]
    %2nd Gaussian, mean = 100-200ms [150ms]
    %CumGauss, mean = 50-400ms [150ms]
    %First Gaussian
    g1g = [0 0.004 0.2];        %Gain of first Gaussian [min mean max]
    g1m = [0.05 0.1 0.125];     %Mean of first Gaussian
    g1s = [0.005 0.03 0.025];   %STD of first Gaussian
    %second Gaussian
    g2g = [0 0.004 0.3];        %Gain of second Gaussian
    g2m = [0.1 0.15 0.2];       %Mean of second Gaussian
    g2s = [0.005 0.03 0.025];   %STD of second Gaussian
    %cumulative gaussian
    cg = [0 0.06 0.5];          %Gain of second Gaussian
    cm = [0.05 0.15 0.4];       %Mean of second Gaussian
    cs = [0.005 0.05 0.2];      %STD of second Gaussian
end

%Assign to vectors
lower = [g1g(1) g1m(1) g1s(1) g2g(1) g2m(1) g2s(1) cg(1) cm(1) cs(1)];
guess = [g1g(2) g1m(2) g1s(2) g2g(2) g2m(2) g2s(2) cg(2) cm(2) cs(2)];
upper = [g1g(3) g1m(3) g1s(3) g2g(3) g2m(3) g2s(3) cg(3) cm(3) cs(3)];

%Set up the fit
ft = fittype('a.*normpdf(t,b,c) + d.*normpdf(t,e,f) +g.*normcdf(t,h,k);',...
    'independent',{'t'},...
    'coefficients',{'a','b','c','d','e','f','g','h','k'});
fo = fitoptions('method','NonlinearLeastSquares','MaxFunEvals',1000,'Lower',lower,'Upper',upper);
set(fo,'Startpoint',guess);

%Fit the model and get the best fitting parameters
[cfun,gof,output] = fit(t',y',ft,fo);
converge = output.exitflag;
coeff=coeffvalues(cfun);
rs = gof.adjrsquare;

%regenerate best fitting model
yf = coeff(1).*normpdf(t,coeff(2),coeff(3)) + coeff(4).*normpdf(t,coeff(5),coeff(6)) +coeff(7).*normcdf(t,coeff(8),coeff(9));

%Get the latency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%There are (at least) two approaches that can be used to get the latency, either
%searching forwards from trial onset to find the first point at which the
%fitted curve reaches 33% or searching backwards from the maximum of the fitted data to find
%this point. These SHOULD give the same value, however when fitting
%modulation, the two Gaussians sometimes pick up noise components and the
%second approach can be more robust. If the modulation tends to increase over time
%then the second approach often gives too late latencies and the first approach is preferred.
%generally we use the first approach (Searching from trial onset), and if
%it picks up too many noise components we examine the fitting limits to try
%and avoid this.

%calculate the 33% of maximum point
[mf,mfix] = max(yf);

approach = 1; %1 = search forwards from start, 2 = search backwards from maximum

%Find the first model value that comes within 1% (peak) of the 33% point. As we're fitting a smooth curve, this should
%generally work, to be safe we iterate using progressively larger % values
%until we find something.
ixs = [1:mfix]; %only search before maximum model value
perc = 0;
l33ix = [];
while isempty(l33ix) 
    perc = perc+0.01;
    searchval = perc.*mf;
    if approach == 1
        l33ix = find(abs(yf(ixs)-mf*0.33)<searchval,1,'first');
    else
        l33ix = find(abs(yf(ixs)-mf*0.33)<searchval,1,'last');
    end
end
lat = t(l33ix);

%Plot out data, fitted model and latency
if fig
    figure;area(t,y),hold on,plot(t,yf,'r')
    hold on
    plot([lat,lat],get(gca,'YLim'),'m')
end

return

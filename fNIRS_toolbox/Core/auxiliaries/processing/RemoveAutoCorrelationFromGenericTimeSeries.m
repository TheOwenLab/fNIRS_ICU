function [yf,OptimalP] = RemoveAutoCorrelationFromGenericTimeSeries(timeSeries,Pmax)
%Removes autocorrelation from Resting State Data
%with prehitenning methodology
% For details Refer to:
%1 - Characterization and correction of the
%false-discovery rates in resting state connectivity
%using functional near-infrared spectroscopy
%
%2 - Autoregressive model based
%algorithmfor correcting motion and
%seriallycorrelated errors in fNIRS

% Time Series Length
n = size(timeSeries,1);

% Take Original Time Series
y = timeSeries;

for P=1:Pmax
    % For a given parameters P we find the coefficients that
    % minimize autoregressive model (AR(P));
    a = aryule(y,P);
    
    % Once we have the parameters a, we can filter the error
    % to find the new non atucorrelated error (vt).
    vt = filter(a,+1,y);
    
    % Next, we can compute the baysian information
    % criterion (BIC(P)).
    
    % Log Likelihood
    LL = -1*(n/2)*log( 2*pi*mean(vt.^2))+...
        -0.5*(1/mean(vt.^2))*sum(vt.^2);
    
    % Baysian information
    BIC(P) = -2*LL+P*log(n);
end

%Optimal is the P that minimizes BIC
[~,OptimalP] = min(BIC);

AR_Parameters = aryule(y,OptimalP); %Find parameters

% Filter y
yf = filter(AR_Parameters,+1,y);

% Remove undefined points
yf(1:Pmax+1) = [];


end

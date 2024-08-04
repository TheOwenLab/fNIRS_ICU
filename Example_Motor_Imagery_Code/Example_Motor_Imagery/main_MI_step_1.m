%%% Step 1: Analyze MI Data
%%%
%%% In this first step, only the Hemoglobin time-series
%%% are computed
clear

% Available Participants:
% We are just including participants that have at least 
% one good short channel
AvailableParticipants = [1 9 10 12 15 17 18 20 21 22 24 26:31 ...
    33 34 35 37:40]; 

% load data (Controls) -
% Set local folder path
load('Controls_MI');

% Short Channels List
SSlist = [8 29 52 66 75 92 112 125];

% Spline values
IndividualSplineValuesControls_MI;

% Find Channels that do not have acceptable SNR
preproc.SNR_threshold = 8;

% Run for Quality Check for all Subjects
for Nsub = AvailableParticipants
    %%% ******* Step 1: NIRS DATA *********
    
    % Take NIRS data from specific Subject and Run
    obj = dataNIRS{Nsub}{1};
    
    % find bad channels
    Baseline = mean(obj.d);
    obj.d = detrend(obj.d)+Baseline;
    obj = obj.MarkBadChannels(preproc);
    BadChan{Nsub} = obj.SD.BadChannels;
    
    % Not needed for MI
    %BadChan{Nsub} = unique(sort([BadChan{Nsub}]));
    
    clear obj Baseline;
end

% Calculate proper dpf for each participant by taking into
% account their ages
Age_healthy_MI;

% extratc Lambda for each participant
for Nsub = AvailableParticipants    
    
   Wavelegnt_Sub{Nsub} = dataNIRS{Nsub}{1}.SD.Lambda; 
   
end

alpha=223.3; 
beta=0.05624;
gamma=0.8493; 
delta=-5.723*10^-7;
epsilon=0.001245; 
zeta=-0.9025; 

for Nsub = AvailableParticipants    
    
    cnt=0;
    for wave_number= Wavelegnt_Sub{Nsub}
        cnt = cnt+1;
        dpf{Nsub}(cnt)=...
            alpha+beta*(Age(Nsub)^gamma)+delta*(wave_number^3)+...
         epsilon*(wave_number^2)+zeta*wave_number; 
     
    end
    
end



%Run preprocessing for all participants
for Nsub = AvailableParticipants
    
    % Get cw_nirs object
    r = dataNIRS{Nsub}{1};
    
    % Compute Optical Density
    dOD = r.Convert2OD();
    
    % Spline Correction
    dOD_spline = ...
        r.SplineCorrection(dOD,SplineValue{Nsub},[]);
    
    % delete dOD to make sure we are using the right dOD
    clear dOD;
    
    % Wavelet Decomposition in every channel
    r.SD.MeasListAct = ones(size(r.d,2),1);
    Wavelet_iqr = 1.5;
    dOD_wavelet = r.WaveletCorrection(dOD_spline,Wavelet_iqr);
    
    % delete dOD_spline to make sure we are using the right dOD
    clear dOD_spline;
    
    % Estimate concentration
    r.dc = r.Convert2Conc(dOD_wavelet,dpf{Nsub});
    
    % Save data for further analysis
    data{Nsub}{1} = r;
    
    clear r;
    
end

% Save Hemoglobin time series for further analysis
save('processed_data_individual_dpf_MI',...
    'data','BadChan','AvailableParticipants','dataNIRS',...
    'SplineValue','Wavelet_iqr');


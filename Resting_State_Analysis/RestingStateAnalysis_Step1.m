clear

% load data
%load Combined_NIRS_Data_Rest.lob -mat
load('FinalData.lob','-mat')

dataNIRS = data;

% Get List of Spline Values
SplineParameterListOfValues;

% List of Good Short Channels
ListOfGoodShortChannels;

% Avaialable Subjects
AvailableSubjects = [1 9 10 11 12 13 14 15 16 17 ...
    21 24 27 28 29 30 31 33 35 37 38 39 40];  

TempAge =[23 23 23 25 22 28 28 22 28 31 ...
    31 28 28 27 22 22 22 22 20 21 21 20 22];

aCon=[AvailableSubjects;TempAge];
Age(AvailableSubjects)=TempAge; 

%%% calculate proper dpf for each participant 

% extratc Lambda for each participant
for Nsub = AvailableSubjects    
    
   Wavelegnt_Sub{Nsub} = dataNIRS{Nsub}{1}.SD.Lambda; 
   
end

alpha=223.3; 
beta=0.05624;
gamma=0.8493; 
delta=-5.723*10^-7;
epsilon=0.001245; 
zeta=-0.9025; 

for Nsub = AvailableSubjects    
    
    cnt=0;
    for wave_number= Wavelegnt_Sub{Nsub}
        cnt = cnt+1;
        dpf{Nsub}(cnt)=...
            alpha+beta*(Age(Nsub)^gamma)+delta*(wave_number^3)+...
         epsilon*(wave_number^2)+zeta*wave_number; 
     
    end
    
end


% Check bad Channels
% We will not exclude bad channels now,
% but we have to be carefull when calculating the correlations

% Short Channels List
SSlist = [8 29 52 66 75 92 112 125];
preproc.SNR_threshold=8;

for Nsub = AvailableSubjects
    
    % Find Bad Long-Channels for each volunteer
    obj_nirs =  dataNIRS{Nsub}{1};
    Baseline = mean(obj_nirs.d);
    obj_nirs.d = detrend(obj_nirs.d)+Baseline;
    obj_nirs = obj_nirs.MarkBadChannels(preproc);
    
    % Save it on a cell for further analysis
    BadChan{Nsub} = obj_nirs.SD.BadChannels;
    
    clear obj_nirs Baseline;
end

BP_f = [0.009 0.08];


% Run Analysis for all volunteers
for Nsub = AvailableSubjects
    
    % Get cw_nirs Data
    r = dataNIRS{Nsub}{1};
    
    % We will work only with the first 6 minutes
    % Therefore, we have the remove the data
    % after 6 minutes. 
    
    Frame_1_min = round(60*r.SD.f);
    Frame_6_min = round(360*r.SD.f);
    
if Nsub<=17
        
    r.t(1:Frame_1_min+1) = [];
    r.d(1:Frame_1_min+1,:) = [];
    r.s(1:Frame_1_min+1) = [];
    r.aux(1:Frame_1_min+1,:) = [];
    
    r.t(Frame_6_min+1:end) = [];
    r.d(Frame_6_min+1:end,:) = [];
    r.s(Frame_6_min+1:end) = [];
    r.aux(Frame_6_min+1:end,:) = [];
    
elseif Nsub>17
    
    r.t(Frame_6_min+1:end) = [];
    r.d(Frame_6_min+1:end,:) = [];
    r.s(Frame_6_min+1:end) = [];
    r.aux(Frame_6_min+1:end,:) = [];
end 
    
    % Compute Optical density
    % Compute dOD
    dOD = r.Convert2OD();
    
    % Spline Correction
    SplineValue = SplineParameter{Nsub}{1};
    dOD_spline = r.SplineCorrection(dOD,SplineValue,[]);
    
    % Wavelet Decomposition in every channel
    r.SD.MeasListAct = ones(size(r.d,2),1);
    Wavelet_iqr = 1.5;
    dOD_wavelet = r.WaveletCorrection(dOD_spline,Wavelet_iqr);
    
    % Estimate concentration
    r.dc = r.Convert2Conc(dOD_wavelet,dpf{Nsub});
    
    % Perform Band-Pass Filtering
    r.dc(:,:,1) = detrend(r.dc(:,:,1));
    r.dc(:,:,2) = detrend(r.dc(:,:,2));
    r.dc(:,:,3) = detrend(r.dc(:,:,3));
    
    r.dc = hmrBandpassFilt( r.dc, r.SD.f, BP_f(1), BP_f(2));
    
    % Remove Border Effects
    r.dc(1:100,:,:) = [];
    
    % Adjust Temporal Information
    r.t = linspace(0,size(r.dc,1)/r.SD.f,size(r.dc,1))';
    
    % Save Data for Further Analysis
    NIRS_DATA{Nsub} = r;
    
    clear r;
    
end

%%
% Perform Regression

% GLM Parameters.
% Perform GLM with all available SC
opt_GLM = 4;

% Use all components
CV = 100;

% Temporal shift is not allowed
flag_time_SC = 0;
flag_time_AD = [];

for Nsub = AvailableSubjects
    
    % Take Data from specific Volunteer
    % raw_data means not regressed in this context
    raw_data = NIRS_DATA{Nsub};
        
    % Perform the Robust Regression with the SC regression and also
    % the physiology regression
    [filtered_dc_SC,Stats_SC] = ...
        raw_data.PerformPhysiologyRegression...
        (SSlist_good{Nsub},[],CV,...
        opt_GLM,flag_time_SC,flag_time_AD);
        
    % Remove autocorrelation from regressed
    % data    
    % prewhiten regressed time series
    r_pw_regressed = raw_data;
    r_pw_regressed.dc = filtered_dc_SC;
    
    r_pw_regressed = r_pw_regressed.RemoveAutoCorrelation();
    r_pw_regressed.dc = r_pw_regressed.dc(51:end,:,:);
        
    % Compute Correlation Matrix of HbO and HbR for each data
    % Regressed - HbO
    C_regressed(:,:,1) = corrcoef(r_pw_regressed.dc(:,:,1));
    % Regressed - HbR
    C_regressed(:,:,2) = corrcoef(r_pw_regressed.dc(:,:,2));
    % Regressed - HbT
    C_regressed(:,:,3) = corrcoef(...
        r_pw_regressed.dc(:,:,1) + r_pw_regressed.dc(:,:,2));
           
    % Save Data for further analysis.
    SC_Regression{Nsub}{1}.stats = Stats_SC;
    SC_Regression{Nsub}{1}.regressed_dc = filtered_dc_SC;
    SC_Regression{Nsub}{1}.regressed_dc_pw = r_pw_regressed.dc;
    SC_Regression{Nsub}{1}.CV = CV;
    SC_Regression{Nsub}{1}.C_regressed = C_regressed;
    
    
    clear r* C_regressed;
    
end

save('Regression_Physiology_6_minutes.stats',...
       'SC_Regression','BadChan');
 
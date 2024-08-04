% Code to perform first level GLM analysis
clear

% load processed data
load('processed_data_individual_dpf_MI.mat');

clear dataNIRS;

% Short Channels List
SSlist = [8 29 52 66 75 92 112 125];
GoodSC_Controls_MI;

% GLM option and Stim duration
opt_GLM = 4;

duration{1} = 30;

% Run analysis for all participants
for Nsub = AvailableParticipants
    
    r =  data{Nsub}{1};
    
    % Low-pass filter the data
    r.dc = r.BPFilter(r.dc, [0.005 0.5]);

    % detrend the data
    r.dc(:,:,1) = detrend(r.dc(:,:,1));
    r.dc(:,:,2) = detrend(r.dc(:,:,2));
    r.dc(:,:,3) = detrend(r.dc(:,:,3));
    
    % Remove data before the first "Rest" trigger
    lst_min = min(find(r.s(:,2)==1)); % In MI, rest is the second trigger
    
    r.s(1:lst_min,:) = [];
    r.dc(1:lst_min,:,:) = [];
    r.aux(1:lst_min,:) = [];
    r.t(1:lst_min,:) = [];
    
    %Remove data after the last "Rest" + 30 seconds
    lst_max = max(find(r.s(:,2)==1))+30*round(r.SD.f);

    r.s(lst_max:end,:) = [];
    r.dc(lst_max:end,:,:) = [];
    r.aux(lst_max:end,:) = [];
    r.t(lst_max:end,:) = [];
    
    clear lst_max lst_min;
 
    
    % Make time vector to start at zero
    r.t = r.t-r.t(1);
    
    % Remove the second type of trigger (rest)
    r.s(:,2) = [];
        
    % Remove bad channels from hemoglobin time series
    % before perfoming the GLM analysis to infer activated 
    %channels
    r.dc(:,BadChan{Nsub},:) = nan;
        
    % Perform GLM with the time series without the SC regression because 
    % the SC regression will be performed as a single step with the GLM.
    [beta_R{Nsub},p_R{Nsub},covb_R{Nsub}] = r.GLM_Prewhitening...
        (duration,GoodSC{Nsub},[],[],[],[],opt_GLM);
    
    % Update data
    data{Nsub}{1} = r;
    
end

save('GLM_processed_MI','beta_R','p_R',...
    'covb_R','BadChan','data',...
    'AvailableParticipants','SSlist');








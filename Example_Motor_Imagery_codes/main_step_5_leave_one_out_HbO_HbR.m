% Function to compute Sensitivity with the leave-one-out
% approach.

clear

% load GLM stats: 1st level analysis
load GLM_processed_MI.mat

% Create Vector of Contrast:
% The vector of contrast is huge to mimic the case in which
% all available regressors were good.
% But only the first entries are important.

C = zeros(1,20);

C(1) = 1;
C(2) = 0;
C(3) = 0;

% List with all Short-channels
SSlist = [8 29 52 66 75 92 112 125];

% Get Stats from the First Level for all participants
[B_k,Cov_k] = ...
    ExtractDataFromFirstLeveL...
    (beta_R,covb_R,C,SSlist,AvailableParticipants,...
    BadChan);

% Compute T-Values for All Subjects
T = B_k./sqrt(Cov_k);

% Perform intra-subject analysis for the defined contrast.
% This steps will evaluate which channels were activated for
% each participants. This result will used to compute the
% sensitivity results

cnt_sub = 0;
at_least_one = [];
at_least_one_sub = [];

for Nsub = AvailableParticipants
    
    % Counter
    cnt_sub = cnt_sub +1;
    
    % Get cw-nirs object to infer the degress of freedom
    r = data{Nsub}{1};
    
    % Compute degrees of freedom as the length of the time-series
    DegreeOfFreedom = size(r.dc,1) - 50;
    
    % Convert T to P
    % HbO
    p_valueHbO = 1-tcdf(abs(T(cnt_sub,:,1)),DegreeOfFreedom);
    
    % HbR
    p_valueHbR = 1-tcdf(abs(T(cnt_sub,:,2)),DegreeOfFreedom);
    
    
    % Compute Act: HbO AND HbR
    Act{cnt_sub} = find(p_valueHbO<0.05 & p_valueHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
    
    % Compute Act: ONLY HbO 
    Act_HbO{cnt_sub} = find(p_valueHbO<0.05 & ...
        T(cnt_sub,:,1)>0);
    
    % Compute Act: ONLY HbR
    Act_HbR{cnt_sub} = find(p_valueHbR<0.05 & ...
        T(cnt_sub,:,2)<0);
    
    
    % Save activated channels following the original index
    Act_2{Nsub} = find(p_valueHbO<0.05 & p_valueHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
       
    
    % Check which participants had at least one activated channel
    if ~isempty(Act{cnt_sub})
        at_least_one = [at_least_one,cnt_sub];
        at_least_one_sub = [at_least_one_sub,Nsub];
    end
    
    clear p_valueHbO p_valueHbR R DegreeOfFreedom
end


% Group Analysis will be performed with only participants
% that had at least one activated channel
% In addition, we will also remove the left out
% participant

sensitivity_bin = zeros(size(AvailableParticipants,2),1);
sensitivity_HbO = zeros(size(AvailableParticipants,2),1);
sensitivity_HbR = zeros(size(AvailableParticipants,2),1);


% Left out participant
for N_sub_out = 1:size(AvailableParticipants,2)
    
    
    % Remove the left out participant from "at_least one"
    
    at_least_one_minus_left_out = at_least_one;
    
    lst_left_out = ...
        find(at_least_one_minus_left_out==N_sub_out);
    
    at_least_one_minus_left_out(lst_left_out) = [];
    
    
    % Perform Group Analysis
    B_k_new = B_k(at_least_one_minus_left_out,:,:);
    Cov_k_new = Cov_k(at_least_one_minus_left_out,:,:);
    
    % Perform second-level analysis for HbO and HbR
    
    for Hb=1:2
        
        [beta_group(:,Hb),p_group(:,Hb),~,t_value(:,Hb)] = ...
            WeightLinearGroupAnalysis ...
            (B_k_new(:,:,Hb)',SSlist,Cov_k_new(:,:,Hb)');
        
    end
    
    
    % Find Activated Channels for the Group
    lst_Act_group = find(beta_group(:,1)>0 & beta_group(:,2)<0 ...
        & p_group(:,1)<0.05 & p_group(:,2)<0.05);
    
    save_act_channels{N_sub_out} = lst_Act_group;
    
    
    clear T r B_k_new Nchan Nsub ...
        p* beta* cnt_sub Cov_k_new Hb covb_R DegreeOfFreedom
    
    
    % Sensitivity for binary condition (i.e., HbO AND HbR)
    if ~isempty(intersect(Act{N_sub_out},lst_Act_group))
        
        sensitivity_bin(N_sub_out) = 1;
        
    end
    
    % Sensitivity for ONLY HbO
    if ~isempty(intersect(Act_HbO{N_sub_out},lst_Act_group))
        
        sensitivity_HbO(N_sub_out) = 1;
        
    end
    
    % Sensitivity for ONLY HbR
    if ~isempty(intersect(Act_HbR{N_sub_out},lst_Act_group))
        
        sensitivity_HbR(N_sub_out) = 1;
        
    end
    
end

sensitivity = ...
    sum(sensitivity_bin)/size(AvailableParticipants,2)

sensitivity_HbO = ...
    sum(sensitivity_HbO)/size(AvailableParticipants,2)

sensitivity_HbR = ...
    sum(sensitivity_HbR)/size(AvailableParticipants,2)


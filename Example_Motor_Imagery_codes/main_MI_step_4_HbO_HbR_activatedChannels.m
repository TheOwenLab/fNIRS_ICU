% Code to extract: 

% Identify activated channels for HbO and HbR. 
%
% It requires the ROI generated in step_3 and GLM results generated in
% step_2. 

clear

% load GLM stats: 1st level analysis
load GLM_processed_MI.mat

% Create Vector of Contrast:
% The vector of contrast is huge to mimic the case in which
% all available regressors were good.
% But only the first entries are important.

C = zeros(1,17);

C(1) = 1;

% List with all Short-channels
SSlist = [8 29 52 66 75 92 112 125];

% Get Stats from the First Level
[B_k,Cov_k] = ...
    ExtractDataFromFirstLeveL...
    (beta_R,covb_R,C,SSlist,AvailableParticipants,...
    BadChan);

% Compute T-Values for All Subjects
T = B_k./sqrt(Cov_k);

% Perform intra-subject analysis for the defined contrast.
% This steps will evaluate which channels were activated for
% each participant. This result will used to compute the
% sensitivity results

cnt_sub = 0;
at_least_one = [];
at_least_one_HbO = [];
at_least_one_HbR = [];

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
    ActHbO{cnt_sub} = find(p_valueHbO<0.05 & ...
        T(cnt_sub,:,1)>0);
    
    
    p_valueHbR = 1-tcdf(abs(T(cnt_sub,:,2)),DegreeOfFreedom);
    ActHbR{cnt_sub} = find(p_valueHbR<0.05 & ...
        T(cnt_sub,:,2)<0);
    
    % Combine HbO and HbR
    Act{cnt_sub} = find(p_valueHbO<0.05 & p_valueHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
    % Check which participants had at least one activated channel
    if ~isempty(Act{cnt_sub})
        at_least_one = [at_least_one,cnt_sub];
    end
    
    if ~isempty(ActHbO{cnt_sub})
        at_least_one_HbO = [at_least_one_HbO,cnt_sub];
    end
    
    if ~isempty(Act{cnt_sub})
        at_least_one_HbR = [at_least_one_HbR,cnt_sub];
    end
    
    clear p_valueHbO p_valueHbR R DegreeOfFreedom
end

% Compute ROI for each chromophore

% 1 HbO - Considering only participants that had HbO response
B_k_new_HbO = B_k(at_least_one_HbO,:,:);
Cov_k_new_HbO = Cov_k(at_least_one_HbO,:,:);

% Perform second-level analysis for HbO 
[beta_group_HbO(:,1),p_group_HbO(:,1),~,t_value_HbO(:,1)] = ...
    WeightLinearGroupAnalysis ...
    (B_k_new_HbO(:,:,1)',SSlist,Cov_k_new_HbO(:,:,1)');

lstActHbO = find(beta_group_HbO(:,1)>0 & p_group_HbO(:,1)<0.05);

% 2 HbR - Considering only participants that had HbO response

B_k_new_HbR = B_k(at_least_one_HbR,:,:);
Cov_k_new_HbR = Cov_k(at_least_one_HbR,:,:);

% Perform second-level analysis for HbO 
[beta_group_HbR(:,1),p_group_HbR(:,1),~,t_value_HbR(:,1)] = ...
    WeightLinearGroupAnalysis ...
    (B_k_new_HbR(:,:,2)',SSlist,Cov_k_new_HbR(:,:,2)');

lstActHbR = find(beta_group_HbR(:,1)<0 & p_group_HbR(:,1)<0.05);

% Code to perform analyse GLM stats
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
% This step will evaluate which channels were activated for
% each participant. This result will be used to compute the
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
    DegreeOfFreedom = size(r.dc,1) -50; 
    
    % Convert T to P
    % HbO
    p_valueHbO = 1-tcdf(abs(T(cnt_sub,:,1)),DegreeOfFreedom);
    
    % HbR
    p_valueHbR = 1-tcdf(abs(T(cnt_sub,:,2)),DegreeOfFreedom);
    
    Act{cnt_sub} = find(p_valueHbO<0.05 & p_valueHbR<0.05 & ...
        T(cnt_sub,:,1)>0 & T(cnt_sub,:,2)<0);
    
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

B_k_new = B_k(at_least_one,:,:);
Cov_k_new = Cov_k(at_least_one,:,:);

% Export Degree of freedom to excel
%Compute_degree_freedom_group_analysis(B_k_new,C);


% Perform second-level analysis for HbO and HbR
for Nchan = 1:size(B_k_new,2)
    
    for Hb=1:2
        
        [beta_group(:,Hb),p_group(:,Hb),~,t_value(:,Hb)] = ...
            WeightLinearGroupAnalysis ...
            (B_k_new(:,:,Hb)',SSlist,Cov_k_new(:,:,Hb)');
        
    end
    
end

% Find activated channels for the group
lst_Act_group = find(beta_group(:,1)>0 & beta_group(:,2)<0 ...
    & p_group(:,1)<0.05 & p_group(:,2)<0.05);

% Find Activated Channels for HbO
lst_Act_group_HbO = find(beta_group(:,1)>0 ...
    & p_group(:,1)<0.05);

% Find Activated Channels for HbR
lst_Act_group_HbR = find(beta_group(:,2)<0 ...
    & p_group(:,2)<0.05);

clear T r B_k B_k_new Nchan Nsub ...
    p* beta* cnt_sub Cov* Hb covb_R DegreeOfFreedom


% Save ROI/Group info for further analysis
save('ROI_group_health_MI','lst_Act_group',...
    'at_least_one_sub','t_value');




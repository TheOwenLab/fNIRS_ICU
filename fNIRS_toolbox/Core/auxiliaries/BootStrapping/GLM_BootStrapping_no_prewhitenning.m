function [betaChan,pValue,covb,Z_score] = GLM_BootStrapping_no_prewhitenning...
    (r,duration,lstSS,SC_exclude,HbO_regressor, HbR_regressor,...
    AdditionalRegressors,opt)
%%%  GLM_WithPrewhiteningRobustFit performs General Linear regression
%%% to recover coefficients that represent the hemodynamic response due
%%% to a task. The implemented method remove the autocorrelation before
%%% performing the regression, and the regression is performed with the
%%% robust fit of MatLab in order to avoid outliers.
%%%
%%% *******************
%%%       INPUT
%%% *********************
%%% r: cw_nirs object.
%%% duration: stim duration is now a cel in which duration{Ntask} has
%%% the stim durations for each trigger from the Ntask collum inside the
%%% stim vector (r.s). If duration{Ntask} is a constant, the code assumes
%%% that each trial had the same duration. If duration{Ntask} is a vector
%%% then it should match the number of trials. In addition, note that the
%%% length of duration (lenght(duration)) should match the number of
%%% colulums from the stim vector (r.s) (i.e, the user has to specify the
%%% duration of each task to properly create a design matrix to perform the
%%% GLM analysis.
%%% lstSS: List of Short channels
%%% SC_exclude: List of Short channel to be excluded
%%% HbO_regressor: Time series to regress HbO
%%% HbR_regressor: Time Series to regress HbR
%%%
%%% opt: 0 - Do not Add SC in the GLM model;
%%%      1 - Add the closest SC for each channel;
%%%      2 - Add every SC in the GLM model;
%%%      3 - Input data to use as regressor
%%%      4 - Add every SC in the GLM model (HbO+HbR) and Perform PCA
%%%      in the SC data to remove collinearity across regressors;
%%%
%%% *******************
%%%       OUTPUT
%%% *********************
%%% betaChan is the Beta Value of each channel that corresponde to the hrf
%%% p-value is the siginificance from a ONE-SIDED t-TEST
%%% covb is the variance for the first regressor.
%%% Stats_nChan is the full regression stats from each channel and Hb.
%%% Examples:
%%%
%%% [beta,value] = GLM_...(r, duration, lstSS, [], [], [] ,1) will perform
%%% the GLM with the closest short channel. If you want to remove a bad SC
%%% you can input on SC_exlude. However, the variable lstSS should always
%%% contain all Short channels. For example, if your short channels are the
%%% channels [1 2 3 4] but you want to remove SC number 4. You should do:
%%% [beta,value] = GLM_...(r, duration, [1 2 3 4],4, [], [] ,1).
%%%
%%%[beta,value] = GLM_...(r, duration, lstSS, [], [], [] ,2) will perform
%%% the regression with all short channels.
%%%
%%% [beta,value] = GLM_...(r, duration, lstSS, [], HbO_regressor, ...
%%% HbR_regressor ,3) will perform the regression with only the specified
%%% time series in HbO_regressor and HbR_regressor.
%%%
%%%[beta,value] = GLM_...(r, duration, lstSS, [], [], [] ,0) will NOT
%%% perform the Short channel regression. The design matrix will only have
%%% the beta model for the hrf.
%%%
%%%
%%%
%%% Sergio Luiz Novi Junior, sergiolnovi@gmail.com
%%%
%%%

% Save Original list of Short Channels
SSlist = lstSS;

if opt==1
    % Let us find the closest SC per channel
    nch = size(r.SD.MeasList,1)/2;
    ml=r.SD.MeasList;
    rhoSD = zeros(nch,1);
    
    % Verify which channels are SC
    if ~isempty(SC_exclude)
        for Nchan=1:length(SC_exclude)
            exclude = find(lstSS==SC_exclude(Nchan));
            lstSS(exclude)=[];
            clear exclude;
        end
    end
    
    % Find Nearest Short channels for the remaning good SC
    [channelspershortch ] = Nearshortchanel_GLM(lstSS,r.SD);
    % Lets flip the cell and create a Lista thar contains
    % the number of the SC for each channel.
    Lista = zeros(nch,1);
    for Nsc=1:size(channelspershortch,2)
        indexes = channelspershortch{Nsc};
        for dummy=1:size(indexes,2)
            Lista(indexes(dummy)) = lstSS(Nsc);
        end
    end
    
end

% Lets not perform the GLM in NAN channels or in the SC
% It is better to avoid the original list of SC

Avoid = SSlist;
%%% Check for channels with NAN
for Nchan=1:size(r.dc,2)
    if ~isempty(find(isnan(r.dc(:,Nchan,1:2))==1))
        Avoid =[Avoid,Nchan];
    end
end

for Hb=1:2
    
    
    for Nchan=1:size(r.dc,2)
        
        
        % Design Matrix for every stimuli
        X=[];
        for Nstim=1:size(r.s,2)
            X = [X,CreateDeasingMatrix(r,r.s(:,Nstim),duration{Nstim})];
        end
        
        if opt==1
            % Add only one SC to the design matrix
            X = [X,r.dc(:,Lista(Nchan),Hb)];
            
        elseif opt==2
            % Add all SCs to the design matrix
            X = [X,r.dc(:,lstSS,Hb)];
        elseif opt==3
            if Hb==1
                % Add external regressor
                X = [X,HbO_regressor];
            else
                % Add external regressor
                X = [X,HbR_regressor];
            end
            
        elseif opt==4
            
            % Add HbO and HbR
            Xshort = [r.dc(:,lstSS,1),r.dc(:,lstSS,2)];
            
            % ** Principal Component Analysis (PCA) **
            % Subtract Mean
            Xshort = Xshort-mean(Xshort);
            
            % Covariance matrix
            covar = cov(Xshort);
            
            % Singular Value Decomposition
            [u,s,v] = svd(covar);
            
            % Rotate SC Data (New Basis of PCA)
            Xshort_r = Xshort*v;
            %Xshort_r = Xshort_r(:,1:2);
            % Check if there additional regressors
            if ~isempty(AdditionalRegressors)
                
                % Perform Shift in the AdditionalRegressors so that it
                % increases corrrelation with the hemoglobin-time series
                % from the long channel.
                
                for NAdditional = 1:size(AdditionalRegressors,2)
                    % Take one signal at a time
                    Regressor = AdditionalRegressors(:,NAdditional);
                    
                    % Perform the Shifting
                    ShiftedAdditionalRegressors(:,NAdditional) = ...
                        PerformRegressorTaskProtolForGLM...
                        (r.dc(:,Nchan,Hb),Regressor,round(20*r.SD.f));
                    
                end
                
                Xshort_r = [Xshort_r,ShiftedAdditionalRegressors];
                
            end
            
            % Update Design Matrix
            X = [X,Xshort_r];
            
        end
        
        if isempty(find(Avoid==Nchan))
            
            y = r.dc(:,Nchan,Hb);
            n = size(y,1); % Number of points of the Time series
            
            % Subtract the mean of y
            y = y - mean(y);
                   
            % Solve the GLM with Robustfit for the original 
            % time series
            [beta, Stats] = ...
                robustfit(X,y,'bisquare',[],'off');
                                    
            % Save beta and p-value from original data
            % with the standard procedure
            betaChan(Nchan,Hb,:) = beta;
            pValue(Nchan,Hb,:) = Stats.p/2;
            covb{Nchan}{Hb} = Stats.covb;
            T(Nchan,Hb,:) = Stats.t;
            
            
            % BootStrapping approach
            for Nboot = 1:500
                
                % Scramble time series with the phase method
                y_scrambled = ...
                    add_phase_shift_2_timeSeries(y,[0 2*pi]);
                
                % subtract the mean of the time series
                y_scrambled = y_scrambled - mean(y_scrambled);
                
                % Solve the GLM again
                [beta_s, Stats_s] = ...
                robustfit(X,y_scrambled,'bisquare',[],'off');
            
                % Save T-value from the scrambled version
                T_scrambled(Nboot,Hb,:) = Stats_s.t;
                
            end
            
            % Compute Z-score
            mean_boot = squeeze(mean(T_scrambled(:,Hb,:)));
            std_boot = squeeze(std(T_scrambled(:,Hb,:)));
            
            Z_score(Nchan,Hb,:) = ...
               (squeeze(T(Nchan,Hb,:)) - mean_boot)./std_boot;                       
            
        end
        
        
    end
end

betaChan(Avoid,:,:)=nan;
pValue(Avoid,:,:)=nan;
%covb(Avoid,:,:)=nan;


end



function[filtered_dc,Stats] = ...
    PhysiologyRegression_GLM(obj,SSlist,AdditionalRegressors,CV,opt,...
    flag_time_SC,flag_time_AD)

% This function performs physiological Regression based
% on Short-channel and with additional Regressors of physiology.
% The additional Regressors must be inputed.
% The regressions will be performed with the robustfit matlab function:
% robustfit(X,y,'bisquare',[],'on');
%
%   INPUT:
%       obj: cw_nirs object.
%
%       AdditionalRegressors: Physiological data to be regressed in
% addition to the short channels. The dimension of this data should match
% the ones from the concentration data (obj.dc).
%
%       opt: 0 - Perform Short-channel regression with the closest SC.
%            1 - Include all Short channels to perform regression.
%            2 - Use the SC with highest correlation with the Long Channel.
%the correlation is performed on the HBO time series.
%            3 - Do not perform SC regression. Use only the additional
% regressor.
%            4 - Combined SC data from HbO and HbR with PCA. The PCA is
%            used to remove collinearity.
%
%
%   OUTPUT:
%       filtered_dc: This is the filtered concentration data after all
% regression, which is basically the residual of the robustfit.


% Define the chosen SC for each long channel based on the distance
if opt<3
    SC_per_LongChannel = CreateListOfShortChannelPerLongChannel_GLM...
        (obj,SSlist,opt);
end

for Hb=1:2
    
    for Nchan=1:size(obj.dc,2)
        
        % Step 1: Perform Regression in a given channel Nchan
        y = obj.dc(:,Nchan,Hb);
        
        % Step 2: Create Design Matrix (X) for regression
        
        X=[];
        
        % 1.1 Short Channel Data
        if opt == 0 || opt == 2
            %Using only the closest Short channel
            Nshort = SC_per_LongChannel(Nchan);
            X = obj.dc(:,Nshort,Hb);
        end
        
        if opt == 1
            % Using all short channels
            X = obj.dc(:,SSlist,Hb);
            % Add HbO and HbR
            %X = [obj.dc(:,SSlist,1),obj.dc(:,SSlist,2)];
            
        end
        
        if opt == 3
            X=[];
        end
        
        % I will separate option 4 in to steps.
        
        if opt == 4
            
            % Add HbO and HbR
            Xshort = [obj.dc(:,SSlist,1),obj.dc(:,SSlist,2)];
            
            
            % Now, I will perform the temporal shifts in the
            %X_short design matrix.
            
            if flag_time_SC == 1
                % I will allow a max lag of 15 seconds (for SC), which
                % could ultimately corresponds to one full cicle of
                % sistemic physiology, such as Mayer Waves
                
                maxLag = round(15*obj.SD.f);
                [y,Xshort,shift_SC,coor_max] = ...
                    AdjustTemporalShift_for_Regression(y,Xshort,maxLag);
                
            end
            
            % ** Principal Component Analysis (PCA) **
            % Subtract Mean
            Xshort = Xshort-mean(Xshort);
            
            % Covariance matrix
            covar = cov(Xshort);
            
            % Singular Value Decomposition
            [u,s,v] = svd(covar);
            
            % Define Number of components
            CovarImportance = 100*diag(s)/sum(diag(s));
            
            if CV<100
                CovarSum = 0;
                Ncomponents=0;
                
                while CovarSum < CV
                    Ncomponents=Ncomponents+1;
                    CovarSum = CovarSum+CovarImportance(Ncomponents);
                end
            else
                Ncomponents = size(CovarImportance,1);
            end
            
            % Rotate SC Data (New Basis of PCA)
            Xshort_r = Xshort*v(:,1:Ncomponents);
            
            % Update Design Matrix
            X = [X,Xshort_r];
            
        end
        
        if ~isempty(AdditionalRegressors)
            % 1.2 Add Additional Regressors
            % if  flag_time_Ad == 0, we are not allowing shifts in the
            % additional regressors.
            if flag_time_AD == 0
                
                X = [X,AdditionalRegressors];
                
                % PCA to remove collinearity
                X = X-mean(X);
                
                % Covariance matrix
                covar = cov(X);
                
                % Singular Value Decomposition
                [u,s,v] = svd(covar);
                
                % Rotate X
                X = X*v;
                
                % if flag_time_Ad == 1, for each channel, we estimate the shift
                % for a maximum of 20 seconds based on Tong's works.
                % The shift is circular and is done for each Hb and each
                % phisiology.
                
            elseif flag_time_AD == 1
                
                % I will allow a lag max of 20 seconds
                % for the Additional Regressor.
                
                maxLag = round(20*obj.SD.f);
                [y,AdditionalRegressors_s,shift_AD,coor_max] = ...
                    AdjustTemporalShift_for_Regression...
                    (y,AdditionalRegressors,maxLag);
                
                % We have to cut the begning and end of the
                % design matrix as we did with the shifted
                % time series and y vector.
                if ~isempty(X)
                    X(1:maxLag,:) = [];
                    X(end-maxLag:end,:) = [];
                end
                % Next, we add the shifted additional regressors
                % to the design matrix
                
                X = [X,AdditionalRegressors_s];
                
            end
        end
        
        % Perform Robust Fit Regression
        [Dummy, StatsDummy] = robustfit(X,y,'bisquare',[],'on');
        
        if opt==4
            StatsDummy.Ncomponents = Ncomponents;
            StatsDummy.GLM_opt = opt;
        end
        
        % Save filtered data (residual)
        filtered_dc(:,Nchan,Hb) = StatsDummy.resid;
        
        % Save Short-Channels Shifts for further analysis
        if exist('shift_SC')
            StatsDummy.shift_SC = shift_SC./obj.SD.f;
            StatsDummy.coor_SC = coor_max;
        end
        
        % Save Additional Regressors Shifts for further analysis
        if exist('shift_AD')
            StatsDummy.shift_AD = shift_AD./obj.SD.f;
            StatsDummy.coor_AD = coor_max;
        end
        
        % Save Stats for further analysis
        Stats{Nchan}{Hb} = StatsDummy;
        
        clear y X StatsDummy;
        
        
        
    end
    
    
    
    
    
end






end




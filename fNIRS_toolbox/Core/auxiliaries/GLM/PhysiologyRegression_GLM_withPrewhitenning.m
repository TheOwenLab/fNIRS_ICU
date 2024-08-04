function[filtered_dc,Stats] = ...
    PhysiologyRegression_GLM_withPrewhitenning(obj,SSlist,AdditionalRegressors,CV,opt)

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
%            3 - Do not perform SC regression. Use only the additional
% regressor.
%            4 - Combined SC data from HbO and HbR with PCA.
%the correlation is performed on the HBO time series.
%
%   OUTPUT:
%       filtered_dc: This is the filtered concentration data after all
% regression, which is basically the residual of the robustfit.


% Define the chosen SC for each long channel based on the distance
if opt<3
    SC_per_LongChannel = CreateListOfShortChannelPerLongChannel_GLM...
        (obj,SSlist,opt);
end

% Define Pmax;
Pmax = round(10*obj.SD.f);


for Hb=1:2
    
    for Nchan=1:size(obj.dc,2)
        
        % Step 1: Create Design Matrix (X) for regression
        
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
        end
        
        if opt == 3
            X=[];
        end
        
        if opt == 4
            
            % Add HbO and HbR
            Xshort = [obj.dc(:,SSlist,1),obj.dc(:,SSlist,2)];
            
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
        
        
        
        % 1.2 Add Additional Regressors
        X = [X,AdditionalRegressors];
        
        % Step 2: Perform Regression in a given channel Nchan
        y = obj.dc(:,Nchan,Hb);
        
        n = size(y,1); % Number of points of the Time series
        change = 1;
        
        % Estimate beta with OLS
        betaOld = inv(X'*X)*X'*y;
        
        % Estimate error
        e = y-X*betaOld;
        
        Pmax = min([Pmax length(y)]);
        maxIte = 0;
        
        
        while change>=0.01 & maxIte<10
            
            for P=1:Pmax
                
                % For a given parameters P we find the coefficients that
                % minimize autoregressive model (AR(P));
                a = aryule(e,P);
                
                % Once we have the parameters a, we can filter the error
                % to find the new non atucorrelated error (vt).
                vt = filter(a,+1,e);
                
                % Next, we can compute the baysian information
                % criterion (BIC(P)).
                
                % Log Likelihood
                LL = -1*(n/2)*log( 2*pi*mean(vt.^2))+...
                    -0.5*(1/mean(vt.^2))*sum(vt.^2);
                
                % Baysian information
                BIC(P) = -2*LL+P*log(n);
                
                
            end
            
            maxIte = maxIte+1;
            % The order (P) of the autoregressive model AR(P) is the value
            % tha minmize BIC.
            [~,OptimalP] = min(BIC);
            
            a = aryule(e,OptimalP); %Find parameters
            
            % Filter y
            yf = filter(a,+1,y);
            
            % Filter X
            Xf = filter(a,+1,X);
            
            % Compute new beta and Stats of beta
            [beta, StatsDummy] = robustfit(Xf,yf,'bisquare',[],'off');
            
            % Update error
            % BE CAREFUL: This error is not the one for the final fit
            % The error of the final fit, i.e. for the prewithened
            % equation is inside Stats. Note Below that Stats.covb
            % is saved in covb.
            % If you want to take a look the error residual, you have
            % to check Stats.
            %
            %
            e = y-X*beta;
            
            % Compute Stopping Criterion
            change = norm(beta-betaOld)./norm((betaOld));
            
            % Update betaOld
            betaOld = beta;
            
        end
        
        if opt==4
            StatsDummy.Ncomponents = Ncomponents;
            StatsDummy.GLM_opt = opt;
        end
    
    % Save filtered data (residual)
    filtered_dc(:,Nchan,Hb) = StatsDummy.resid;
    
    % Save Stats for further analysis
    Stats{Nchan}{Hb} = StatsDummy;
    
    clear y X StatsDummy;
        
    end
    
    
    
    %         % Perform Robust Fit Regression
    %         [Dummy, StatsDummy] = robustfit(X,y,'bisquare',[],'on');
    %
    %         if opt==4
    %             StatsDummy.Ncomponents = Ncomponents;
    %             StatsDummy.GLM_opt = opt;
    %         end
    %
    %         % Save filtered data (residual)
    %         filtered_dc(:,Nchan,Hb) = StatsDummy.resid;
    %
    %         % Save Stats for further analysis
    %         Stats{Nchan}{Hb} = StatsDummy;
    
    
    
    
    
end




end




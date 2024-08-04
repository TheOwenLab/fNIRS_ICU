function[betaAverage,p,VarChan,tvalue] = WeightLinearGroupAnalysis...
    (betaGroup,AvoidChan,covb)

%Second Level Analysis: Weighted Linear Regression
% *** INPUT:
%   1 - betaGroup is a n by m matrix with betas from HbO or HbR for the
% whole group of volunteers. n is the number of channels and m is the
% number of volunteers.
%   2 - AvoidChan is a vector with the list of channels to be avoided.
%If there is no channel to be avoided, enter an empty vector ([]);
%   3 - covb is a n by m matrix with the variance of each channel of each
%volunteer. n is the number of channels and m is the number of
%volunteer. The variances come from first level analysis.
%
%*** OUTPUT
%   1 - betaAverage - Vector with the average beta for each channel.
%   2 - p - p-value of each channel computed with one-sided t-test.
%   3 - VarChan - Variance of each average beta of each channel
%
%
%*** Model:
% B = Xg*Bg + eg. Xg=[1 1 ... 1]'; B=[B_1 B_2...B_n];
%
%
% Sergio Luiz Novi Junior, sergiolnovi@gmail.com
% May,2020

% Beta of the group is computed separated for each channel
for Nchan=1:size(betaGroup,1)
    
    % Select beta for each channel
    y = betaGroup(Nchan,:);
    y=y';
    
    % Remove nan entries from y and cov_b
     lst_nan = find(isnan(y)==1);
     y(lst_nan) = [];
    
    % Create the model (X_g). vector with ones.
    X = ones(size(y,1),1);
    
    % Run the analysis with the chosen channel is not selected to be
    % avoided
    
    if isempty(find(AvoidChan == Nchan))
        
        % Estimate beta without solving the error heterosdacity
        beta = inv((X'*X))*X'*y;
        e = y-beta*X;
        
        % Stats.covb is the first estimation of the error \sigma_g
        Stats.covb = inv((X'*X))*var(e);
        
        
        betaOld=beta;
        change=1;
        itr=0;
        
        if ~isempty(covb)
            
            % remove nan entries
             covb_Nchan = covb(Nchan,:);            
             covb_Nchan(lst_nan) = [];
            
            % Criteria to stop the loop
            while change>=0.01 & itr<5
                
                % Compute W_g
                W = eye(size(y,1))./sqrt(covb_Nchan+ Stats.covb);
                
                % Multiply both sides of the equation
                Xf = W*X;
                yf = W*y;
                % Compute New beta
                beta = inv((Xf'*Xf))*Xf'*yf;
                % Compute Beta Change
                change = abs(betaOld-beta)./abs(betaOld);
                betaOld=beta;
                % Update the error variace
                e = yf-beta*Xf;
                Stats.covb = inv((Xf'*Xf))*var(e);
                itr=itr+1;
                
            end
        end
        
        % Save computed beta.
        betaAverage(Nchan)=beta;
        % Compute t-test and p-value
        t = beta/sqrt(Stats.covb);
        p(Nchan) = (1-tcdf(abs(t),size(y,1)-1));
        VarChan(Nchan)=Stats.covb;
        tvalue(Nchan)=t;
        
    end
    
end




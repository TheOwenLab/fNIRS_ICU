function[B_k,Cov_k] = ...
    ExtractDataFromFirstLeveL(beta_use,covb_use,C,SSlist,AvailableSubjects,...
    BadChan)


for Nsub=1:size(AvailableSubjects,2)
    for Nchan=1:size(beta_use{AvailableSubjects(Nsub)},1)

            

        for Hb=1:2
            if isempty(find([SSlist,BadChan{AvailableSubjects(Nsub)}']==Nchan))
                Cov_aux = covb_use{AvailableSubjects(Nsub)}{Nchan}{Hb};
                Beta_aux = squeeze(beta_use{AvailableSubjects(Nsub)}(Nchan,Hb,:));

                % Correct for different dimension due to 
                % different number of regressors per subject
                Cs = C(1:size(Beta_aux,1));

                B_k(Nsub,Nchan,Hb) = Cs*Beta_aux;
                Cov_k(Nsub,Nchan,Hb) = Cs*Cov_aux*Cs';
                
            else 
                B_k(Nsub,Nchan,Hb) = nan;
                Cov_k(Nsub,Nchan,Hb) = nan;
            end

        end
    end
end




end
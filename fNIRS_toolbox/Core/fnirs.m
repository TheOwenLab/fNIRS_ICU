classdef fnirs
    % functional NIRS results estimates class
    properties
        t               % HRF time vector (Ntps x 1)
        concs           % hemoglobin concentration data around tHRF (Ntps x Nchn x 3 x Nevents)
        duration        % duration of stimulus
        hrf             % HRF data matrix (Ntps x Nchn)
    end
    
    methods
        function obj = RemoveBadTrials(obj,Tolerance)
            % Get %change vector at each time step
            dconcs_dt = abs( diff(obj.concs)./(1 + abs(obj.concs(1:end-1,:,:,:))));
            % Remove trials with big motion artifact. Here, motion artifact is defined
            % as sudden temporal changes (i.e., %changes in the derivative that are
            % higher than (mean + tol*std)
            for chn=1:size(obj.concs,2)
                for trial = 1:size(obj.concs,4)
                    for Hb = 1:2
                        lst2 = find( abs(dconcs_dt(:,chn,Hb,trial)) >= ...
                            mean(dconcs_dt(:,chn,Hb,trial)) + Tolerance*std(dconcs_dt(:,chn,Hb,trial)) );
                        if ~isempty(lst2) % there are sudden changes
                            obj.concs(:,chn,:,trial) = NaN; % make Hb for this channel & trial NaN
                        else
                            obj.concs(:,chn,Hb,trial) = detrend(obj.concs(:,chn,Hb,trial)); % just remove linear trends in the trial (due to slight motion); this is typically the 2nd most common source of MA
                        end
                    end
                end
            end
            obj.concs(:,:,3,:) = obj.concs(:,:,1,:) + obj.concs(:,:,2,:); %re-calculate HbT
        end
        function obj = AverageTrials(obj)
            obj.hrf = nanmean(obj.concs,4);
            % Renormalize to zero mean before baseline
            baseline = find(obj.t <= 0);
            for Hb = 1:3
                obj.hrf(:,:,Hb) = obj.hrf(:,:,Hb) - ones(size(obj.hrf,1),1)*...
                    mean(obj.hrf(baseline,:,Hb),1);
            end
        end
        
    end
end
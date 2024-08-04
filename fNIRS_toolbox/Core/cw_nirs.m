classdef cw_nirs
    % NIRS dataset class
    properties
        t               % time vector (Ntps x 1)
        d               % data matrix (Ntps x Nchn)
        s               % stimulation vector (Ntps x 1)
        aux             % auxiliary data matrix (Ntps x 4)
        SD              % source-detector positioning info structure
        StimTriggers    % triggers marked as stimulation onset
        dc              % hemoglobin concentration matrix (Ntps x Nchn/2 x 3)
        fnirs           % cell with N fnirs objects (N is the # of different stim types)
        fcnirs          % fcnirs object from fcnirs class with connectivity results from data
    end
    
    methods
        function obj = MarkBadChannels(obj,preproc)
            if isfield(preproc,'SNR_threshold') & ~isempty(preproc.SNR_threshold)
                obj.SD = enPruneChannelsLOB(obj.d,obj.SD,ones(size(obj.t,1),1),[-10 10^7],preproc.SNR_threshold,[0 100],0);
                obj.SD.BadChannels = find(obj.SD.MeasListActAuto==0);
                obj.SD.MeasListAct = obj.SD.MeasListActAuto;
            end
        end
        
        function obj = DeletePointsPrior2FirstTrigger(obj,Nseconds)
            
            
            % Take the first trigger as reference
            triggers = find(obj.s==1);
            FirstTrigger = triggers(1);
            
            % Compute the limit for deleting points
            deleteUp2Point = FirstTrigger-round(Nseconds*obj.SD.f);
            
            % Delete points from s, d, aux and t
            obj.s(1:deleteUp2Point,:)=[];
            obj.d(1:deleteUp2Point,:)=[];
            obj.aux(1:deleteUp2Point,:)=[];
            obj.t(end-deleteUp2Point+1:end,:)=[];
            
        end
        
        function SCI = ComputeScalpCouplingIndex(obj,freq)
            % Compute the SCI index to verify if a given
            % channel has an appropriate coupling.
            % The goodness of the coupling is based on
            % the correlation of the largest and smallest
            % wavelenghts in the heart-rate frequencies
            %
            % Input: 
            %   freq: [HighPassFilter LowPassFilter]. For this method,
            % the recommended range is [0.5 2.5].
            % Output:
            %   SCI: vector with the quality of each channel. 
            % a good treshold is 0.7
            %
            % Reference: Cortical activation patterns 
            % correlate with speech understanding after 
            % cochlear implantation, Olds and Polonini et al. 
            
            % Filter data in the desired frequencies
            d_filtered = obj.BPFilter(obj.d,freq);          
            d_filtered = d_filtered./max(d_filtered);
            
            % Find the index for the lowest and highest wavelengths
            IndexL = find(obj.SD.MeasList(:,4)==1);
            IndexH = find(obj.SD.MeasList(:,4)...
                == max(obj.SD.MeasList(:,4)));
 
            for Nchannel=1:size(IndexL,1)
                C = corrcoef([d_filtered(:,IndexL(Nchannel)) ...
                d_filtered(:,IndexH(Nchannel))]); 
                SCI(Nchannel) = C(1,2);    
            end
            
            
        end
        
        function dOD = Convert2OD(obj)
            dm = nanmean(abs(obj.d),1);
            nTpts = size(obj.d,1);
            dOD = -log(abs(obj.d)./(ones(nTpts,1)*dm));
        end
        
        function d_spline = SplineCorrection(obj, data, SplineThreshold, timings)
            d_spline = data;
            if isempty(timings)
                d_spline = MotionCorrectionMARA(data,obj.SD,round(2.5*obj.SD.f),SplineThreshold);
            else
                for i = 1:size(timings,1)
                    clear lst lst_nframes
                    % We have to give as INPUT to MotionCorrectionMARA
                    % 2.5 seconds before the begining of onset. In this,
                    % way we can correct the whole stimulus block.
                    n_frames = round(2.5*obj.SD.f);
                    lst_n_frames = timings(i,1)-n_frames:1:timings(i,1)-1;
                    lst = [lst_n_frames timings(i,:)];
                    d_spline(lst,:) = MotionCorrectionMARA(data(lst,:),obj.SD,round(2.5*obj.SD.f),SplineThreshold);
                end
            end
        end
        
        function d_wavelet = WaveletCorrection(obj, data, Wavelet_iqr)
            d_wavelet = hmrMotionCorrectWavelet(data,obj.SD,Wavelet_iqr);
        end
        
        function d_filter = BPFilter(obj, data, FrequencyRange)
            d_filter = hmrBandpassFiltLOB(data, obj.SD.f, FrequencyRange(1), FrequencyRange(2));
        end
        
        function [dc_PCA, svs, nSV] = PCAFilter(obj,preproc)
            tInc = ones(length(obj.dc),1)';
            dc = obj.dc;
            dc = permute(dc,[1 3 2]);
            [dc_PCA,svs,nSV] = enPCAFilterLOB(dc,obj.SD,tInc,preproc.PCAFilter_nSV);
            dc_PCA = permute(dc_PCA,[1 3 2]);
        end
               
        function dc = Convert2Conc(obj, data, ppf)
            %Hb = hmrOD2Conc(data, obj.SD, ppf);
            Hb = hmrOD2ConcLOB(data,obj.SD,ppf);
            Hb = Hb*1e6;
            dc = permute(Hb,[1 3 2]);
        end
               
        function obj = DiscardBadChannels(obj)
            % Make Hb data from bad channels NaN so that one does not
            % analyze bad data
            for chn = 1:size(obj.SD.BadChannels,1)
                if obj.SD.BadChannels(chn) > size(obj.SD.MeasList,1)/2
                    obj.dc(:,obj.SD.BadChannels(chn)+size(obj.SD.MeasList,1)/2,:) = NaN;
                else
                    obj.dc(:,obj.SD.BadChannels(chn),:) = NaN;
                end
            end
        end
        
        function [trials, timings] = GetEventData(obj, fnirs)
            % Get event trials
            trials = find(obj.s(:,fnirs.StimType) == 1 );
            % Get timing blocks for each trial
            for i = 1:length(trials)
                timings(i,:) = trials(i) - round(fnirs.t_baseline*obj.SD.f):...
                    trials(i) + round(fnirs.t_recovery*obj.SD.f);
                % Check whether there is datapoints for last timing
                if timings(i,end) > length(obj.t)
                    timings(i,:) = [];
                    trials(i) = [];
                    disp('NOTE: not enough data points for last trial. Discarding last trial...')
                end
                
            end
        end
        
        function [tHRF, concs] = CombineTrials(obj, fnirs)
            % Recalculate trials and timings so that function can stand
            % alone.
            [trials, timings] = GetEventData(obj, fnirs);
            % Combine Hb data onto a 4-D data matrix
            for trial = 1:size(timings,1)
                concs(:,:,:,trial) = obj.dc(timings(trial,:),:,:);
                % Normalize each trial to zero mean before stimulus onset, for
                % each hemoglobin concentration
                for Hb = 1:3
                    concs(:,:,Hb,trial) = concs(:,:,Hb,trial) - ones(size(concs,1),1)*...
                        mean(concs(1:round(fnirs.t_baseline*obj.SD.f),:,Hb,trial),1);
                end
            end
            % Create temporal vector around the trial timing
            tHRF = -fnirs.t_baseline:1/obj.SD.f:fnirs.t_recovery;
            if length(tHRF) ~= size(concs,1)
                tHRF(end+1:size(concs,1)) = tHRF(end)+1/obj.SD.f;
            end
            tHRF = tHRF';
        end
        
        function obj = RemoveAutoCorrelation(obj)
            %Removes autocorrelation from Resting State Data
            %with prehitenning methodology
            % For details Refer to:
            %1 - Characterization and correction of the
            %false-discovery rates in resting state connectivity
            %using functional near-infrared spectroscopy
            %
            %2 - Autoregressive model based
            %algorithmfor correcting motion and
            %seriallycorrelated errors in fNIRS
            
            
            Pmax = 50;
            
            % Time Series Length
            n = size(obj.dc,1);            
            
            for Hb=1:2
            
             for Nchannel=1:size(obj.dc,2)                
                
              if isempty(find(isnan(obj.dc(:,Nchannel,Hb))==1))
                   
                clear y yf a vt 
                % Take Original Time Series                  
                y = obj.dc(:,Nchannel,Hb);
                
                for P=1:Pmax
                  % For a given parameters P we find the coefficients that
                  % minimize autoregressive model (AR(P));
                  a = aryule(y,P);
                   
                  % Once we have the parameters a, we can filter the error
                  % to find the new non atucorrelated error (vt).
                  vt = filter(a,+1,y);
                    
                  % Next, we can compute the baysian information
                  % criterion (BIC(P)).
                    
                  % Log Likelihood
                  LL = -1*(n/2)*log( 2*pi*mean(vt.^2))+...
                     -0.5*(1/mean(vt.^2))*sum(vt.^2);
                    
                  % Baysian information
                  BIC(P) = -2*LL+P*log(n);
                 end
                
                 %Optimal is the P that minimizes BIC
                 [~,OptimalP] = min(BIC);
    
                 AR_Parameters = aryule(y,OptimalP); %Find parameters
                
                 % Filter y
                 yf = filter(AR_Parameters,+1,y);
    
                 % Update DC
                 obj.dc(:,Nchannel,Hb)=yf;
                
                 % Save OptimalP for double checking
                 obj.SD.Optimal_P(Nchannel,Hb) = OptimalP;

                end
                
              end
                
            end
            
       
            
        end
        
        function [betaChan,pValue,covb] = GLM_Prewhitening...
            (obj,duration,lstSS,SC_exclude,HbO_regressor, HbR_regressor,...,
            AdditionalRegressors,opt)
            % Perform GLM analysis with PreWhitening
            [betaChan,pValue,covb] = GLM_WithPrewhiteningRobustFit...
            (obj,duration,lstSS,SC_exclude,HbO_regressor, HbR_regressor,...
            AdditionalRegressors,opt);
        end
        
        function [filtered_dc,Stats] = ...
                PerformPhysiologyRegression...
                (obj,SSlist,AdditionalRegressors,CV,opt,...
                flag_time_SC,flag_time_AD)
                   
            % This function performs physiological Regression based 
            % on Short-channel and with additional Regressors of physiology. 
            % The additional Regressors must be inputed.  
            % The regressions will be performed with the robustfit 
            % matlab function:
            % robustfit(X,y,'bisquare',[],'on');
            % FOR FULL DOCUMENTATION: "help PhysiologyRegression_GLM". 
            [filtered_dc,Stats] = ...
                PhysiologyRegression_GLM...
                (obj,SSlist,AdditionalRegressors,CV,opt,...
                flag_time_SC,flag_time_AD);        
        
        end
        
        function s = CreateFakeTriggers2RestingStateData...
                (obj,tfirst,tlast,Rest)
            
            % Add fake triggers in order to perform simulations
            % Input:
            %   obj - cw_nirs Object
            %   tfirst - First trial in seconds. 
            %   tlast - Last trial in seconds. If you do not want to use it, 
            % insert tlast as empty [].
            %   Rest - Vector with the maximum and minimum period of rest 
            %[min max]
            %
            % Output: 
            %   s - Vector of stimulus (only ONSET).            
            
            f = obj.SD.f;
            if ~isempty(obj.dc)
                s = zeros(size(obj.dc,1),1);
            else
                s = zeros(size(obj.d,1),1);
            end
            % First Trial
            trial = tfirst; 
    
            
            if isempty(tlast)
              tlast=0;
            end
            
            % Add trigger of the first trial
            s(round(trial*f))=1;    
    
            while trial <= (size(obj.d,1)/f)-tlast
        
                trial = randi([(trial+Rest(1)) (trial+Rest(2))]);
        
                % Make sure that the last trial has a proper rest 
                % time before  the end of the task
                if trial <= (size(obj.d,1)/f) - Rest(2)
                    s(round(trial*f))=1;
                end
                
            end
    
            
            
        end
                
        function obj = CreateSemiSyntheticSimulation(obj,Params)
            % Simulates hemodynamic response on resting state data. 
            % The HRF is simulated in the hemoglobin data, thereby
            % the cw_nirs object needs to have the dc calculated. 
            %
            % INPUTS: 
            %       cw_nirs object with dc calculated.
            %       Params - Params is a structure that contains
            % info for the simulations. The fields of Params are: 
            %   Params.opt: 
            %    opt = 0, the simulated hrf of each channel
            % will all have the same amplitude. This amplitude
            % must be specified in the fieds:.HbO and .HbR (see below). 
            %    opt = 1, the simulated hrf will have amplitude depending
            % on the mean absolute deviation of each hemoglobin. The CNR
            % of each hemoglobin needs to be especified in the Params.HbO
            % and Params.HbR.
            
            %   Params.HbO - Amplitude of HbO in micro MOL for opt=0, or
            % CNR for opt = 1;
            %   Params.HbR - Amplitude of HbR in micro MOL for opt=0, or
            % CNR for opt = 1; (For activation input negative for HbR.)
            %   Params.tfirst - time (seconds) for the first trigger.
            %   Params.tlast - time (secons) fot the last trigger. If you 
            % do not want to set tlast, you can input it as empty ([]). 
            %   Params.Rest - Vector with the minimum and maximum rest time
            %between trials [min max]. It should be in seconds.  
            %   Params.StimDuration - Duration of the task in seconds. 
            %   Params.ChosenChannels - Vector with the channels to include
            %   the simulated HRF. 
            %
            
            % Create vector of stimulus
            s_new =  obj.CreateFakeTriggers2RestingStateData...
                (Params.tfirst,Params.tlast,Params.Rest);
            obj.s = s_new;
            
            % Create HRF based on the new vector of stimulus 
            X = CreateDeasingMatrix(obj,s_new,Params.StimDuration);            
            
           % All channels are simulated with the same amplitude
           if Params.opt==0 
               
            % Add hrf in Resting State data with the proper amplitudes 
            
            %HbO            
            obj.dc(:,Params.ChosenChannels,1) = ...
               obj.dc(:,Params.ChosenChannels,1) + ...
               Params.HbO*X;
            
            %HbR                     
            obj.dc(:,Params.ChosenChannels,2) = ...
               obj.dc(:,Params.ChosenChannels,2) + ...
               Params.HbR*X;
           
           end
           
           % Amplitude of each channel is based on its median
           %absolute deviation.
           if Params.opt==1               
               
               
               % Add hrf in Resting State data with the proper amplitudes 

               % HbO
               AmpHbO = 1.4826*Params.HbO*mad(obj.dc(:,:,1));
               X_HbO = X*AmpHbO;
               
               obj.dc(:,Params.ChosenChannels,1) = ...
               obj.dc(:,Params.ChosenChannels,1) + ...
               X_HbO(:,Params.ChosenChannels);
               
               % HbR
               
               AmpHbR = 1.4826*Params.HbR*mad(obj.dc(:,:,2));
               X_HbR = X*AmpHbR;
               
               obj.dc(:,Params.ChosenChannels,2) = ...
               obj.dc(:,Params.ChosenChannels,2) + ...
               X_HbR(:,Params.ChosenChannels);
                              
           end
           
           
           
           
        end
                
                
end


end


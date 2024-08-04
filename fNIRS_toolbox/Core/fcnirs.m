classdef fcnirs
    % functional connectivity (network-based) NIRS class
    properties
        SD               % SD structure from standard cw_nirs
        CorrMtx          % Similarity Matrix from cw_nirs data (Correlation matrix, correntropy matrix, etc.)
        
        CorrMtxEachTrial % Similarity Matrix of each trial. If it is 
                         % Resting state data then CorrMtxEach Trial is 
                         % is equal to CorrMtx   
        
        p                % Correlation threshold vector
        NetworkProp
        
        concs            % Concentration per Trial (if Resting State 
                         % Concs is the same as dc from cw_nirs
        
        r                % Correlation threshold (single value)
        
        G                % Binary Graph based on HbO and HbR
        
        GraphProperties  % Graph metrics from G
        
    end
    
    methods
        
        function obj = GlobalNetworkProperties(obj)
            % Create matrices for common properties
            av_k = zeros(size(obj.CorrMtx,3),length(obj.p));  % average degree
            clustering = zeros(size(obj.CorrMtx,3),length(obj.p));  % clustering coefficient
            std_k = zeros(size(obj.CorrMtx,3),length(obj.p));  % clustering coefficient
            % Calculate properties from basic functions
            for Hb = 1:size(obj.CorrMtx,3)
                for threshold = 1:length(obj.p)
                    [A,av_k(Hb,threshold),K1(:,threshold)] = ...
                        adjacency_matrix_av_degree(obj.p(threshold),obj.CorrMtx(:,:,Hb));
                    Distances = distance_bin(A);
                    [L(Hb,threshold),L_std(Hb,threshold)] = average_pathlength(Distances);
                    D(Hb,threshold) = max(max(Distances));
                    clustering(Hb,threshold) = mean(clustering_coef_bu(A));
                    links(Hb,threshold) = link_number(A);
                end
                Pk(:,:,Hb) = K1;
                std_k(Hb,:) = std(K1);
                clear D_* A_* map i x;
            end
            % Renormalize by real number of valid channels
            N_nodes = size(obj.CorrMtx,1) - length(obj.SD.BadChannels) - 1;
            % Save properties in obj
            obj.NetworkProp.K_avg = av_k'/N_nodes;
            obj.NetworkProp.CC = clustering';
            obj.NetworkProp.Pathlength = L';
            obj.NetworkProp.Std_K = std_k';
            obj.NetworkProp.Dist_K = Pk;
            obj.NetworkProp.N_links = links';
                    
        end
        
        function obj = ComputeCorrMtx(obj,lag_opt,max_lag)
            % Compute CorrMtx.
            %   INPUT: 
            %        lag_opt = 0 -> perform corrcoef
            %        lag_opt = 1 -> perform xcorr 
            %        max_lag is the maximum allowed lag in seconds.     
            %  1 - If there are more than one trial, CorrMtx will be
            % the average correlation matrix of each trial. The averages
            % are compute in the Fisher space.
            %  2 - If there is only one single trial, CorrMtx will just be
            % the correlation matrix opf this specif trial.
            
            % Convert max lag to frames
            if lag_opt==1
               
                max_lag = round(obj.SD.f*max_lag);
                
            end
            
            
            for hb=1:3
                
                % Run Correlation for all trials
                for Ntrial=1:size(obj.concs,4)
                    
                    % Compute Correlation per trial
                    data = obj.concs(:,:,hb,Ntrial);
                    
                    if lag_opt==0
                        C = corrcoef(data);
                    end
                    
                    if lag_opt==1
                        
                        for Nline = 1:size(data,2)
                            
                           for Ncollum = 1:size(data,2) 
                            
                                % save data into x and y to facilitate
                                % nomenclature                                
                                x = data(:,Nline);
                                y = data(:,Ncollum);
                            
                                auxCorr = xcorr(x,y,max_lag,'coeff');
                            
                                [value index] = max(abs(auxCorr));
                            
                                C(Nline,Ncollum) = auxCorr(index);
                           end                           
                        end
                    end
                    
                    % Apply fisher tranform
                    Cf(:,:,Ntrial) = atanh(C);
                    
                    % Save CorrMtx of each trial
                    C_eachTrial(:,:,Ntrial,hb) = C;
                    
                    clear data
                    clear C;
                end
                clear Ntrial;
                %if size(size(obj.concs),2)==4
                % Average Correlation from all trials
                Cf_av = mean(Cf,3);
                clear Cf;
                
                % Back to Correlation space (Inverse Fisher)
                Cf_av =  tanh(Cf_av);
                
                % Set main diagonal to zero
                for Nchn =1:size(Cf_av,1)
                    Cf_av(Nchn,Nchn)=0;
                end
                clear Nchn
                
                % Remove NAN by assigning it to zero
                lst = find(isnan(Cf_av)==1);
                Cf_av(lst)=0;
                clear lst;
                
                CorrMtx(:,:,hb) = Cf_av;                
                clear Cf_av;
            end
                        
            obj.CorrMtx = CorrMtx;            
            obj.CorrMtxEachTrial = C_eachTrial;
            
        end
        
        function obj = ExtractBinaryGraphHbOHbR(obj)
            % Create unweighted graph based on correlation 
            % matrices from HbO and HbR (fcnirs.CorrMtx).
            % links will be assigned to 1 if both HbO and HbR weights 
            % are higher than a correlation threshold
            % set in the property obj.r and 0 otherwise.
            % 
            % INPUTS: obj - fcnirs objects
            %          
            
            % Step 1: Binarize Correlation Matrix based on a fixed 
            % threshold independently for each hemoglobin. 
            
            for hb=1:2
                [A(:,:,hb),~] = ...
                    adjacency_matrix(obj.r,squeeze(obj.CorrMtx(:,:,hb)));
            end
            
            % Step 2: Create Binary Graph G based on the adjacency matrix 
            % from HbO and HbR. 
            obj.G = adjacency_matrix(2,sum(A,3));
                        
        end
     
        function obj = ComputeGraphProperties(obj)
            % Compute several network properties, such as 
            % degree, clustering coefficient and betweeness centrality
            % from the Binary Graph G.
                       
            % Node degree                   
            obj.GraphProperties.degree = sum(obj.G);   
            
            % Clustering Coefficient
            obj.GraphProperties.clustering = clustering_coef_bu(obj.G);
            
            % Betweenness Centrality 
            obj.GraphProperties.betweenness = ...
                centrality(graph(obj.G),'betweenness'); 
            
        end
        
        
        
    end
end
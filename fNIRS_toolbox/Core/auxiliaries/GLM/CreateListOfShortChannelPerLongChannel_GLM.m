function[SC_per_LongChannel] = CreateListOfShortChannelPerLongChannel_GLM...
    (obj,lstSS,opt)
% INPUT:
%   obj: cw_nirs obj
%   lstSS: list of short channels
%   opt: 0 - Use the closest short channel
%        1 - Do nothing
%        2 - Use the SC with highest correlation with the Long channel
%
% OUTPUT:
%   SC_per_LongChannel: Vector with the correspondent SC per long channel.

    
% Variable to save the SC of each long channel 
SC_per_LongChannel = zeros(size(obj.dc,2),1);

%Opt = 0
if opt==0
    % Find Nearest Short channels for each long channel
    [channelspershortch] = Nearshortchanel_GLM(lstSS,obj.SD);
    
    % Lets flip the cell and create a vector (SC_per_LongChannel) with
    % the number of the SC per Long channel
    
    
    for Nsc=1:size(channelspershortch,2)
        indexes = channelspershortch{Nsc};
        for dummy=1:size(indexes,2)
            SC_per_LongChannel(indexes(dummy)) = lstSS(Nsc);
        end
    end
    
end


if opt==2
    
    % Compute Pearson correlation with among the SCs and Long channels
    % Then choose the SC with highest correlation
    % CORRELATION IS BEING BASED ON HBO
    ShortCorrelation = zeros(size(obj.dc,2),size(lstSS,2));
    
    for NLong = 1:size(obj.dc,2)
        
        for Nshort=1:size(lstSS,2)
            
            aux = corrcoef(obj.dc(:,NLong,1),obj.dc(:,lstSS(Nshort),1));
            ShortCorrelation(NLong,Nshort) = abs(aux(1,2));
            
        end
        
        [~,index] = max(ShortCorrelation(NLong,:));
        SC_per_LongChannel(NLong,1) = lstSS(index);
        
    end
end



end
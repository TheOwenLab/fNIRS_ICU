function[lstList] = findShortChannels_GLM(SD,threshold)

% Create a list of short channels
% Input: 
%   SD - matrix
%   threshold - threshold (source-detector distance)
%   for finding the short channels  
%
% Output:
%   lstList: List of short channels.



%Number of channels
nCh = size(find(SD.MeasList(:,4)==1),1);

for Nchannel=1:nCh
    
    Source = SD.MeasList(Nchannel,1);
    Detector = SD.MeasList(Nchannel,2);
    
    
    Vec = SD.DetPos_3d(Detector,:)- SD.SrcPos_3d(Source,:);
            
    pho(Nchannel) = sqrt(Vec*Vec');
        
end


    lstList = find(pho<threshold);

end
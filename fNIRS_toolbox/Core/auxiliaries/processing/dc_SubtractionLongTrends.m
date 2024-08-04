function[dc_subtracted,dc_long] = dc_SubtractionLongTrends(dc,f,Method,period)

% Function to Remove long term trends from the hemoglobin concentration 
% changes. 
%
% INPUT: 
%   dc: hemoglobin concentration changes
%   f: frequency of acquisition 
%   Method: 1 - smoothing averaged;
%           2 - Spline interpolation (to be added);
%           3 - High Band Pass filter (to be added).
%   period: period in seconds to perform smoothing average; 
%           if period is empty then period = 180 seconds.
%
% OUTPUT:
%   dc_subtracted: hemoglobin changes after subtracttion of long term 
%                  trends.
%   dc_long: hemoglobin long term trends.
%
%
%
%
% Campinas January 16th, 
% S. L. Novi (sergiolnovi@gmail.com)
% Laboratory of Biomedical Optics 
% University of Campinas


for Nchn=1:size(dc,2) %Run over all channels
    for hb=1:2 % HbO and HbR
        
        % Smoothing method
        if Method==1
            if isempty(period)
                period= 180;
            end
            dc_long(:,Nchn,hb) = smooth(dc(:,Nchn,hb),round(f*period));
        end
        
    end
end

% Subtract long term trends from dc
dc_subtracted = dc(:,:,1:2) - dc_long(:,:,1:2);
% Compute HbT
dc_subtracted(:,:,3) = dc_subtracted(:,:,1) + dc_subtracted(:,:,2);

% Computed dc_long HbT
dc_long(:,:,3) = dc_long(:,:,1) + dc_long(:,:,2); 

end

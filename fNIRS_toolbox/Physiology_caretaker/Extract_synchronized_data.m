function[CO2_data,vital_data] = ...
    Extract_synchronized_data(CO2_data, vital_data, trigger_phys)


% 1 - CO2 Data

% Find window to be kept
lst = min(find(CO2_data.t>=trigger_phys));

CO2_data.data(1:lst) = [];
CO2_data.t(1:lst) = [];

clear lst;


% 2 - Vital Data


% Find when trigger happened
lst = ...
    min(find(...
    vital_data.Raw_data.t(:,1)>=trigger_phys));

% Map the trigger in ms
TimePoint_in_ms = ...
    vital_data.Raw_data.t(lst,2);

% Extract the correct position for the resampled data
lst = ...
    min(find(vital_data.t>=TimePoint_in_ms));

vital_data.data(1:lst,:) = [];
vital_data.t(1:lst) = [];


end
% Example of how to extract the physiological data 
% acquired with the caretaker system (exported by the tablet)

% folder that has all the data
folder_path = './sampleData/';

% 1 - Extract Trigger
[trigger] = Read_extract_trigger_CareTaker(folder_path);

% 2 - end-tidal CO2
[Co2_data, time] = Read_extract_etCO2_CareTaker(folder_path,4,1);

% 3 - Blood pressure, Heart-rate, and respiration 
[vital_data] = ...
    Read_extract_data_vitals(folder_path,4,1);

close all

% Map the trigger point to the blood pressure data
trigger_pressure_aux = ...
    min(find(vital_data.Raw_data.t(:,1)>=trigger));

time_trigger_pressure = ...
    vital_data.Raw_data.t(trigger_pressure_aux,2);

pos_trigger_vital = min(find(vital_data.t>=time_trigger_pressure));

s_pressure = zeros(size(vital_data.t));
s_pressure(pos_trigger_vital) = 1;

plot(vital_data.t/(1000),vital_data.data(:,3),'-k');
hold on
plot(vital_data.t/(1000),80*s_pressure,'-r');
ylim([60 90])

% lst = min(find(time>=trigger));
% 
% s = zeros(size(time));
% s(lst) = 42;
% 
% 
% plot(time,Co2_data,'-b');
% hold on;
% plot(time,s,'-k');




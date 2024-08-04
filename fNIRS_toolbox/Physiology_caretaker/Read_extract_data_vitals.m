function[resampled_data] = ...
    Read_extract_data_vitals(folder_path,f_phys,plot_flag)
%
% Code to clean and organize data inside "data-vitals".
% Data-vitals includes blood pressure, heart rate, and respiration.
% 
%
% [Systolic, Diastolic, MAP, HeartRate, Respiration] = ...
% Read_extract_data_vitals(folder_path)
%
% INPUT:
%       Folder path in which the data is located: 'path/'
%       f_phys - desired frequency for physiological parameters
%
% OUTPUT:
%       resampled_data: resampled_data is a structute that 
% has the physiolgical parameters and time vector (in ms), both at the
% desired frequency defined by f_phys. Also, resampled_data
% has the original (raw data) to recover the time stamps and double
% check the resampled data. 
%
% Optional Input:
%   [Systolic, Diastolic, MAP, HeartRate, Respiration] = ...
%    Read_extract_data_vitals(folder_path,f_phys,1)
%   plot the updated and raw data to compare them.

if nargin<3
    plot_flag = 0;
end

% Get the name of the data
files = dir([folder_path '*_vitals_*.csv']);

% Open data file
file_vitals = fopen([folder_path files.name],'r');

% We will read line by line

% For the data-vitals, the first line is the heading:
% Date,Time,Systolic (mmHg), ...
% The data starts at the second line

% Read first line
LineContent = fgetl(file_vitals);

% Save headings
head_dummy = split(LineContent,',');
Raw_data.label_data{1} = head_dummy{3};
Raw_data.label_data{2} = head_dummy{4};
Raw_data.label_data{3} = head_dummy{5};
Raw_data.label_data{4} = head_dummy{6};
Raw_data.label_data{5} = head_dummy{7};

clear head_dummy;

% Read second line
LineContent = fgetl(file_vitals);

% Create variables extract the raw data

% vitals raw data (dummy variable)
vitals_data_dummy = [];

% Raw time vector (dummy variable)
t_dummy = [];

while isnumeric(LineContent)==0
    
    % Split line content based on the commas
    S_Line = split(LineContent,',');
    
    % 1 - Get vital numerical values
    %
    % cell 3 - Systolic (mmHg)
    % cell 4 - Diastolic (mmHg)
    % cell 5 - MAP (mmHg)
    % cell 6 - HeartRate (bpm)
    % cell 7 - Respiration (Bpm)
    aux_vital_dummy = [str2num(S_Line{3}),str2num(S_Line{4}),...
        str2num(S_Line{5}),str2num(S_Line{6}),...
        str2num(S_Line{7})];
    
    % Concatenate the data
    vitals_data_dummy = [vitals_data_dummy;aux_vital_dummy];
    
    clear aux_vital_dummy;
    
    % 2 - Get temporal information: hour and mileseconds stamps
    
    % Split Time Stamp
    S_time_stamp = split(S_Line{2},':');
    
    % Convert hour to seconds
    Time_in_seconds = ...
        3600*str2num(S_time_stamp{1}) + ...
        60*str2num(S_time_stamp{2}) +...
        str2num(S_time_stamp{3});
    
    % Get the samps in ms
    Time_in_ms = str2num(S_Line{10});
    
    % concatenate temporal data
    t_dummy = [t_dummy;[Time_in_seconds,Time_in_ms]];
    
    % Read line at each iteration
    LineContent = fgetl(file_vitals);
    
end

% close the opened files
fclose all;

clear Time_in_ms Time_in_seconds...
    LineContent files file_vitals ...
    S_Line S_time_stamp folder_path;

% Save Raw data
Raw_data.t = t_dummy;
Raw_data.data = vitals_data_dummy;


% Create a regular vector time with constant frequency
% The frequency will be sampled a f_phys
NumberOfPoints = round(f_phys*(t_dummy(end,2) - t_dummy(1,2))/(1000));

time = linspace(t_dummy(1,2),...
    t_dummy(end,2),...
    NumberOfPoints);

time = time';

% Compute regular data_vitals for the time-points in "time"
%
vitals_data_new = interp1(t_dummy(:,2),vitals_data_dummy,time);

% Organize all data into a single structure
resampled_data.label_data = Raw_data.label_data;
resampled_data.t = time;
resampled_data.data = vitals_data_new;
resampled_data.f = f_phys;
resampled_data.Raw_data = Raw_data;

% Plot data to double check
if plot_flag==1
    
    for Nfigure = 1:size(vitals_data_new,2)
                
        figure(Nfigure)
        plot(t_dummy(:,2),vitals_data_dummy(:,Nfigure),'-k');
        hold on;
        plot(time,vitals_data_new(:,Nfigure),'--m');
        legend('raw data','Corrected Data');
        title(Raw_data.label_data{Nfigure});
                
    end
end



end
function[time_1,time_2] = Read_extract_trigger_CareTaker(folder_path)
% Code to extratc the trigger used in the experiment
% This function was built with the assumption there is only
% the start trigger. It may need to be revised in the future
% for protocols that have more than one trigger.
%
% INPUT:
%   folder_path: Folder path in which the data is located: 'path/'
%
% OUTPUT:
%   time_1: Time in seconds of the trigger (number).
%   time_2: Time in the formart: Hour: Minute: Second (str).


% Get the name of the data
files = dir([folder_path '*events*.csv']);

% Open data file
file_events = fopen([folder_path files.name],'r');

% We will read line by line. The important data is in the 4th line

% 1st line: Heading
LineContent = fgetl(file_events);

% 2dn line:
LineContent = fgetl(file_events);

% 3rd line:
LineContent = fgetl(file_events);

% 4th line:
LineContent = fgetl(file_events);

% Split line content based on the commas
S_Line = split(LineContent,',');

% Convert data to seconds

% Split Time Stamp
S_time_stamp = split(S_Line{2},':');

% Trigger in seconds
time_1 = ...
    3600*str2num(S_time_stamp{1}) + ...
    60*str2num(S_time_stamp{2}) +...
    str2num(S_time_stamp{3});

% Original trigger time format
time_2 = S_Line{2};

end
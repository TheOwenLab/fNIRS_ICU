function[Co2_data, time] = ...
    Read_extract_etCO2_CareTaker(folder_path,Frequency,plot_flag)
%
% Code to clean and organize CO2 data
% This extracts the data acquired with the capnograph
% that is connected with the care taker system.
%
% [Co2_data, time] = Read_extract_etCO2_data(folder_path)
%
% INPUT:
%       Folder path in which the data is located: 'path/'
%       Frequency - Specify desired frequency of CO2 outpu
%
%
% OUTPUT:
%       Co2_data: end-tidal CO2 data
%       time: vector in seconds
%
% Optional Input:
%   [Co2_data, time] = ...
%    Read_extract_etCO2_data(folder_path,1)
%
%   plot the CO2 data from each each "step".

if nargin<2
    plot_flag = 0;
end

% Get the name of the data
files = dir([folder_path '*etco2_waveform*.csv']);

% Open data file
file_CO2 = fopen([folder_path files.name],'r');

% We will read line by line

% For the CO2 data, the first line is the heading:
% Date, Time, ETCO2 (mmHg)
% The data starts at the second line

% Read first line
LineContent = fgetl(file_CO2);
% Read second line
LineContent = fgetl(file_CO2);


% Create variables extract the raw data

% CO2 raw data (dummy variable)
Co2_data_dummy = [];

% Raw time vector (dummy variable)
t_dummy = [];

while isnumeric(LineContent)==0
    
    % Split line content based on the commas
    S_Line = split(LineContent,',');
    
    % 1 - Get Co2 numerical values
    Co2_data_dummy = [Co2_data_dummy;str2num(S_Line{3})];
    
    % 2 - Convert data to seconds
    
    % Split Time Stamp
    S_time_stamp = split(S_Line{2},':');
    
    Time_in_seconds = ...
        3600*str2num(S_time_stamp{1}) + ...
        60*str2num(S_time_stamp{2}) +...
        str2num(S_time_stamp{3});
    
    % concatenate data
    t_dummy = [t_dummy;Time_in_seconds];
    
    % Read line at each iteration
    LineContent = fgetl(file_CO2);
    
end

% close the opened files
fclose all;

% Correct for repetitive time points
[t_dummy_unique,IA,~] = unique(t_dummy);

Co2_data_dummy_unique = Co2_data_dummy(IA);

% Create a regular vector time with constant frequency

if isempty(Frequency)
    time = ...
        linspace(min(t_dummy_unique),...
        max(t_dummy_unique),...
        size(t_dummy_unique,1));
    
else    
    time = ...
        linspace(min(t_dummy_unique),...
        max(t_dummy_unique),...
        round(Frequency*(max(t_dummy_unique)-min(t_dummy_unique))));
    
end

time = time';

% Compute regular CO2_data for the time-points in "time"
Co2_data = interp1(t_dummy_unique,Co2_data_dummy_unique,time);

% Plot data to double check
if plot_flag==1
    
    figure
    plot(t_dummy,Co2_data_dummy,'-r');
    hold on;
    plot(t_dummy_unique,Co2_data_dummy_unique,'--b');
    plot(time,Co2_data,'--m');
    legend('raw data','Unique data','Corrected Data');
    
end


% Check the frequency of acquisition
cnt=1;
for N = 2:size(time,1)
    
    f_dummy(cnt) = 1/(time(N)-time(N-1));
    
    cnt = cnt+1;
end






end
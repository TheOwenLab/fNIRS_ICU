function probefile = GetProbeFileName
%
% This function will make a window pop up asking for a probe layout file.
% It can read .sd/.layout/nirsInfo.mat files. Once selected, the layout
% file path will be passed as output.
% This function is currently used to ready probe layouts for both DCS and
% NIRS reading functions.

uiwait( msgbox({'Missing a probe configuration for this dataset.'...
    'Please select the probe configuration file (.layout / nirsInfo.mat / .sd) for the dataset'}, ...
    'Warning','warn') );
[probelayout,probepath] = uigetfile('*.layout;*.mat;*.sd','Choose the file containing SD info');
probefile = [probepath probelayout];
end
% ========================================================================
% Script to Check the Quality of Short Channels in fNIRS Data
% ========================================================================
% This script loads raw fNIRS data and plots the signal quality of 
% predefined short-separation (SS) channels. For each SS channel, 
% it shows:
%   - The power spectral density (using Welch's method, pwelch)
%   - The raw time series signal
%
% HOW TO USE:
% 1. Ensure that fNIRS-main and it's sub directories are added to the path
% 2. Make sure mat file containing dat is available in the working 
%    directory.
% 3. Update the 'SSlist' to reflect the indices of the short channels 
%    you want to inspect.
% 4. Set 'Nsub' to the subject index (from dataNIRS) you want to analyze.
% 5. Run the script — each SS channel will open in its own figure with
%    two subplots (PSD and raw time series).
%
% ========================================================================

clear
close all

% ------------------------------------------------------------------------
% Load dataset
% ------------------------------------------------------------------------
% The .mat file must contain the variable 'dataNIRS' structured as:
%   dataNIRS{subject}{session}

load Controls_MI

% ------------------------------------------------------------------------
% Define list of Short Channels (SC)
% ------------------------------------------------------------------------
% These are channel indices identified as short-separation detectors
SSlist = [8 29 52 66 75 92 112 125];

% ------------------------------------------------------------------------
% Define available subject indices
% ------------------------------------------------------------------------
% Not all subjects are usable — this list excludes poor-quality datasets
AvailableSubjects = [1 9 10 12 15 17 18 20 21 22 24 26:31 33 34 35 37:40];

% ------------------------------------------------------------------------
% Select subject to analyze
% ------------------------------------------------------------------------
% Subject's short channels are inspcted one at a time
% Change this number to inspect a different subject
Nsub = 40;

% Extract session data for chosen subject
r = dataNIRS{1,Nsub}{1,1};

% ------------------------------------------------------------------------
% Loop over each Short Channel and plot PSD + raw signal
% ------------------------------------------------------------------------

for Nchan = SSlist

    figure(Nchan) % open a new figure per channel
        
    % --- Power spectral density ---
    subplot(1,2,1)
    pwelch(r.d(:,Nchan+3*129), size(r.d,1), [], [], r.SD.f);
    title(['PSD - SC ' num2str(Nchan)])
        
    % --- Raw time series ---
    subplot(1,2,2)
    plot(r.d(:,Nchan+3*129));
    xlim([1 size(r.d,1)])
    title(['Time Series - SC ' num2str(Nchan)])
end

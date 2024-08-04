function DataInfo = GroupData(varargin)

% ----------------------------------------------------------
% ------ PART 1: CONVERT NIRx TO .nirs DATA & ARRANGE DATA FOR ANALYSIS 
% ------         (DO IT ONCE ONLY) -----------
% ----------------------------------------------------------

if nargin > 0   % filename was passed
    files{1} = varargin{1};
    probe_file{1} = varargin{2};
    if nargin > 2
        MC_file = varargin{3};
    end
else    % find file to organize
    disp('Choose the .HDR/.NIRS file you want to analyze in the popup window')
    files = uigetfile_n_dir;
    probe_file = {};
end

DataInfo = cw_nirs;
for i = 1:length(files)
    clear DataTemp
    % Check if .HDR or .NIRS file
    ext = strfind(files{i},'.');
    if strcmp(files{i}(ext:end),'.hdr')
        % Convert NIRx file to .nirs
        if isempty(probe_file)
            disp('Pick the .layout OR *nirsInfo.mat file for data: ')
            disp(files{i})
            probe_file = uigetfile_n_dir;
        end
        nirs_file{i} = ReadNIRx2TechEn(files{i},probe_file{1});
    
    elseif strcmp(files{i}(ext:end),'.nirs')
        nirs_file{i} = files{i};
    
    else
        disp('Not a valid fNIRS file.')
        break
    end
    % Load .nirs data
    DataTemp = load(nirs_file{i},'-mat');
    
    % Check if it is originated from BrainSight systems
    if isfield(DataTemp,'brainsight')    
        FixBrainSightNIRSformat;
    end
    % Delete useless data points
    DataTemp = DiscardDataPoints(DataTemp);
    
    % Check stimulation vector for run
    foo = find(DataTemp.s(:)==1);
    if ~isempty(foo)
        [DataTemp.s,DataTemp.StimTriggers] = ConvertTrigger2Stim(DataTemp.t,DataTemp.s); 
    end
    
    % Load MC info if available
%    disp('Pick the Monte Carlo simulation for the probe, if available')
%    MC_file = uigetfile_n_dir;
%    if ~isempty(MC_file)
%        load(MC_file{1},'-mat')
%        DataTemp.fwMC = fwMC;
%    end
    

    % Load .nirs file onto structure
    DataInfo.t = DataTemp.t;
    DataInfo.d = DataTemp.d;
    DataInfo.s = DataTemp.s;
    DataInfo.aux = DataTemp.aux;
    DataInfo.SD = DataTemp.SD;
    if isfield(DataTemp,'StimTriggers')
        DataInfo.StimTriggers = DataTemp.StimTriggers;
    else
        DataInfo.StimTriggers = [];
    end
    
    close all
    disp('Done! Have a great day!')
end




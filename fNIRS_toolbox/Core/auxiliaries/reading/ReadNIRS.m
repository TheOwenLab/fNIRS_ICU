function varargout = ReadNIRS(Instrument, filename, probefile)
%
% This function reads files from any NIRS instrument used in Mesquita's lab
% and collaborators. It needs 3 input arguments:
%  1) Instrument: type of NIRS instrument used to acquire data. At the
%  moment it can be: 'CW6 (TechEn)', 'NIRScout (NIRx)', 'Imagent (ISS)'
%  2) filename: name of file(s) to be analyzed (use cell for multiple
%  files);
%  2) probefile: SD structure from a .sd/.layout/nirsInfo.mat file
%
% If no input argument (Instrument/filename/probename is passed as function 
% then a window will pop up and guide the user on how to properly load data.
% The combined data will be stored in an object of either cw- or fd-class. 
% If no output argument is passed then a .lob file will be saved on the 
% input folder containing an object called data.
%
% Created by: R. Mesquita (2018/3/14)
% Last Modified by:
%   R. Mesquita (2018/3/14): to read both NIRx and TechEn's instruments


% Asks for instrument type to load appropriate file types
% The old version based on questdialog allowed only three options.
if ~exist('Instrument','var')
    %Instrument = GetNIRSInstrument;
    display('Which system did you use:');
    display('(1) CW6 (TechEn)');
    display('(2) NIRScout (NIRx)');
    display('(3) Imagent (ISS)');
    display('(4) Brainsight'); 
    
    
    prompt = 'Choose the corresponding from 1 to 4: ';
    SystemNumber = input(prompt);
    
    %List of Systems:
    lstSystem{1} = 'CW6 (TechEn)'; lstSystem{2} = 'NIRScout (NIRx)';
    lstSystem{3} = 'Imagent (ISS)'; lstSystem{4} = 'Brainsight';
    
    Instrument = lstSystem{SystemNumber};
end


% Handle each instrument separately
switch Instrument
    case 'CW6 (TechEn)'
        % Asks for a file if none is passed as argument
        if exist('filename','var')
            data = GetTechEnData(filename);
        else
            [data, filename] = GetTechEnData;
        end
        saved_file = [filename(1:end-5) '.lob'];
        
    case 'NIRScout (NIRx)'
        % Asks for a file if none is passed as argument
        if exist('filename','var') & exist('probefile','var')
            data = GetNIRxData(filename, probefile);
        else
            [data, filename] = GetNIRxData;
        end
        saved_file = [filename(1:end-4) '.lob'];
        
    case 'Imagent (ISS)'
        [data, filename] = ReadISS2TechEn; 
        if iscell(filename)
            saved_file = [filename{1}(1:end-4) '.lob'];
        else
            saved_file = [filename(1:end-4) '.lob'];
        end    
        
    case 'Brainsight'
        [data,filename] = GetBrainsightData;
        saved_file = [filename(1:end-5) '.lob'];       
        
end


%     % Delete useless data points
%     data = DiscardDataPoints(data);

    % Check stimulation vector for run
    foo = find(data.s(:)==1);
    if ~isempty(foo)
        [data.s,data.StimTriggers] = ConvertTrigger2Stim(data.t,data.s); 
    end


if nargout > 0; varargout{1} = data;    % Save data on an argument
else    % Save data on file
    save(saved_file,'-mat','data')
    uiwait( msgbox({'We saved your data in a nicer .lob file',...
        'in the same folder where your data came from',...
        'Good luck in your analysis!'},'Success!') );
end


end

% -------------------------------------

function Instrument = GetNIRSInstrument
Instrument = questdlg('Which NIRS instrument do you want to load data from?', ...
    'NIRS Instrumentation', ...
    'CW6 (TechEn)', 'NIRScout (NIRx)', 'Imagent (ISS)', 'CW6 (TechEn)',...
    'Brainsight');
end

function [data,varargout] = GetBrainsightData

    filename = uigetfile('*.nirs','Choose the .nirs file to read');
    load(filename,'-mat');
    
    varargout{1} = filename;
    
    if exist('brainsight','var')
        % Correction 1 - Import Acquisition Frequency
        save([filename '.orig']);
        SD.f = brainsight.parameters.samplingRate;
        clear brainsight
                       
        % Correction 2 - Get stimulus to the stim vector (s)
        % Ask which auxiliary channel was used
        display(['There are ' num2str(size(aux,2)) ' axuliary channels.']);
        prompt = 'Which one did you use for the triggers? ';
        aux_chan = input(prompt);
        
        % Correction 3 - Repeated triggers due to long pedal pressing
        % We will remove triggers that are too close to each other
        display('Some triggers may be repeated.');
        prompt = 'What is the minimum expected time(s) between two triggers ';
        interval_triggers = input(prompt);
        
        
        % subtract 2 seconds to be on the safest side
        if (interval_triggers)>5
            interval_triggers = interval_triggers - 2;
        end
        
        % Convert trigger time to frames
        interval_triggers = round(SD.f*interval_triggers);
        
        % Create stim vector (s)
        lst = find(aux(:,aux_chan)<1);
        
        if ~isempty(lst)
            
            if size(lst,1)>1
                
                  distance = [nan;diff(lst)];    
                  remove_triggers = find(distance<interval_triggers);
                  lst(remove_triggers) = [];
            end
        end
                
        s = zeros(size(t));
        s(lst) = 1;
        
        clear lst aux_chan ml
        
        % Correction 3 - Fix SD.MeasList and data matrix (d)
        [~,order] = sort(SD.MeasList(:,4));
        
        d_new = d(:,order); % Correct order of data (d)
        d = d_new;
        
        ml_new = SD.MeasList(order,:);
        SD.MeasList = ml_new;
        
        clear ml_new d_new order sorted
        
        % Correction 4 - Fix Det and Source Pos
        
        % Load 2D Probe
        probeFile = uigetfile('*.probe','Load 2D Probe');
        load(probeFile,'-mat')
        % Replace 2D information
        SD.SrcPos = coor(1:SD.nSrcs,:);
        SD.DetPos = coor(SD.nSrcs+1:end,:);
        clear coor;
        
        % Convert 3D values to cm
        
        SD.DetPos3D = SD.DetPos3D/10;
        SD.SrcPos3D = SD.SrcPos3D/10;
        SD.SpatialUnit = 'cm';
        
    end
    save(filename)
    data = cw_nirs;
    data.t = t; data.d = d; data.s = s; data.aux = aux; data.SD = SD;

    
end


function [data, varargout] = GetTechEnData(filename)
if ~exist('filename','var')
    uiwait( msgbox({'Missing a file to read the data from.'...
        'Please select the .nirs file to read'}, ...
        'Warning','warn') );
    [nirsfile, filepath] = uigetfile('*.nirs','Choose the .nirs file to read');
    filename = [filepath nirsfile];
end

% Unwrap cw-nirs info from CW6 files
load(filename,'-mat')
data = cw_nirs;
data.t = t; data.d = d; data.s = s; data.aux = aux; data.SD = SD;
clear aux d m* s* S* t*

% Send filename back if not passed as input argument
if nargout > 0; varargout{1} = filename; end
end

function [data, varargout] = GetNIRxData(filename, probefile)
if ~exist('filename','var')
    uiwait( msgbox({'Missing a folder to read the data from.'...
        'Please select the FOLDER where all NIRS data is located'}, ...
        'Warning','warn') );
    DataFolder = uigetfile_n_dir; % get directory list
    if isempty(DataFolder)
        return
    else
        for i = 1:length(DataFolder) % for each directory
            NIRx_file = dir([DataFolder{i} filesep '*.hdr']); % get file list
            
            for files = 1:length(NIRx_file) % for each file
                filename = [DataFolder{i} filesep NIRx_file(files).name];
                
                % Asks for a probe file if none is passed as argument
                if ~exist('probefile','var'); probefile = GetProbeFileName; end
            end
        end
    end
end
% Unwrap cw-nirs info from NIRx files
data = RunNIRx2nirs(filename, probefile);

if nargout > 0; varargout{1} = filename; end
end




function data = RunNIRx2nirs(HDRfile,probefile)

% Part 1: Open and Extract Header Information
[SD,mask,markers] = GetHDRInfo(HDRfile);


% Part 2: Add in the optical geometry information
SD = ReadLayoutFile(probefile,SD);
ml = SD.MeasList;

% Part 3: Read and Organize Intensity from NIRx Wavelength files
[t,d] = GetIntensityData(HDRfile,mask,SD);


% Part 4: Organize stimulation & auxiliary info from markers
[s,aux] = GetStim(markers,length(t));

% Part 5: create CW-NIRS object
data = cw_nirs;
data.t = t; data.d = d; data.s = s; data.aux = aux; data.SD = SD;
clear aux d m* s* S* t*

end

function [SD,mask,markers] = GetHDRInfo(HDRfile)

fi = fopen(HDRfile,'rt');

% Locate experiment protocol details from header
position='aaaaaaa';
while sum( position(1:7)~='Sources' ) ~= 0
    clear position
    position = fgetl(fi);
    if length(position) < 7
        position='aaaaaaa';
    end
end
SD.nSrcs = str2num(position(9:end));

position = fgetl(fi);
SD.nDets = str2num(position(11:end));

WavelengthFlag=0;

while(WavelengthFlag==0)
    
    position = fgetl(fi);
    
    if length(position)>11
        if position(1:11) == 'Wavelengths'         
            
            SD.Lambda = str2num(position(14:end-1));           
            WavelengthFlag=1;
            
        end
    end
    
end

position='aaaaaaaaaaaa';
while sum( position(1:12)~='SamplingRate' ) ~= 0
    clear position
    position = fgetl(fi);
    if length(position) < 12
        position = 'aaaaaaaaaaaa';
    end
end
SD.f = str2num(position(14:end));
%SD.Lambda = [760 850];
%if Nfiles == 1
%    disp(' ')
%    disp('IMPORTANT: We are assuming your NIRx system has the following wavelengths: 760 & 850 nm')
%    disp('          (in this order!) If that is not true, please contact someone to fix this.')
%end
%disp(['Starting to work on file ' num2str(Nfiles)])

% Get Event markers from .HDR (it will be used later to build the s
% variable)
while position(1:3)~='Eve'
    clear position
    position = fgetl(fi);
    if length(position) < 3
        position='aaa';
    end
end
position = fgetl(fi);

cnt = 1;
while position(1)~='#'
    markers(cnt,:) = str2num(position);
    position = fgetl(fi);
    cnt=cnt+1;
end
if ~exist('markers')
    markers = [];
end

% Store the SD map distribution
for i=1:5
    position = fgetl(fi);
end
cnt = 1;
while position(1)~='#'
    mask(cnt,:) = str2num(position);
    position = fgetl(fi);
    cnt=cnt+1;
end
fclose(fi);

% Create the MeasList Matrix
SD.MeasList = [];
cnt = 1;
for i=1:SD.nSrcs,
    for j=1:SD.nDets,
        if mask(i,j)~=0
            SD.MeasList(cnt,1) = i;
            SD.MeasList(cnt,2) = j;
            SD.MeasList(cnt,3) = 1;
            cnt = cnt+1;
        end
    end
end

MeasListAux = SD.MeasList;
% SD.MeasList = [SD.MeasList ones(size(SD.MeasList,1),1);...
%     SD.MeasList ones(size(SD.MeasList,1),1)*2];


SD.MeasList = [MeasListAux ones(size(SD.MeasList,1),1)];

for Nlambda = 2:size(SD.Lambda,2)
   SD.MeasList = [SD.MeasList;...
        [MeasListAux Nlambda*ones(size(MeasListAux,1),1)]];        
end

end


function SD = ReadLayoutFile(probefile,SD);

% Find the SD info for Homer. The .layout file was created by our lab to
% add in structural information (by simple mapping or by digitalization
% with a digitizer).
% IF THERE IS A DIGITIZER, their x,y,z coordinates should be in a
% structure called 'RS_MRI'. The structure contains 2 variables:
% RS_MRI.Map (with 2D topographic info) and RS_MRI.Map3d (with 3D
% coordinates).
% IF THERE IS NO DIGITIZER, then it will look for a variable called Map
% (dimension: Mx3). The first s lines are the x,y,z coordinates of the
% S sources, and the last d lines are the x,y,z coordinates of the D
% detectors. Therefore, M = S + D!
%
% But now NIRx can make you export information on probe through a
% *_nirsInfo.mat file. The below works for both ways.
% (adapted by RM, 1/19/2017).
%
% (adapted by SL,  23/03/2017) - 1. Comented lines 219-223
%                                2. Added a condition to delete
%                                    the SD.SrcPOS and DetPos in the
%                                    second data.
%

% Commented after change on Sep 25, 2017 (not sure if we still need this):

%%% The variable RS_MRI.Map is usually 2d. It makes conflict.
% if NIRxDir > 1
%     SD = rmfield( SD, 'SrcPos' );
%     SD = rmfield( SD, 'SrcPos_3d' );
%     SD = rmfield(SD , 'DetPos');
%     SD = rmfield(SD , 'DetPos_3d');
% end


% There is a .layout file in the folder
if ~isempty(probefile)
    load(probefile,'-mat');
    
    if  exist('nirsInfo')% There is at least a nirsInfo.mat file in the folder
        SD.SrcPos = nirsInfo.probeInfo.probes.coords_s2;
        SD.SrcPos(:,3) = 0;
        SD.DetPos = nirsInfo.probeInfo.probes.coords_d2;
        SD.DetPos(:,3) = 0;
        SD.SrcPos_3d = nirsInfo.probeInfo.probes.coords_s3;
        SD.DetPos_3d = nirsInfo.probeInfo.probes.coords_d3;
    elseif exist('probeInfo')
        SD.SrcPos = probeInfo.probes.coords_s2;
        SD.SrcPos(:,3) = 0;
        SD.DetPos = probeInfo.probes.coords_d2;
        SD.DetPos(:,3) = 0;
        SD.SrcPos_3d = probeInfo.probes.coords_s3;
        SD.DetPos_3d = probeInfo.probes.coords_d3;
    else
        
        if isstruct(RS_MRI) % There is digitizer information
            for i=1:SD.nSrcs,
                SD.SrcPos(i,:) = RS_MRI.Map(i,:);
                SD.SrcPos_3d(i,:) = RS_MRI.Map3d(i,:);
            end
            
            for i=1:SD.nDets,
                SD.DetPos(i,:) =  RS_MRI.Map(i + (SD.nSrcs),:);
                SD.DetPos_3d(i,:) = RS_MRI.Map3d(i + (SD.nSrcs),:);
            end
            SD.SrcPos(:,3) = 0;
            SD.DetPos(:,3) = 0;
            
            
        elseif exist('Map') % work with 2D projection
            for i=1:SD.nSrcs,
                SD.SrcPos(i,:) = Map(i,:);
            end
            for i=1:SD.nDets,
                SD.DetPos(i,:) =  RS_MRI.Map(i + (SD.nSrcs),:);
            end
            SD.SrcPos(:,3) = 0;
            SD.DetPos(:,3) = 0;
            
        else
            disp(' ')
            disp('There is an error with your .layout file. Please fix it before you continue...')
        end
    end
end
SD.SpatialUnit = 'cm';


% Check if the number of sources/detectors match the number described in
% the header file. Sometimes there are extra sources/detectors that are not
% really used; one should change the SD.nSrcs/SD.nDets then.
MeasListMax = max(SD.MeasList);
if MeasListMax(1) ~= SD.nSrcs
    disp(['    ATTENTION: The system says you used ' num2str(SD.nSrcs) ' sources.'])
    disp(['               However, we noticed that only ' num2str(MeasListMax(1)) ' are connected.'])
    disp(['               We will discard the remanining ' num2str(abs(SD.nSrcs - MeasListMax(1))) ' sources... OK?'])
    SD.nSrcs = MeasListMax(1);
end
if MeasListMax(2) ~= SD.nDets
    disp(['    ATTENTION: The system says you used ' num2str(SD.nDets) ' detectors.'])
    disp(['               However, we noticed that only ' num2str(MeasListMax(2)) ' are connected.'])
    disp(['               We will discard the remanining ' num2str(abs(SD.nSrcs - MeasListMax(1))) ' detectors... OK?'])
    SD.nDets = MeasListMax(2);
end

end

function [t,d] = GetIntensityData(HDRfile,mask,SD)

% Import Data from all Wavelengths and save on d 
d=[];
for Nlambda=1:size(SD.Lambda,2)
    
    % Import data from wavelenth Nlambda
    fi = fopen([HDRfile(1:end-4) '.wl' num2str(Nlambda)],'r');
    
    % Concatenate data from wavelength Nlambda
    data{Nlambda} = fscanf(fi,'%g',[size(mask,1)*size(mask,2) inf]);
    data{Nlambda} = data{Nlambda}';
    
    % Close File
    fclose(fi);
    
    % Trash SD pairs that were not used, and keep only the ones listed on
    % MeasList.
    % We will frist assign as NAN the useless data then exclude from
    % variable d
    limpeza = reshape(mask',1,size(mask,1)*size(mask,2));
    rk = find(limpeza==0);
    limpeza(rk) = NaN;
    limpeza = ones(size(data{Nlambda},1),1)*limpeza;
    data{Nlambda} = data{Nlambda}.*limpeza;
    
    % Create matrix similar to .nirs
    d = [d data{Nlambda}];
end

% Remove NANs from the data d
lst = find(isnan(d(1,:)));
d(:,lst)=[];

% Create the time varialbe based on the acquisition frequency
t = linspace(0,length(d)/SD.f,length(d));
t = t';

end

function [s,aux] = GetStim(markers,Tpts)

% Create the stim matrix based on the event markers
if ~isempty(markers)
    s = zeros(Tpts,max(markers(:,2)));
    for i=1:max(markers(:,2))
        lst = find( markers(:,2) == i );
        if ~isempty(markers)
            s(markers(lst,3),i)=1;
        end
    end
else
    s = zeros(Tpts,4);
end

% Create auxiliary variables
aux = zeros(Tpts,4);

end

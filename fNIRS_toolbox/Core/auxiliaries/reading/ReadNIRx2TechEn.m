function varargout = ReadNIRx2TechEn(filename, probefile)

% This function helps you to convert NIRS data from NIRx systems onto
% TechEn's file format (the same used by HomER).
%
% You can call the function either with or without arguments. With the
% latter, a window will appear for you to choose the folder(s) where your
% NIRx data is(are). (NOTE: You can choose more than one folder by holding
% the 'Shift' key to start a batch conversion.) In each folder where you
% have data you will need to have either a .layout file or a *nirsInfo.mat
% file containing information from the source-detector geometry. If you
% have more than one set of data in a folder, you can use the same
% .layout/nirsInfo.mat file for all files in that folder.
%
% If you are running a batch analysis, you can pass the .hdr filename and
% the .layou/*nirsInfo.mat filename as input arguments, and the function
% will process the data for a given filename.
%
% In either case, the output of the function is a .nirs file for each NIRx
% file that will be saved with the same name as the NIRx file, on the same
% folder. An output argument will send the path to the .nirs file;
%
% Example of usage for batch command:
%
% HDRfile = '/Users/rickson/data/NIRS-2017-07-07_006.hdr';
% PROBEfile = '/Users/rickson/data/WholeHead64Channels.layout';
% ReadNIRx2TechEn(HDRfile,PROBEfile)                (OR:)
% nirs_file = ReadNIRx2TechEn(HDRfile,PROBEfile);
%
% Created by: R. Mesquita on May 12, 2012.
%
% Modified on:
%   - Jan 20, 2017: to read and analyze several NIRS files simultaneously
%   (R. Mesquita)
%   - Jan 22, 2017: to be able to pass a nirsInfo.mat file (created by NIRx
%   analysis software) as the probe info (R. Mesquita)
%   - Sep 25, 2017: to run the function with arguments so that it can go to
%   batch processing (R. Mesquita)
%   - Sep 29, 2017: to read .wl1 and .wl2 based on the variable 'mask'
%   rather than SD.nDets/nSrcs. This will make data1/data2 have 16x16=256
%   columns; columns not used will be deleted by variable 'limpeza' (R.
%   Mesquita & S. Novi)



% Asks for a file if none is passed as argument
if ~exist('filename')
    NIRxFolder = uigetfile_n_dir; % get directory list
    % Get probe layout
    [probelayout,pathprobe] = uigetfile('*.layout;*.mat','Pick the .layout/SDInfo.mat file containing SD info');
    probefile = [pathprobe probelayout];
    
    for i=1:length(NIRxFolder) % for each directory
        NIRx_file= dir([NIRxFolder{i} filesep 'NIR*.hdr']); % get file list
        
        for files=1:length(NIRx_file) % for each file
            disp(['Running file ' num2str(i) ' of ' num2str(length(NIRx_file)) ' in folder ' num2str(i)])
            RunNIRx2nirs([NIRxFolder{i} filesep NIRx_file(files).name],probefile)
        end
    end
    
else
    % Simply convert the filename passed throught the function
    RunNIRx2nirs(filename,probefile)
    varargout{1} = [filename(1:end-4) '.nirs'];
end
disp('Conversion to .nirs successfully completed. Please run analysis now.')





function RunNIRx2nirs(HDRfile,probefile)

% Part 1: Open and Extract Header Information
[SD,mask,markers] = GetHDRInfo(HDRfile);


% Part 2: Add in the optical geometry information
SD = ReadLayoutFile(probefile,SD);
ml = SD.MeasList;

% Part 3: Read and Organize Intensity from NIRx Wavelength files
[t,d] = GetIntensityData(HDRfile,mask,SD);


% Part 4: Organize stimulation & auxiliary info from markers
[s,aux] = GetStim(markers,length(t));


% Part 5: Save data onto new file
file = HDRfile;
clear ma* probefile HDRfile
save([file(1:end-4) '.nirs'],'-MAT')




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

position='aaaaaaaaaaaa';
while sum( position(1:12)~='SamplingRate' ) ~= 0
    clear position
    position = fgetl(fi);
    if length(position) < 12
        position = 'aaaaaaaaaaaa';
    end
end
SD.f = str2num(position(14:end));
SD.Lambda = [760 850];
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
SD.MeasList = [SD.MeasList ones(size(SD.MeasList,1),1);...
    SD.MeasList ones(size(SD.MeasList,1),1)*2];


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
    
    
    function [t,d] = GetIntensityData(HDRfile,mask,SD)
    
    % Import data from both wavelengths
    fi = fopen([HDRfile(1:end-4) '.wl1'],'r');
    fid = fopen([HDRfile(1:end-4) '.wl2'],'r');
    
    
    % Concatenate data from wavelength 1
    data1 = fscanf(fi,'%g',[size(mask,1)*size(mask,2) inf]);
    data1 = data1';
    
    % Concatenate data from wavelength 2
    data2 = fscanf(fid,'%g',[size(mask,1)*size(mask,2) inf]);
    data2 = data2';
    
    % Close files
    fclose(fi);
    fclose(fid);
    
    % Trash SD pairs that were not used, and keep only the ones listed on
    % MeasList.
    limpeza = reshape(mask',1,size(mask,1)*size(mask,2));
    rk = find(limpeza==0);
    limpeza(rk) = NaN;
    limpeza = ones(size(data1,1),1)*limpeza;
    data1 = data1.*limpeza;
    data2 = data2.*limpeza;
    
    % Create data matrix to HomER
    d = [data1 data2];
    lst = find(isnan(d(1,:)));
    d(:,lst)=[];
    
    % Create the time varialbe based on the acquisition frequency
    t = linspace(0,length(d)/SD.f,length(d));
    t = t';
    
    
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

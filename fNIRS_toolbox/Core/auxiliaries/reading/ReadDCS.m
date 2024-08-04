function varargout = ReadDCS(filename, probefile, Channels, InstrOrigin)
%
% This function reads files from DCS instruments created in Yodh's &
% Mesquita's labs. It needs 3 input arguments:
%  1) filename: name of .dat files to be analyzed (use cell for multiple
%  .dat files);
%  2) probefile: SD structure from a .sd/.layout/nirsInfo.mat file
%  3) Channels: 8, 16 (number of APDs in the DCS instrument)
%  4) InstrOrigin: since labs have different way of store data, the user
%  needs to say from where data come from ('Penn','Unicamp')
%
% If no filename/probename is passed as function argument a window will pop
% up and guide the user on how to properly load data.
% The combined data will be stored in an object of dcs class. If no output
% argument is passed then a .lob file will be saved on the input folder
% containing an object called data.
%
% Created by: R. Mesquita (2018/2/20)
% Last Modified by:
%   R. Mesquita (2018/2/23): to read both UNICAMP's and Penn' data, and to
%   collect info on # of channels and type of DCS via GUI if no argument is
%   passed.


progress = 0; % flag for waitbar while loading files

% Asks for a file if none is passed as argument
if ~exist('filename','var'); filename = GetDCSFileName; progress = 1; end

% Asks for how many channels the instrument has if info not provided
if ~exist('Channels','var'); Channels = GetNumberChannels; end

% Asks for a probe file if none is passed as argument
if ~exist('probefile','var'); probefile = GetProbeFileName; end

% Asks for the origin of the data if not provided
if ~exist('InstrOrigin','var'); InstrOrigin = GetLabData; end

% Read SD info
load(probefile,'-mat')

% Read .dat files
Nfiles = length(filename); % total number of files to read
data = dcs; % initialize object
if progress; h = waitbar(0,['Reading DCS files...']); end % create waitbar
for i = 1:Nfiles
    if progress; waitbar(i / Nfiles); end % update waitbar
    fi = fopen(filename{i},'rt'); % open file for reading as text
    switch InstrOrigin % read file according to lab
        case 'Unicamp'
            [time(i,:),d(i,:),taus(i,:),g2(i,:,:),s(i,:),rho(i,:)] = ...
                ReadDCSUnicamp(fi,Channels);
            t(i) = etime(time(i,1:6),time(1,1:6));
            SD.rho = rho(1,:)';
        case 'Penn'
            [time(i,:),d(i,:),taus(i,:),g2(i,:,:),s(i,:)] = ...
                ReadDCSPenn(fi,Channels);
            t(i) = etime(time(i,1:6),time(1,1:6));
    end
    st = fclose(fi);
    clear fi tmp st
end

% Populate dcs object
data.t = t'; data.time = time; data.s = s; data.aux = zeros(length(data.t),4);
data.taus = taus; data.g2 = g2; data.d = d(:,1:Channels);
data.SD = SD;

if progress; close(h); end  % close waitbar

if nargout > 0; varargout{1} = data;    % Save data on an argument
else    % Save data on file
    save([filename{1}(1:end-4) '.lob'],'-mat','data')
    uiwait( msgbox({'We saved your data in a nicer .lob file',...
        'in the same folder where your data came from',...
        'Good luck in your analysis!'},'Success!') );
end

end



function filename = GetDCSFileName
uiwait( msgbox({'Missing a folder to read the data from.'...
    'Please select the FOLDER where all DCS data is located'}, ...
    'Warning','warn') );
DataFolder = uigetfile_n_dir; % get directory list
if isempty(DataFolder)
    return
else
    for folder = 1:length(DataFolder)
        files = dir( [ DataFolder{folder} filesep '*_flow*.dat'] );
        for fl = 1:length(files)
            filename{fl} = [DataFolder{folder} filesep files(fl).name];
        end
    end
    filename = natsortfiles(filename);
end
end

function Channels = GetNumberChannels
Channels = questdlg('How many channels did the DCS instrument have?', ...
    'Instrument Details', ...
    '4','8','16', '16');
Channels = str2num(Channels);
end



function InstrOrigin = GetLabData
InstrOrigin = questdlg('From which lab does the dataset come from?', ...
    'Finding Data Source', ...
    'Unicamp','Penn', 'Unicamp');
end

function [time,d,taus,g2,s,rho] = ReadDCSUnicamp(fi,Channels)
% Getting exact time
tmp = fgets(fi);
tmp = tmp(2:end);
numtime = datevec(tmp, 'mm/dd/yyyy HH:MM:SS');
% Organize time stamps according to file
time = [numtime(3) numtime(1) numtime(2) numtime(4) numtime(5) numtime(6)];
% Get Intensity from APDs
while isempty(strfind(fgets(fi), '%Duration'))
end
clear tmp
tmp = fgets(fi);
d = str2num(tmp);
% Get rho given in instrument
while isempty(strfind(fgets(fi), '%Rho'))
end
clear tmp
tmp = fgets(fi);
rho = str2num(tmp);
% Get temporal autocorrelation values
while isempty(strfind(fgets(fi), '%Delay'))
end
clear tmp
tmp = fscanf(fi,'%g',[Channels+1 inf]);
taus = tmp(1,:)';
g2 = tmp(2:end,:)';
% Get marks in the frame
while isempty(strfind(fgets(fi), '%Marks'))
end
clear tmp
tmp = fgets(fi);
mark = str2num(tmp);
s = mark(1);
clear tmp mark

end

function [time,d,taus,g2,s] = ReadDCSPenn(fi,Channels);
% Getting exact time
tmp = fgets(fi);
tmp = tmp(2:end);
numtime = datevec(tmp, 'mm/dd/yyyy HH:MM:SS PM');
% Organize time stamps according to file
time = [numtime(3) numtime(1) numtime(2) numtime(4) numtime(5) numtime(6)];
% Get Intensity from APDs
while isempty(strfind(fgets(fi), '%Duration'))
end
clear tmp
tmp = fgets(fi);
d = str2num(tmp);
% Get temporal autocorrelation values
while isempty(strfind(fgets(fi), '%Delay'))
end
clear tmp
tmp = fscanf(fi,'%g',[Channels+1 inf]);
taus = tmp(1,1:end-1)';
g2 = tmp(2:end,1:end-1)';
% Get marks in the frame
s = tmp(1,end);
end
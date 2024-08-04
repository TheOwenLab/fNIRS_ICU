    % Correction 1 - Import Acquisition Frequency
    DataTemp.SD.f = DataTemp.brainsight.parameters.samplingRate;
    clear brainsight
    
    % Correction 2 - Get stimulus to the stim vector (s)
    % Ask which auxiliary channel was used
    display(['There are ' num2str(size(DataTemp.aux,2)) ' axuliary channels.']);
    prompt = 'Which one did you use for the triggers? ';
    aux_chan = input(prompt);
    
    % Create stim vector (s)
    lst = find(DataTemp.aux(:,aux_chan)<1);
    DataTemp.s = zeros(size(DataTemp.t));
    DataTemp.s(lst) = 1;
    clear lst aux_chan ml
    
    % Correction 3 - Fix SD.MeasList and data matrix (d)
    [~,order] = sort(DataTemp.SD.MeasList(:,4));
    
    d_new = DataTemp.d(:,order); % Correct order of data (d)
    DataTemp.d = d_new;
    
    ml_new = DataTemp.SD.MeasList(order,:);
    DataTemp.SD.MeasList = ml_new;
    
    clear ml_new d_new order sorted
    
    % Correction 4 - Fix Det and Source Pos
    
    % Load 2D Probe
    probeFile = uigetfile('*.probe','Load 2D Probe');
    load(probeFile,'-mat')
    % Replace 2D information
    DataTemp.SD.SrcPos = coor(1:DataTemp.SD.nSrcs,:);
    DataTemp.SD.DetPos = coor(DataTemp.SD.nSrcs+1:end,:);
    clear coor;
    
    % Convert 3D values to cm
    
    DataTemp.SD.DetPos3D = DataTemp.SD.DetPos3D/10;
    DataTemp.SD.SrcPos3D = DataTemp.SD.SrcPos3D/10;
    DataTemp.SD.SpatialUnit = 'cm';
       

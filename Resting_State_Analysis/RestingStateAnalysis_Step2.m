clear

% Load .stats data
load Regression_Physiology_6_minutes.stats -mat

AvailableSubjects = [1 9 10 11 12 13 14 15 16 17 ...
    21 24 27 28 29 30 31 33 35 37 38 39 40]; 

for Nsub=AvailableSubjects
    
    CorrMtx_SC_Regression{Nsub} = ...
        SC_Regression{Nsub}{1}.C_regressed;
    
end

save('CorrMtx_SC_6_minutes','CorrMtx_SC_Regression','BadChan')

% Threshold prior to group averages
% -1 indicates that there is no threshold
Threshold = -1;

% Channels to be analyzed. We must remove the SC
SSlist = [8 29 52 66 75 92 112 125];

cnt=1;
for Nsub = AvailableSubjects
    Channels_remove{cnt} = (unique([SSlist BadChan{Nsub}']));
    cnt=cnt+1;
end

% Convert CorrMtrix to another format
Corr_SC_Phys = [];

% Define the hemoglobin to be analyzed
% HbO: 1, HbR: 2, HbT: 3
Hb = 3;

for Nsub = AvailableSubjects
    Corr_SC_Phys = cat(3,Corr_SC_Phys,CorrMtx_SC_Regression{Nsub}(:,:,Hb));
end

% Exclude Bad Channels before group averaging them
for Nsub = 1:size(AvailableSubjects,2)
    Corr_SC_Phys(Channels_remove{Nsub},:,Nsub)=nan;
    Corr_SC_Phys(:,Channels_remove{Nsub},Nsub)=nan;
end

clear BadChan cnt SSlist

for Nsub = 1:size(AvailableSubjects,2)
    % Threshold each correlation matrix based on the
    % threshold of each volunteer (not used)
    
    % For the SC_Phys
    clear aux
    aux = Corr_SC_Phys(:,:,Nsub);
    
    lst_SC_Phys = find(aux < Threshold);
    aux(lst_SC_Phys) = nan;
    Corr_SC_Phys(:,:,Nsub) = aux;
    
    clear lst_SC_Phys aux
end

clear lst* aux*

clear CorrMtx_SC_Phys_Regression

% Tranform Corr to Z-space
Corr_SC_Phys_z = atanh(Corr_SC_Phys);

save('Corr_SC_Phys_z','Corr_SC_Phys_z')



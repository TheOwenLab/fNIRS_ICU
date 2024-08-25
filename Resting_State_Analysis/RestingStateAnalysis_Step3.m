close all
clear 
clc

load('Corr_SC_Phys_z.mat') 
Corr_group=Corr_SC_Phys_z; 
clear Corr_SC_Phys_z

% uncomment the network of interest 

%seed = [30 94]; 
%Name = 'Motor'; 

%seed = [3 50 72 111]; 
%Name = 'DMN'; 

%seed = [10 79 50 110]; 
%Name = 'FPC'; 

%seed = [44 116]; 
%Name = 'Aud';

%% average across seeds 

Hb=3; % HbT
TempMotor=Corr_group(seed,:,:);
TempMotor=nanmean(TempMotor,3);

for i = 1:length(seed)
    TempMotor(i,seed(i))=nan; 
end

seed_SC = nanmean(TempMotor,1);

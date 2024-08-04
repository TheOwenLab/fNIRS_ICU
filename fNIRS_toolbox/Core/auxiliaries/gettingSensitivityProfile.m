function[intensity] = gettingSensitivityProfile(fwMC,Imap)

% This function use the density profile from the MonteCarlo profile
% To create a distribution from the vector Imap. 
% Input: 1 - fwMC: MonteCarlo structure of each probe;
%        2 - Imap: Vector with the values to be ploted of each channel;
% Created by S. L. Novi

%%% Converting the activation map to the same profile as the MC.
Adot = zeros(size(fwMC.Adot));
for iCH=1:size(Adot,1)
    Adot(iCH,:) = fwMC.Adot(iCH,:)./max(fwMC.Adot(iCH,:));
    Adot(iCH,:) = Adot(iCH,:).*Imap(iCH);
end

%%% Intensity Vector
intensity = zeros(1,size(fwMC.Adot,2));

for Nvox = 1:size(fwMC.Adot,2)    
        [value indice] = max(abs(Adot(:,Nvox))); % Find higher absolut Value
        % Save the maximum value in the intensity vector
        intensity(1,Nvox) = Adot(indice,Nvox);
end
intensity = intensity';
end
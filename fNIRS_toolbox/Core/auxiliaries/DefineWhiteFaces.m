function wFaces = DefineWhiteFaces(fwMC,intensity,lambda)

lst = find(abs(intensity)<lambda*max(abs(intensity)));
intensity(lst) = nan; 
wFaces = [];

for Nvox = 1:size(lst,1)
    for j=1:3
         aux = find(fwMC.mesh.faces(:,j)==lst(Nvox));
         wFaces = [fwMC.mesh.faces(aux,:);wFaces];
    end
end

end
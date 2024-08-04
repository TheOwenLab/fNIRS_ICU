function[f] = CreateHrf_GLM_analysis
% This function creates a hemodynamic response function


th = 1.1; %'largura' da gama, em segundos de acordo com a literatura
f=@(x)(1/(th.*120)).*((x/th).^5).*exp(-1*(x)./th).*heaviside(x);

end
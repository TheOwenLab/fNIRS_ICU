function [X] = CreateDeasingMatrix(r,s,duration)
% Code for creating the design matrix
% r is the cw_nirs object
% s is a vector with the triggers. s here has to be the triggers for only
%one type of stimulus.
% duration can be a vector or a scalar. If it is scalar, we assume that every
% trial had the same stimulus duration. If it is a vector then it should
% match the number of trials

% Check and correct the dimension of duration
if size(duration,1)>size(duration,2)
    duration=duration';
end

%%% Create the canonical hrf
f = CreateHrf_GLM_analysis;
 
% Create  box cars based on each situmulus duration
caixa = s;
lst=find(s==1);

StimVecDuration = ones(1,size(lst,1));

if size(duration,2)==1    
    StimVecDuration(1,:)=duration;
else
    StimVecDuration = duration.*StimVecDuration;    
end

for N=1:length(lst)
    caixa(lst(N):lst(N)+round(r.SD.f*StimVecDuration(N)))=1;
end

% Convolve the hrf with the box cars
X = conv(f(r.t),caixa);
X = X(1:length(r.t));

%Normalize for the peack to be 1.
X = X./max(X);

end
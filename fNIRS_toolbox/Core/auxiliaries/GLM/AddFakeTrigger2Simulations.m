function[s] = AddFakeTrigger2Simulations(obj,tfirst,tlast,Rest)
    % Add fake triggers in order to perform simulations
    % Input AddFakeTrigger2Simulations(obj,tfirst,tlast,Rest):
    %
    %   obj - cw_nirs Object
    %   tfirst - First trial in seconds. 
    %   tlast - Maximum time from the end of the data the
    %   last trial in seconds.
    %   Rest - Vector with the maximum and minimum period of rest [min max]
    %
    % Create Vector of Stimulus
    %
    % Sergio Novi, sergiolnovi@gmail.com
    f = obj.SD.f;
    s = zeros(size(obj.d,1),1);
    
    trial = tfirst; % First Trial
    s(round(trial*f))=1;    
    
    while trial <= (size(obj.d,1)/f)-tlast
        trial = randi([(trial+Rest(1)) (trial+Rest(2))]);
        s(round(trial*f))=1;
    end
    

end
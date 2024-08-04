function[s,duration] = createStimVectorForRestin(s,f,TaskDuration)


% Take original triggers
Task_triggers =  find(s==1);

% Compute Resting triggers
Rest_triggers = Task_triggers + round(TaskDuration*f);

% Create Stim vector from Resting 
s_aux_new = zeros(size(s,1),1);
s_aux_new(Rest_triggers)=1;

% Combine Rest and Task periods
s = [s,s_aux_new];

   
duration = ...
    (Task_triggers(2:end) - Rest_triggers(1:end-1))/f;    

duration(end+1) = mean(duration);
end
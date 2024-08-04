function  [s,stim,varargout] = ConvertTrigger2Stim(t,s,varargin)

% This function checks the 's' vector in a .nirs file so that the user can
% modify it accordingly for functional processing. It reads the .nirs file
% as argument, plots the stim vector for each type of trigger, and asks the
% user which triggers are actually stimulation onset. It will then remake
% the 's' vector and save a new .nirs file (with the same filename).
%
%
% Created by: R. Mesquita on May 12, 2012.
%
% Modified on:


if nargin > 2 % vector with trials to be removed was passed
    for trigger = 1:nargin-1
        lst = find( s(:,trigger) ~= 0 );
        s(:,trigger) = zeros(length(t),1);
        s(lst(varargin{1}),trigger) = 1;
    end
    
else % user will interactively tell which trials should be removed
    close
    for trigger = 1:size(s,2)        
        h = figure;
        plot(t,s(:,trigger))
        title(['Trigger # ' num2str(trigger)])
        xlabel('Time (s)')
        ylabel('Trigger')
        lst = find( s(:,trigger) ~= 0 );
        for i=1:length(lst)
            text(t(lst(i))+1,s(lst(i),trigger)-.1,num2str(i))
        end
        stim = input('Which trigger(s) should we accept as stimulation onset? (e.g., [1 3]) ');
        s(:,trigger) = zeros(length(t),1);
        s(lst(stim),trigger) = 1;
        hold on, plot(t,s(:,trigger),'k','LineWidth',2)
    end
    hold off
end

end